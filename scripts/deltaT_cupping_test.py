"""Sub-mm detector cupping test for CONRAD's fan backprojector normalization.

Reconstructs a 360-degree fan FBP of two CONRAD phantoms while varying ONLY the
detector element size deltaT in {1.0, 0.5, 0.25} mm (recon grid + object fixed).
CONRAD's fan backprojector omits the detector integration measure deltaT from its
normalization (normalizationFactor = maxBetaIndex/PI, both CPU and CL). At
deltaT != 1 mm this produces radial cupping and a wrong DC scale. The proposed fix
divides the normalization by deltaT.

This script runs ONE phase per process ("before" or "after") because a rebuilt Java
class cannot be reloaded in a live JVM -- the caller rebuilds conrad_ext between
phases and re-invokes with a fresh interpreter. Results (recon arrays + metrics) are
written to an .npz per phase; a final "compare" phase builds side-by-side figures.

Chain per view (RamLakKernel -- dimensionally correct 1/deltaS, NOT SheppLogan):
    FanBeamProjector2D.projectRayDrivenCL -> CosineFilter
    -> RamLakKernel((int)(maxT/deltaT), deltaT)
    -> FanBeamBackprojector2D.backprojectPixelDrivenCL

Usage:
    python scripts/deltaT_cupping_test.py before  <out_dir>
    python scripts/deltaT_cupping_test.py after   <out_dir>
    python scripts/deltaT_cupping_test.py compare <out_dir>
"""
from __future__ import annotations
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

# --- fixed geometry (only deltaT varies) ---
FOCAL = 750.0          # source->isocenter (CONRAD focalLength) [mm]
MAXT = 256.0           # detector length [mm], FIXED across deltaT
N_PIX = 256            # recon grid
VOXEL_MM = 1.0         # recon voxel [mm]
N_VIEWS = 500
DELTA_TS = [1.0, 0.5, 0.25]
PHANTOMS = ["UniformCircle", "MickeyMouse"]


def _geo(deltaT):
    n_det = int(round(MAXT / deltaT))
    gamma_max = float(np.arctan((MAXT / 2.0 - 0.5 * deltaT) / FOCAL))
    maxBeta = 2.0 * np.pi                      # full 360-degree scan
    deltaBeta = maxBeta / N_VIEWS
    return dict(focal=FOCAL, maxBeta=maxBeta, deltaBeta=deltaBeta,
                maxT=MAXT, deltaT=deltaT, n_det=n_det,
                imgN=N_PIX, spacing=VOXEL_MM, voxel_mm=VOXEL_MM,
                gamma_max=gamma_max)


def _make_phantom(name):
    """Return a CONRAD phantom as a numpy (N,N) array."""
    import conrad_ct as cc
    if name == "UniformCircle":
        P = cc._cls("edu.stanford.rsl.tutorial.phantoms", "UniformCircleGrid2D")
        g = P(N_PIX, N_PIX)
    elif name == "MickeyMouse":
        P = cc._cls("edu.stanford.rsl.tutorial.phantoms", "MickeyMouseGrid2D")
        g = P(N_PIX, N_PIX)
    elif name == "SheppLogan":
        P = cc._cls("edu.stanford.rsl.tutorial.phantoms", "SheppLogan")
        g = P(N_PIX)
    else:
        raise ValueError(name)
    return cc.grid2d_to_np(g)


def _recon(phantom_np, geo):
    """Project + FBP with RamLakKernel using the CONRAD API (GPU projector+backproj)."""
    import conrad_ct as cc
    focal, maxT, deltaT = geo["focal"], geo["maxT"], geo["deltaT"]
    maxBeta, deltaBeta = geo["maxBeta"], geo["deltaBeta"]

    # forward projection (GPU ray-driven)
    g = cc.np_to_grid2d(phantom_np)
    g.setSpacing(geo["spacing"], geo["spacing"])
    g.setOrigin(-(geo["imgN"] - 1) / 2.0 * geo["spacing"],
                -(geo["imgN"] - 1) / 2.0 * geo["spacing"])
    FP = cc._cls("edu.stanford.rsl.tutorial.fan", "FanBeamProjector2D")
    fp = FP(focal, maxBeta, deltaBeta, maxT, deltaT)
    try:
        sino_np = cc.grid2d_to_np(fp.projectRayDrivenCL(g))
    except Exception:
        sino_np = cc.grid2d_to_np(fp.projectRayDriven(g))

    sino = cc.np_to_grid2d(sino_np)
    sino.setSpacing(deltaT, deltaBeta)

    # cosine + RamLak ramp (dimensionally correct 1/deltaS), per view
    Cos = cc._cls("edu.stanford.rsl.tutorial.fan", "CosineFilter")
    RL = cc._cls("edu.stanford.rsl.tutorial.filters", "RamLakKernel")
    cos = Cos(focal, maxT, deltaT)
    ram = RL(int(round(maxT / deltaT)), deltaT)
    n_views = int(sino.getSize()[1])
    for th in range(n_views):
        cos.applyToGrid(sino.getSubGrid(th))
    for th in range(n_views):
        ram.applyToGrid(sino.getSubGrid(th))

    # distance-weighted pixel-driven backprojection (GPU, spacing-aware)
    BP = cc._cls("edu.stanford.rsl.tutorial.fan", "FanBeamBackprojector2D")
    bp = BP(focal, deltaT, deltaBeta, int(geo["imgN"]), int(geo["imgN"]))
    try:
        bp.setSpacing(float(geo["voxel_mm"]))
        recon = cc.grid2d_to_np(bp.backprojectPixelDrivenCL(sino))
    except Exception:
        recon = cc.grid2d_to_np(bp.backprojectPixelDriven(sino))
    return recon


def _metrics(recon, phantom):
    """DC (mean over object interior), center/edge ratio, min/max.

    Center = small central disk; edge annulus = near the object boundary. The
    object radius is read from the phantom mask (>0.5*max) so both phantoms work;
    for the uniform circle this is the clean cupping read.
    """
    N = recon.shape[0]
    yy, xx = np.mgrid[0:N, 0:N]
    cx = cy = (N - 1) / 2.0
    r = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2)

    mask = phantom > 0.5 * phantom.max()
    if mask.sum() == 0:
        mask = phantom > 0
    # object radius = max radius of masked pixels (robust to Mickey ears via 95th pct)
    obj_r = float(np.percentile(r[mask], 98))

    center_disk = r < (0.15 * obj_r)
    edge_annulus = (r > 0.80 * obj_r) & (r < 0.92 * obj_r)
    interior = r < (0.85 * obj_r)

    dc = float(recon[interior].mean())
    c_val = float(recon[center_disk].mean())
    e_val = float(recon[edge_annulus].mean())
    ratio = c_val / e_val if e_val != 0 else float("nan")
    return dict(dc=dc, center=c_val, edge=e_val, center_edge=ratio,
                vmin=float(recon.min()), vmax=float(recon.max()),
                obj_r=obj_r)


def run_phase(phase, out_dir):
    import conrad_backend
    conrad_backend.setup()
    os.makedirs(out_dir, exist_ok=True)
    results = {}
    phantoms_np = {name: _make_phantom(name) for name in PHANTOMS}

    for name in PHANTOMS:
        for dt in DELTA_TS:
            geo = _geo(dt)
            recon = _recon(phantoms_np[name], geo)
            m = _metrics(recon, phantoms_np[name])
            key = f"{name}_{dt}"
            results[key + "__recon"] = recon.astype(np.float32)
            results[key + "__metrics"] = np.array(
                [m["dc"], m["center"], m["edge"], m["center_edge"],
                 m["vmin"], m["vmax"], m["obj_r"]], dtype=np.float64)
            print(f"[{phase}] {name:14s} dt={dt:<5} "
                  f"DC={m['dc']:.5f} c/e={m['center_edge']:.4f} "
                  f"min={m['vmin']:.4f} max={m['vmax']:.4f}")

    for name in PHANTOMS:
        results[name + "__phantom"] = phantoms_np[name].astype(np.float32)

    np.savez_compressed(os.path.join(out_dir, f"{phase}.npz"), **results)
    print(f"[{phase}] wrote {os.path.join(out_dir, phase + '.npz')}")


METRIC_NAMES = ["dc", "center", "edge", "center_edge", "vmin", "vmax", "obj_r"]


def _load(out_dir, phase):
    return dict(np.load(os.path.join(out_dir, f"{phase}.npz")))


def compare(out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    before = _load(out_dir, "before")
    after = _load(out_dir, "after")

    def metrics(d, name, dt):
        v = d[f"{name}_{dt}__metrics"]
        return {k: v[i] for i, k in enumerate(METRIC_NAMES)}

    # --- per-phantom BEFORE|AFTER grid ---
    fig_paths = []
    for name in PHANTOMS:
        # common window/level per phantom from the AFTER dt=1 interior (a stable ref)
        ref = after[f"{name}_1.0__recon"]
        obj_r = metrics(after, name, 1.0)["obj_r"]
        N = ref.shape[0]
        yy, xx = np.mgrid[0:N, 0:N]
        r = np.sqrt((xx - (N - 1) / 2.0) ** 2 + (yy - (N - 1) / 2.0) ** 2)
        interior = r < 0.85 * obj_r
        mu = float(ref[interior].mean())
        vlo, vhi = 0.0, mu * 1.8

        fig, axes = plt.subplots(len(DELTA_TS), 2, figsize=(7.5, 3.3 * len(DELTA_TS)))
        for i, dt in enumerate(DELTA_TS):
            for j, (phase, d) in enumerate([("BEFORE", before), ("AFTER", after)]):
                rec = d[f"{name}_{dt}__recon"]
                m = metrics(d, name, dt)
                ax = axes[i, j]
                ax.imshow(rec, cmap="gray", vmin=vlo, vmax=vhi)
                ax.set_title(f"{phase}  dt={dt} mm\n"
                             f"c/e={m['center_edge']:.3f}  DC={m['dc']:.4f}",
                             fontsize=10)
                ax.axis("off")
        fig.suptitle(f"{name}  (fixed object & recon grid; window shared per phantom)",
                     fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        p = os.path.join(out_dir, f"compare_{name}.png")
        fig.savefig(p, dpi=130)
        plt.close(fig)
        fig_paths.append(p)
        print("wrote", p)

    # --- radial profile figure (uniform circle: the clean cupping read) ---
    if "UniformCircle" in PHANTOMS:
        name = "UniformCircle"
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
        for ax, (phase, d) in zip(axes, [("BEFORE", before), ("AFTER", after)]):
            for dt in DELTA_TS:
                rec = d[f"{name}_{dt}__recon"]
                N = rec.shape[0]
                prof = rec[N // 2, :]
                ax.plot(np.arange(N) - N / 2.0, prof, label=f"dt={dt}")
            ax.set_title(f"{name}  horizontal profile  ({phase})")
            ax.set_xlabel("x [px]"); ax.set_ylabel("mu")
            ax.grid(alpha=0.3); ax.legend()
        fig.tight_layout()
        p = os.path.join(out_dir, "compare_UniformCircle_profiles.png")
        fig.savefig(p, dpi=130)
        plt.close(fig)
        fig_paths.append(p)
        print("wrote", p)

    # --- summary table ---
    print("\n" + "=" * 92)
    hdr = (f"{'phantom':14s} {'dt':>5s} | {'c/e BEFORE':>11s} {'c/e AFTER':>10s} "
           f"| {'DC BEFORE':>10s} {'DC AFTER':>9s}")
    print(hdr)
    print("-" * 92)
    rows = []
    for name in PHANTOMS:
        for dt in DELTA_TS:
            mb, ma = metrics(before, name, dt), metrics(after, name, dt)
            print(f"{name:14s} {dt:>5} | {mb['center_edge']:>11.4f} "
                  f"{ma['center_edge']:>10.4f} | {mb['dc']:>10.5f} {ma['dc']:>9.5f}")
            rows.append((name, dt, mb, ma))
    print("=" * 92)
    return fig_paths, rows


if __name__ == "__main__":
    phase = sys.argv[1] if len(sys.argv) > 1 else "before"
    out_dir = sys.argv[2] if len(sys.argv) > 2 else "."
    if phase == "compare":
        compare(out_dir)
    else:
        run_phase(phase, out_dir)
