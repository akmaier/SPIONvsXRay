"""Forward-projection method comparison for the SPION fan-beam pipeline.

Compares three base-material forward projectors on the bone-free phantom and the
current SPION geometry (fan_geometry(n_pix=512), deltaT=1), NOISE-FREE EID:

  (1) hard     : hard 0/1 indicator-grid rasterisation projected with CONRAD's
                 FanBeamProjector2D.projectRayDrivenCL (the current/broken method,
                 reproduced inline here exactly as committed in b4a2f28).
  (2) aa       : anti-aliased fractional-coverage grid (conrad_project method="aa").
  (3) analytic : exact closed-form disk-chord Radon (conrad_project method="analytic").

Deliverables (written to the scratchpad dir passed as OUTDIR, or a default):
  - hu_regression_fix.png : iron dHU vs c_Fe, three curves + ideal line
  - recon_forward_methods.png : noise-free EID recon (3 forwards side by side)
  - printed per-insert dHU table + edge-rim overshoot table

Run:  .venv/bin/python scripts/forward_method_compare.py
"""
from __future__ import annotations
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import conrad_backend
import conrad_ct
import conrad_phantom
import conrad_project as cp

CG = conrad_backend.class_getter

OUTDIR = os.environ.get(
    "OUTDIR",
    "/private/tmp/claude-501/-Users-maier-Documents-SPIONvsXRay/"
    "f77329d1-6870-4df8-b0b7-dd7a968fa1bb/scratchpad",
)


def project_hard(inserts, geo):
    """HARD 0/1 indicator-grid forward projection (reproduces committed b4a2f28).

    Rasterise each material as a hard 0/1 mask on the detector-matched grid
    (1 px = deltaT mm, pixel index convention x=(i-n/2)*vox matching the projector
    origin), inserts overriding the body, and project with projectRayDrivenCL.
    """
    conrad_backend.setup()
    n = int(round(geo["maxT"] / geo["deltaT"]))
    vox = geo["deltaT"]
    yy, xx = np.mgrid[0:n, 0:n]
    x = (xx - n / 2.0) * vox
    y = (yy - n / 2.0) * vox
    assigned = np.zeros((n, n), bool)
    masks = {}
    for ins in inserts:                                  # inserts override body
        cx, cy = ins["center_mm"]; r = ins["radius_mm"]; nm = ins["name"]
        m = ((x - cx) ** 2 + (y - cy) ** 2 <= r * r) & ~assigned
        masks.setdefault(nm, np.zeros((n, n), np.float32))[m] = 1.0
        assigned |= m
    body = ((x * x + y * y) <= conrad_phantom.BODY_RADIUS_MM ** 2) & ~assigned
    masks.setdefault("water", np.zeros((n, n), np.float32))[body] = 1.0
    FP = CG("edu.stanford.rsl.tutorial.fan").FanBeamProjector2D
    fp = FP(geo["focal"], geo["maxBeta"], geo["deltaBeta"], geo["maxT"], geo["deltaT"])
    base = {}
    for nm, mask in masks.items():
        g = conrad_ct.np_to_grid2d(mask)
        g.setSpacing(vox, vox)
        g.setOrigin(-(n * vox) / 2.0, -(n * vox) / 2.0)
        base[nm] = conrad_ct.grid2d_to_np(fp.projectRayDrivenCL(g)) * vox
    return base


def recon_for(base, geo):
    det = cp.detector_sinograms(base, add_noise=False)
    return conrad_ct.fbp(det["eid"], geo)


def edge_rim_pct(recon, geo, body_r_mm=conrad_phantom.BODY_RADIUS_MM):
    """Body-edge overshoot %: ring just inside r=body_r vs interior mean.

    Interior mean = water region well away from inserts and the edge (r < 30 mm,
    which for the bone-free phantom is pure water between the insert circle at
    50 mm and the centre). Rim = annulus [body_r-4, body_r-1] mm. Overshoot% =
    100*(rim_mean - interior_mean)/interior_mean.
    """
    N = recon.shape[0]
    sp = geo.get("voxel_mm", 1.0)
    yy, xx = np.mgrid[0:N, 0:N]
    r = np.hypot(xx - N / 2.0, yy - N / 2.0) * sp
    interior = (r < 30.0)
    rim = (r >= body_r_mm - 4.0) & (r <= body_r_mm - 1.0)
    mu_int = float(recon[interior].mean())
    mu_rim = float(recon[rim].mean())
    return 100.0 * (mu_rim - mu_int) / mu_int, mu_int, mu_rim


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(OUTDIR, exist_ok=True)
    conrad_backend.setup()

    scene, inserts = conrad_phantom.build_phantom(with_bone=False)
    geo = conrad_ct.fan_geometry(n_pix=512)
    print(f"geometry: deltaT={geo['deltaT']} maxT={geo['maxT']} "
          f"n_det={int(round(geo['maxT']/geo['deltaT']))} imgN={geo['imgN']} "
          f"voxel_mm={geo['voxel_mm']}")

    # --- three forwards ---
    base_hard = project_hard(inserts, geo)
    base_aa, _ = cp.project_base_materials(inserts, geo, method="aa")
    base_ana, _ = cp.project_base_materials(inserts, geo, method="analytic")

    recon_hard = recon_for(base_hard, geo)
    recon_aa = recon_for(base_aa, geo)
    recon_ana = recon_for(base_ana, geo)

    meas_hard = cp.measure_inserts(recon_hard, geo, inserts)
    meas_aa = cp.measure_inserts(recon_aa, geo, inserts)
    meas_ana = cp.measure_inserts(recon_ana, geo, inserts)

    # SPION inserts only (bone-free phantom has none, but guard anyway), sorted by c_Fe
    def spions(meas):
        return [m for m in meas if m.get("c_form") is not None]
    mh, ma, mn = spions(meas_hard), spions(meas_aa), spions(meas_ana)
    order = np.argsort([m["c_fe"] for m in mn])
    mh = [mh[i] for i in order]; ma = [ma[i] for i in order]; mn = [mn[i] for i in order]
    c_fe = np.array([m["c_fe"] for m in mn])
    hu_hard = np.array([m["iron_delta_hu"] for m in mh])
    hu_aa = np.array([m["iron_delta_hu"] for m in ma])
    hu_ana = np.array([m["iron_delta_hu"] for m in mn])

    # --- table ---
    print("\nPer-insert IRON dHU (local-annulus, c0-corrected; noise-free EID):")
    hdr = f"{'c_Fe[mg/ml]':>12} {'dHU_hard':>10} {'dHU_aa':>10} {'dHU_analytic':>13}"
    print(hdr); print("-" * len(hdr))
    table_lines = [hdr, "-" * len(hdr)]
    for cf, hh, aa, an in zip(c_fe, hu_hard, hu_aa, hu_ana):
        line = f"{cf:12.4f} {hh:10.2f} {aa:10.2f} {an:13.2f}"
        print(line); table_lines.append(line)

    # --- monotonicity check ---
    def mono(v):
        return bool(np.all(np.diff(v) > 0))
    print(f"\nmonotonic (strictly increasing) : hard={mono(hu_hard)}  "
          f"aa={mono(hu_aa)}  analytic={mono(hu_ana)}")

    # --- edge-rim ---
    rim_hard = edge_rim_pct(recon_hard, geo)
    rim_aa = edge_rim_pct(recon_aa, geo)
    rim_ana = edge_rim_pct(recon_ana, geo)
    print("\nBody-edge overshoot %% (rim [76,79] mm vs interior r<30 mm, noise-free EID):")
    print(f"  hard     : {rim_hard[0]:+.2f} %%  (interior mu={rim_hard[1]:.5f}, rim mu={rim_hard[2]:.5f})")
    print(f"  aa       : {rim_aa[0]:+.2f} %%  (interior mu={rim_aa[1]:.5f}, rim mu={rim_aa[2]:.5f})")
    print(f"  analytic : {rim_ana[0]:+.2f} %%  (interior mu={rim_ana[1]:.5f}, rim mu={rim_ana[2]:.5f})")

    # --- line plot ---
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    # ideal line: through the origin with the analytic high-c slope (fit top 3 points)
    hi = c_fe >= c_fe[len(c_fe) // 2]
    slope = np.polyfit(c_fe[hi], hu_ana[hi], 1)[0]
    cx = np.linspace(0, c_fe.max() * 1.02, 50)
    ax.plot(cx, slope * cx, "k--", lw=1.2, alpha=0.7,
            label=f"ideal linear (analytic high-c slope {slope:.1f} HU per mg/ml)")
    ax.plot(c_fe, hu_hard, "o-", color="#d62728", label="hard 0/1 grid (current)")
    ax.plot(c_fe, hu_aa, "s-", color="#1f77b4", label="anti-aliased grid (aa)")
    ax.plot(c_fe, hu_ana, "^-", color="#2ca02c", label="analytic disk-chord")
    ax.axhline(0, color="grey", lw=0.6)
    ax.set_xscale("symlog", linthresh=0.05)
    ax.set_xticks(c_fe)
    ax.set_xticklabels([f"{v:.3g}" for v in c_fe], rotation=45, fontsize=8)
    ax.set_xlabel("iron concentration c_Fe [mg Fe/ml]")
    ax.set_ylabel("iron dHU [HU]  (c0-corrected, local annulus)")
    ax.set_title("Iron contrast regression: forward-projection method\n"
                 "(bone-free phantom, noise-free EID, fan geometry deltaT=1)")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    p1 = os.path.join(OUTDIR, "hu_regression_fix.png")
    fig.savefig(p1, dpi=140); plt.close(fig)

    # --- recon images ---
    stack = np.concatenate([recon_hard, recon_aa, recon_ana])
    lo, hi_w = np.percentile(stack, [20, 99.7])
    fig, ax = plt.subplots(1, 3, figsize=(15, 5.2))
    for a, rc, ttl in zip(ax, (recon_hard, recon_aa, recon_ana),
                          ("hard 0/1 grid (current)", "anti-aliased (aa)", "analytic")):
        im = a.imshow(rc, cmap="gray", vmin=lo, vmax=hi_w)
        a.set_title(ttl); a.axis("off")
    fig.suptitle(f"Noise-free EID reconstruction, bone-free phantom "
                 f"(window [{lo:.4f},{hi_w:.4f}])", y=0.99)
    fig.tight_layout()
    p2 = os.path.join(OUTDIR, "recon_forward_methods.png")
    fig.savefig(p2, dpi=140); plt.close(fig)

    print(f"\n[ok] wrote {p1}")
    print(f"[ok] wrote {p2}")
    return dict(c_fe=c_fe, hu_hard=hu_hard, hu_aa=hu_aa, hu_ana=hu_ana,
                rim=(rim_hard, rim_aa, rim_ana), table="\n".join(table_lines))


if __name__ == "__main__":
    main()
