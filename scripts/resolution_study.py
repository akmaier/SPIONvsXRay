"""Characterization study: fan-beam FBP recon of a UNIFORM DISC vs (recon voxel
size) x (detector element size deltaT), with the PHYSICAL object and PHYSICAL
detector held fixed.

CHARACTERIZATION ONLY. Modifies no source code, applies no fixes. Uses the CONRAD
API exactly as src/conrad_ct.py does (`_cls`, `np_to_grid2d`, `grid2d_to_np`) and
CONRAD's own CosineFilter / Ram-Lak / Shepp-Logan / FanBeamBackprojector2D.

Forward projection is ANALYTIC (exact chord through the disc) so the physical
object and detector are truly fixed regardless of recon grid or detector sampling.

Outputs (PNGs) go to the scratchpad dir set in OUT_DIR below.
"""
from __future__ import annotations
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Make src importable (conrad_ct, conrad_backend live there).
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "src"))

import conrad_backend  # noqa: E402
from conrad_ct import _cls, np_to_grid2d, grid2d_to_np  # noqa: E402

OUT_DIR = ("/private/tmp/claude-501/-Users-maier-Documents-SPIONvsXRay/"
           "f77329d1-6870-4df8-b0b7-dd7a968fa1bb/scratchpad")
os.makedirs(OUT_DIR, exist_ok=True)

# ----------------------------- fixed physical setup -----------------------------
FOCAL = 750.0          # SID [mm]
MAXT = 256.0           # physical detector length [mm]  (FIXED)
R_DISC = 80.0          # disc radius [mm]                (FIXED)
DISC_VAL = 1.0         # disc value                      (FIXED)
N_VIEWS = 720
DELTA_BETA = 2.0 * np.pi / N_VIEWS
FOV_MM = 256.0

VOXELS = [1.0, 0.5, 0.25]        # rows
DELTA_TS = [1.0, 0.5, 0.25]      # cols
KERNELS = {
    "ramlak": "RamLakKernel",
    "shepplogan": "SheppLoganKernel",
}


# ----------------------------- analytic forward projection ----------------------
def analytic_sinogram(deltaT: float) -> np.ndarray:
    """Exact fan-beam sinogram [N_VIEWS, n_det] of a uniform disc (center 0, R_DISC).

    For view beta the source is a=(focal cos b, focal sin b); the detector point for
    bin center t is P=(t sin b, -t cos b). The line integral of the disc along ray
    a->P is the chord 2*sqrt(R^2 - d^2), d = perpendicular distance of disc center 0
    to the ray (0 when d>=R). Bin centers: t=(i+0.5)*deltaT - maxT/2.
    """
    n_det = int(round(MAXT / deltaT))
    i = np.arange(n_det)
    t = (i + 0.5) * deltaT - MAXT / 2.0            # [n_det] bin centers [mm]
    betas = np.arange(N_VIEWS) * DELTA_BETA        # [N_VIEWS]
    cb, sb = np.cos(betas), np.sin(betas)          # [N_VIEWS]

    # source a and detector point P per (view, bin)
    ax = (FOCAL * cb)[:, None]                      # [V,1]
    ay = (FOCAL * sb)[:, None]
    px = t[None, :] * sb[:, None]                   # [V,n_det]
    py = -t[None, :] * cb[:, None]

    dx = px - ax
    dy = py - ay
    seglen = np.hypot(dx, dy)
    # perpendicular distance from origin (disc center) to the line through a,P:
    # |cross(P-a, a)| / |P-a|   (cross of 2D vectors -> scalar)
    cross = np.abs(dx * ay - dy * ax)
    d = cross / seglen
    inside = d < R_DISC
    chord = np.zeros_like(d)
    chord[inside] = 2.0 * np.sqrt(R_DISC**2 - d[inside] ** 2)
    return (DISC_VAL * chord).astype(np.float32)


# ----------------------------- CONRAD FBP chain --------------------------------
def reconstruct(sino_np: np.ndarray, voxel: float, deltaT: float,
                kernel_name: str) -> np.ndarray:
    """Fan-beam FBP via CONRAD classes, per the study spec."""
    n_pix = int(round(FOV_MM / voxel))
    n_det = int(round(MAXT / deltaT))

    sino = np_to_grid2d(sino_np)                 # Grid2D(width=n_det, height=N_VIEWS)
    sino.setSpacing(deltaT, DELTA_BETA)

    Cos = _cls("edu.stanford.rsl.tutorial.fan", "CosineFilter")
    Kern = _cls("edu.stanford.rsl.tutorial.filters", kernel_name)
    cos = Cos(FOCAL, MAXT, deltaT)
    kern = Kern(int(round(MAXT / deltaT)), deltaT)

    n_views = int(sino.getSize()[1])
    for th in range(n_views):
        cos.applyToGrid(sino.getSubGrid(th))
    for th in range(n_views):
        kern.applyToGrid(sino.getSubGrid(th))

    BP = _cls("edu.stanford.rsl.tutorial.fan", "FanBeamBackprojector2D")
    bp = BP(FOCAL, deltaT, DELTA_BETA, n_pix, n_pix)
    used_gpu = False
    if conrad_backend.opencl_available():
        try:
            bp.setSpacing(float(voxel))
            rec = grid2d_to_np(bp.backprojectPixelDrivenCL(sino))
            used_gpu = True
        except Exception as e:
            print(f"  [GPU FAILED n_pix={n_pix} n_det={n_det}] {e!r} -> CPU fallback")
            rec = grid2d_to_np(bp.backprojectPixelDriven(sino))
    else:
        rec = grid2d_to_np(bp.backprojectPixelDriven(sino))
    return rec, used_gpu


# ----------------------------- ROI measurements --------------------------------
def radius_grid(n_pix: int, voxel: float) -> np.ndarray:
    """Per-pixel physical radius [mm] from the recon center."""
    c = (n_pix - 1) / 2.0
    yy, xx = np.mgrid[0:n_pix, 0:n_pix]
    return np.hypot(xx - c, yy - c) * voxel


def measure(rec: np.ndarray, voxel: float) -> dict:
    n_pix = rec.shape[0]
    r = radius_grid(n_pix, voxel)
    interior = r < 0.6 * R_DISC
    dc = float(rec[interior].mean())
    center = r < 0.12 * R_DISC
    edge = (r > 0.75 * R_DISC) & (r < 0.85 * R_DISC)
    ce = float(rec[center].mean() / rec[edge].mean())
    # horizontal profile through center row, x-axis in physical mm
    mid = n_pix // 2
    prof = rec[mid, :].astype(float)
    xmm = (np.arange(n_pix) - (n_pix - 1) / 2.0) * voxel
    return dict(dc=dc, ce=ce, prof=prof, xmm=xmm, n_pix=n_pix)


# ----------------------------- run the grid ------------------------------------
def main():
    conrad_backend.setup()
    # cache analytic sinograms per deltaT (independent of voxel / kernel)
    sino_cache = {dt: analytic_sinogram(dt) for dt in DELTA_TS}

    results = {k: {} for k in KERNELS}   # results[kern][(voxel,deltaT)] = dict
    for kern_key, kern_name in KERNELS.items():
        print(f"=== kernel {kern_key} ({kern_name}) ===")
        for voxel in VOXELS:
            for dt in DELTA_TS:
                sino = sino_cache[dt]
                rec, gpu = reconstruct(sino, voxel, dt, kern_name)
                m = measure(rec, voxel)
                m["rec"] = rec
                m["gpu"] = gpu
                results[kern_key][(voxel, dt)] = m
                print(f"  voxel={voxel:>4}  dt={dt:>4}  n_pix={m['n_pix']:>4}  "
                      f"DC={m['dc']:.4f}  c/e={m['ce']:.4f}  gpu={gpu}")

    make_image_grids(results)
    make_line_plots(results)
    make_profile_plots(results)
    print_tables(results)
    return results


# ----------------------------- plotting ----------------------------------------
def make_image_grids(results):
    for kern_key in KERNELS:
        res = results[kern_key]
        # shared grayscale window across the whole grid
        vmin = min(r["rec"].min() for r in res.values())
        vmax = max(r["rec"].max() for r in res.values())
        fig, axes = plt.subplots(3, 3, figsize=(12, 12))
        for ri, voxel in enumerate(VOXELS):
            for ci, dt in enumerate(DELTA_TS):
                ax = axes[ri, ci]
                m = res[(voxel, dt)]
                ax.imshow(m["rec"], cmap="gray", vmin=vmin, vmax=vmax)
                ax.set_title(f"vox={voxel} dt={dt}\nDC={m['dc']:.3f} c/e={m['ce']:.3f}",
                             fontsize=10)
                ax.axis("off")
        fig.suptitle(f"Recon grid ({kern_key})  rows=voxel {VOXELS}  cols=deltaT "
                     f"{DELTA_TS}\nshared window [{vmin:.3f},{vmax:.3f}]", fontsize=13)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        path = os.path.join(OUT_DIR, f"res_grid_{kern_key}.png")
        fig.savefig(path, dpi=110)
        plt.close(fig)
        print("wrote", path)


def make_line_plots(results):
    for kern_key in KERNELS:
        res = results[kern_key]
        fig, ax = plt.subplots(1, 3, figsize=(18, 5))
        # (a) DC vs deltaT, one line per voxel
        for voxel in VOXELS:
            ys = [res[(voxel, dt)]["dc"] for dt in DELTA_TS]
            ax[0].plot(DELTA_TS, ys, "o-", label=f"voxel={voxel}")
        ax[0].set_xlabel("deltaT [mm]"); ax[0].set_ylabel("DC (interior mean)")
        ax[0].set_title(f"DC vs deltaT ({kern_key})"); ax[0].legend()
        ax[0].invert_xaxis(); ax[0].grid(alpha=0.3)
        # (b) center/edge vs deltaT, one line per voxel
        for voxel in VOXELS:
            ys = [res[(voxel, dt)]["ce"] for dt in DELTA_TS]
            ax[1].plot(DELTA_TS, ys, "o-", label=f"voxel={voxel}")
        ax[1].set_xlabel("deltaT [mm]"); ax[1].set_ylabel("center/edge ratio")
        ax[1].set_title(f"center/edge vs deltaT ({kern_key})"); ax[1].legend()
        ax[1].invert_xaxis(); ax[1].grid(alpha=0.3)
        # (c) DC vs voxel, one line per deltaT
        for dt in DELTA_TS:
            ys = [res[(voxel, dt)]["dc"] for voxel in VOXELS]
            ax[2].plot(VOXELS, ys, "s-", label=f"deltaT={dt}")
        ax[2].set_xlabel("voxel [mm]"); ax[2].set_ylabel("DC (interior mean)")
        ax[2].set_title(f"DC vs voxel ({kern_key})"); ax[2].legend()
        ax[2].invert_xaxis(); ax[2].grid(alpha=0.3)
        fig.tight_layout()
        path = os.path.join(OUT_DIR, f"res_lines_{kern_key}.png")
        fig.savefig(path, dpi=120)
        plt.close(fig)
        print("wrote", path)

    # kernel comparison: DC vs deltaT and c/e vs deltaT at voxel=0.5
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    for kern_key in KERNELS:
        res = results[kern_key]
        dc = [res[(0.5, dt)]["dc"] for dt in DELTA_TS]
        ce = [res[(0.5, dt)]["ce"] for dt in DELTA_TS]
        ax[0].plot(DELTA_TS, dc, "o-", label=kern_key)
        ax[1].plot(DELTA_TS, ce, "o-", label=kern_key)
    ax[0].set_title("DC vs deltaT @ voxel=0.5 (kernel compare)")
    ax[0].set_xlabel("deltaT [mm]"); ax[0].set_ylabel("DC")
    ax[1].set_title("center/edge vs deltaT @ voxel=0.5 (kernel compare)")
    ax[1].set_xlabel("deltaT [mm]"); ax[1].set_ylabel("center/edge")
    for a in ax:
        a.invert_xaxis(); a.legend(); a.grid(alpha=0.3)
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "res_lines_kernel_compare.png")
    fig.savefig(path, dpi=120)
    plt.close(fig)
    print("wrote", path)


def make_profile_plots(results):
    for kern_key in KERNELS:
        res = results[kern_key]
        fig, ax = plt.subplots(1, 2, figsize=(14, 5))
        # (i) 3 deltaT at fixed voxel=0.5
        for dt in DELTA_TS:
            m = res[(0.5, dt)]
            ax[0].plot(m["xmm"], m["prof"], label=f"deltaT={dt}")
        ax[0].axvline(-R_DISC, color="k", ls=":", lw=0.8)
        ax[0].axvline(R_DISC, color="k", ls=":", lw=0.8)
        ax[0].axhline(1.0, color="gray", ls="--", lw=0.8)
        ax[0].set_title(f"Center profile, voxel=0.5, vary deltaT ({kern_key})")
        ax[0].set_xlabel("x [mm]"); ax[0].set_ylabel("recon value"); ax[0].legend()
        ax[0].grid(alpha=0.3)
        # (ii) 3 voxel at fixed deltaT=0.5
        for voxel in VOXELS:
            m = res[(voxel, 0.5)]
            ax[1].plot(m["xmm"], m["prof"], label=f"voxel={voxel}")
        ax[1].axvline(-R_DISC, color="k", ls=":", lw=0.8)
        ax[1].axvline(R_DISC, color="k", ls=":", lw=0.8)
        ax[1].axhline(1.0, color="gray", ls="--", lw=0.8)
        ax[1].set_title(f"Center profile, deltaT=0.5, vary voxel ({kern_key})")
        ax[1].set_xlabel("x [mm]"); ax[1].set_ylabel("recon value"); ax[1].legend()
        ax[1].grid(alpha=0.3)
        fig.tight_layout()
        path = os.path.join(OUT_DIR, f"res_profiles_{kern_key}.png")
        fig.savefig(path, dpi=120)
        plt.close(fig)
        print("wrote", path)


def print_tables(results):
    for kern_key in KERNELS:
        res = results[kern_key]
        print(f"\n===== {kern_key} : DC (rows=voxel, cols=deltaT {DELTA_TS}) =====")
        header = "voxel\\dt " + " ".join(f"{dt:>9}" for dt in DELTA_TS)
        print(header)
        for voxel in VOXELS:
            row = " ".join(f"{res[(voxel,dt)]['dc']:>9.4f}" for dt in DELTA_TS)
            print(f"{voxel:>7} {row}")
        print(f"----- {kern_key} : center/edge (rows=voxel, cols=deltaT) -----")
        print(header)
        for voxel in VOXELS:
            row = " ".join(f"{res[(voxel,dt)]['ce']:>9.4f}" for dt in DELTA_TS)
            print(f"{voxel:>7} {row}")


if __name__ == "__main__":
    main()
