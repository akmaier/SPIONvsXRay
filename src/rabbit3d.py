"""3D RabbitCT proof-of-concept: cone-beam forward projection with the real
extracted geometry, and the projected contrast of an inserted SPION tumor.

Uses the 496 projection matrices from src/rabbitct.py. A ray-driven digital
radiograph (DRR) is formed by, for each detector pixel, casting the ray from the
decoded source through the pixel and trilinearly sampling the rabbit volume along
it (Beer-Lambert line integral). Inserting an 8 cm^3 iron tumor and differencing
the DRRs shows the tumor's projected signal. A full polychromatic cone-beam FDK
detectability study is the documented next step; this validates the geometry and
the 3D SPION phantom end to end.
"""
from __future__ import annotations
import numpy as np
from scipy.ndimage import map_coordinates

import rabbitct
import conrad_backend

ROOT = conrad_backend.REPO_ROOT
RCTD = str(ROOT / "data" / "rabbitct" / "rabbitct_512-v2.rctd")
VOLF = str(ROOT / "data" / "rabbitct" / "reference_256.vol")


def source_and_dirs(P, us, vs, det_ds=1):
    """Decode source position and per-pixel ray directions (world mm) from a
    projection matrix P (world[mm]->detector[px])."""
    Pp = P.copy()
    Pp[:2, :] = Pp[:2, :] / det_ds
    M = Pp[:, :3]
    Minv = np.linalg.inv(M)
    C = -Minv @ Pp[:, 3]                                  # source (world mm)
    uu, vv = np.meshgrid(us, vs)
    pix = np.stack([uu.ravel(), vv.ravel(), np.ones(uu.size)], 0)  # 3xN
    dirs = Minv @ pix                                    # 3xN world directions
    dirs = dirs / np.linalg.norm(dirs, axis=0, keepdims=True)
    return C, dirs, uu.shape


def drr(vol, P, R_vox, det_shape, det_ds=4, n_steps=384, half_mm=310.0):
    """Ray-driven cone-beam DRR (line integral) of `vol` for projection matrix P."""
    L = vol.shape[0]
    Sy = det_shape[1] // det_ds
    Sx = det_shape[0] // det_ds
    us = np.arange(Sx); vs = np.arange(Sy)
    C, dirs, shp = source_and_dirs(P, us, vs, det_ds)
    # march t across the volume-containing sphere around isocenter
    d0 = np.linalg.norm(C)
    ts = np.linspace(d0 - half_mm, d0 + half_mm, n_steps)
    dt = ts[1] - ts[0]
    acc = np.zeros(dirs.shape[1], np.float32)
    off = (L - 1) / 2.0
    for t in ts:
        pts = C[:, None] + t * dirs                      # 3xN world (x,y,z)
        ix = pts[0] / R_vox + off
        iy = pts[1] / R_vox + off
        iz = pts[2] / R_vox + off
        acc += map_coordinates(vol, [iz, iy, ix], order=1, mode="constant", cval=0.0)
    return (acc * dt * 1e-3).reshape(shp)                # arb. line-integral units


def insert_tumor(vol, R_vox, center_idx, radius_mm, delta_raw):
    """Add an iron tumor (raw-unit bump) inside soft tissue."""
    L = vol.shape[0]
    zz, yy, xx = np.mgrid[0:L, 0:L, 0:L]
    r = radius_mm / R_vox
    m = ((zz - center_idx[0])**2 + (yy - center_idx[1])**2 + (xx - center_idx[2])**2) < r**2
    m &= (vol > 300) & (vol < 1300)                      # keep it in soft tissue
    out = vol.copy().astype(np.float32)
    out[m] += delta_raw
    return out, int(m.sum())


def main():
    import os, matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    mats, hd = rabbitct.read_matrices(RCTD)
    vol = rabbitct.load_reference_volume(VOLF, 256)
    R = 512 * hd["R_L"] / 256.0                          # 256^3 voxel size (mm)
    det = (hd["S_x"], hd["S_y"])

    # tumor in the neck soft tissue; exaggerate x40 for a visible PoC contrast
    tvol, nvox = insert_tumor(vol, R, (150, 150, 128), 12.4, 40 * 3.4)
    outdir = str(ROOT / "results" / "rabbit"); os.makedirs(outdir, exist_ok=True)

    fig, ax = plt.subplots(2, 3, figsize=(12, 7))
    for j, k in enumerate([0, 124, 248]):
        base = drr(vol, mats[k], R, det)
        tum = drr(tvol, mats[k], R, det)
        lo, hi = np.percentile(base, [2, 99])
        ax[0, j].imshow(base, cmap="gray", vmin=lo, vmax=hi)
        ax[0, j].set_title(f"DRR view {k}"); ax[0, j].axis("off")
        d = tum - base
        mx = np.abs(d).max() + 1e-9
        ax[1, j].imshow(d, cmap="magma", vmin=0, vmax=mx)
        ax[1, j].set_title(f"tumor projected signal (view {k})"); ax[1, j].axis("off")
    fig.suptitle(f"RabbitCT cone-beam DRR (real geometry) + SPION tumor "
                 f"({nvox} voxels, 8 cm$^3$)", fontsize=12)
    fig.tight_layout(); fig.savefig(f"{outdir}/rabbit_drr_tumor.png", dpi=120); plt.close(fig)
    print(f"[rabbit3d] tumor {nvox} voxels; wrote {outdir}/rabbit_drr_tumor.png")


if __name__ == "__main__":
    main()
