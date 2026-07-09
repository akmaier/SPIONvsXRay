"""Reconstruction gallery for the dashboard: the ED phantom (7 iron inserts on a
ring + bone rod) reconstructed for EID and PCD, noise-free and noisy, at a tight
iron window so the low-contrast inserts are visible.
"""
from __future__ import annotations
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import conrad_backend
import conrad_ct
import conrad_phantom
import conrad_project as cp
import run_detectability as rf


def noise_free_sino(acc, det):
    if det == "EID":
        return -np.log(np.clip(acc["S_det"], 1e-6, None) / acc["S_air"])
    w = rf._pcd_weights(acc)
    M = np.zeros_like(acc["S_det"]); Ma = 0.0
    for b, (cd, ca) in enumerate(zip(acc["C_det"], acc["C_air"])):
        M += w[b] * cd; Ma += w[b] * ca
    return -np.log(np.clip(M, 1e-6, None) / max(Ma, 1e-6))


def main():
    # bone-free display phantom: the 7 iron inserts show cleanly without the
    # bone-rod streaks that swamp them at an iron-level window (see the paper's
    # discussion). The quantitative factorial always keeps the bone source.
    scene, inserts = conrad_phantom.build_phantom(with_bone=False)
    geo = conrad_ct.fan_geometry(n_pix=512)
    base, geo = cp.project_base_materials(inserts, geo)
    acc = rf.polychromatic_accumulators(base)
    sp = geo["voxel_mm"]

    recons = {}
    for det in ("EID", "PCD"):
        recons[(det, "noise-free")] = conrad_ct.fbp(noise_free_sino(acc, det), geo, bh_correction=False)
        recons[(det, "noisy")] = conrad_ct.fbp(rf.line_integral(acc, det, 3), geo, bh_correction=False)

    # per-DETECTOR HU normalization (EID and PCD reconstruct at different effective
    # mu, so each is referenced to its own central water background)
    N = recons[("EID", "noise-free")].shape[0]
    yy, xx = np.mgrid[0:N, 0:N]
    c = N / 2.0
    bg = ((xx - c) ** 2 + (yy - c) ** 2) < (8 / sp) ** 2

    def to_hu(img):
        mu_w = float(np.median(img[bg]))
        return 1000.0 * (img - mu_w) / mu_w

    WIN = 20  # HU: iron window; no bone rod in the display phantom
    outdir = str(conrad_backend.REPO_ROOT / "results" / "recon")
    os.makedirs(outdir, exist_ok=True)

    fig, ax = plt.subplots(2, 2, figsize=(9, 9))
    order = [("EID", "noise-free"), ("PCD", "noise-free"), ("EID", "noisy"), ("PCD", "noisy")]
    for a, key in zip(ax.ravel(), order):
        a.imshow(to_hu(recons[key]), cmap="gray", vmin=-WIN, vmax=WIN)
        a.set_title(f"{key[0]} — {key[1]}", fontsize=13)
        a.axis("off")
    fig.suptitle(f"7 iron inserts (c$_{{Fe}}$ = 0–1.09 mg/ml) on a ring — bone omitted for display\n"
                 f"window [−{WIN}, {WIN}] HU  ·  70 000 photons/pixel  ·  beam-hardening off",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(f"{outdir}/recon_gallery.png", dpi=120, bbox_inches="tight"); plt.close(fig)
    print(f"wrote {outdir}/recon_gallery.png")

    # labeled single panel (noise-free EID) with insert concentrations
    fig, a = plt.subplots(figsize=(6.0, 6.0))
    a.imshow(to_hu(recons[("EID", "noise-free")]), cmap="gray", vmin=-WIN, vmax=WIN)
    for ins in inserts:
        cx, cy = ins["center_mm"]
        col, row = cx / sp + c, cy / sp + c
        if ins.get("c_fe") is not None:
            a.text(col, row + 14, f"{ins['c_fe']:.2f}", color="#ffd23b", fontsize=8,
                   ha="center", va="center")
    a.set_title("Iron inserts labeled by c$_{Fe}$ [mg/ml] (EID, noise-free)"); a.axis("off")
    fig.tight_layout(); fig.savefig(f"{outdir}/recon_labeled.png", dpi=120, bbox_inches="tight"); plt.close(fig)
    print(f"wrote {outdir}/recon_labeled.png")


if __name__ == "__main__":
    main()
