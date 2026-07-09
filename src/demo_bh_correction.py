"""Beam-hardening correction comparison for the photon-counting detector:
no correction vs. single water precorrection (on the combined signal) vs.
bin-dependent water precorrection in the sinogram domain.

Key point (per-bin, sinogram-domain): each energy bin is narrower and hardens
differently, so we calibrate a water precorrection polynomial for *each bin's*
spectrum, apply it to that bin's line integral, then weighted-sum the corrected
bin sinograms and reconstruct ONCE. Correcting in the sinogram domain and adding
before FBP costs a single reconstruction (not one per bin), because after the
per-bin nonlinear correction the combination and FBP are both linear.
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
import materials
import run_detectability as rf


def _bin_poly(E, s, mu_w, lo, hi):
    w_eff = np.where((E >= lo) & (E < hi), s, 0.0)
    return conrad_ct.water_precorrection_poly(E, w_eff, mu_w)


def main():
    scene, inserts = conrad_phantom.build_phantom(with_bone=True)   # realistic BH source
    geo = conrad_ct.fan_geometry(n_pix=512)
    base, geo = cp.project_base_materials(scene, geo)
    acc = rf.polychromatic_accumulators(base)
    E, s, edges = acc["E"], acc["s"], acc["edges"]
    mu_w = materials.linear_attenuation("water", E)
    w = rf._pcd_weights(acc)
    nb = len(edges) - 1

    # noise-free per-bin line integrals
    p_bin = [-np.log(np.clip(acc["C_det"][b], 1e-6, None) / acc["C_air"][b]) for b in range(nb)]

    # per-bin water precorrection polynomials (each bin's own spectrum)
    bin_polys = [_bin_poly(E, s, mu_w, edges[b], edges[b + 1]) for b in range(nb)]

    # All three use the SAME sinogram-domain combination p = sum_b w_b * p_b
    # (sum of per-bin log line integrals). Calibrate the single-correction poly for
    # THAT combined water response (not the count-domain one) so the only difference
    # between "single" and "per-bin" is where the correction is applied.
    L = np.linspace(0.0, 30.0, 400)
    p_bin_L = []
    for b in range(nb):
        m = (E >= edges[b]) & (E < edges[b + 1])
        sb = s[m] / s[m].sum()
        p_bin_L.append(np.array([-np.log(np.sum(sb * np.exp(-mu_w[m] * l))) for l in L]))
    p_comb_L = sum(w[b] * p_bin_L[b] for b in range(nb))
    mu_ref = float((p_comb_L[1] - p_comb_L[0]) / (L[1] - L[0]))       # combined low-dose slope
    poly_comb = np.polyfit(p_comb_L, mu_ref * L, 4)

    p_none = sum(w[b] * p_bin[b] for b in range(nb))                 # combine, no correction
    p_single = np.polyval(poly_comb, p_none)                         # one correction on combined
    p_perbin = sum(w[b] * np.polyval(bin_polys[b], p_bin[b]) for b in range(nb))  # per-bin, summed

    recons = {"no correction": conrad_ct.fbp(p_none, geo),
              "single water corr.": conrad_ct.fbp(p_single, geo),
              "per-bin corr. (sinogram)": conrad_ct.fbp(p_perbin, geo)}

    # HU vs the central water background; cupping = center/edge in a water annulus
    N = next(iter(recons.values())).shape[0]
    sp = geo["voxel_mm"]
    yy, xx = np.mgrid[0:N, 0:N]; c = N / 2.0
    r = np.sqrt((xx - c) ** 2 + (yy - c) ** 2)
    core = r < (6 / sp)
    edge_ring = (r > (60 / sp)) & (r < (74 / sp))   # water near body edge (body r=80mm)

    def to_hu(img):
        return 1000.0 * (img - np.median(img[core])) / np.median(img[core])

    outdir = str(conrad_backend.REPO_ROOT / "results" / "recon")
    os.makedirs(outdir, exist_ok=True)
    fig, ax = plt.subplots(2, 3, figsize=(13, 8.4))
    for j, (name, img) in enumerate(recons.items()):
        hu = to_hu(img)
        cup = float(np.median(img[core]) / np.median(img[edge_ring]))
        ax[0, j].imshow(hu, cmap="gray", vmin=-30, vmax=30)
        ax[0, j].set_title(f"{name}\ncenter/edge = {cup:.3f}", fontsize=11); ax[0, j].axis("off")
        # radial profile of HU
        prof = np.array([hu[(r >= k) & (r < k + 4)].mean() for k in range(0, 160, 4)])
        ax[1, j].plot(np.arange(len(prof)) * 4 * sp, prof, lw=1.4)
        ax[1, j].axhline(0, color="#888", lw=0.8); ax[1, j].set_ylim(-40, 40)
        ax[1, j].set_xlabel("radius [mm]"); ax[1, j].set_ylabel("HU"); ax[1, j].grid(alpha=0.3)
        ax[1, j].set_title("radial HU profile")
    fig.suptitle("PCD beam-hardening correction: none vs. single water vs. per-bin (sinogram domain)\n"
                 "window [-30,30] HU; flat profile = no cupping", fontsize=12)
    fig.tight_layout()
    fig.savefig(f"{outdir}/bh_correction_compare.png", dpi=120, bbox_inches="tight"); plt.close(fig)
    for name, img in recons.items():
        print(f"  {name:28s} center/edge = {float(np.median(img[core])/np.median(img[edge_ring])):.4f}")
    print("wrote", f"{outdir}/bh_correction_compare.png")


if __name__ == "__main__":
    main()
