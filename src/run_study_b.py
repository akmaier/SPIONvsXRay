"""Study B — heterogeneous (vascular) tumor uptake vs. the homogeneous model.

SPEC 5.9: the same delivered iron mass is confined to 150 um vessels occupying
10% of the tumor volume (=> 10x local concentration). Because the CT voxel
(0.5 mm) cannot resolve 150 um vessels, each voxel partial-volume-averages vessel
+ tissue. With mass conserved the mean tumor iron per voxel is unchanged, so to
first order Study B == Study A. This script quantifies the second-order effects the
homogeneous model misses, all against the real 90 kVp spectrum and CONRAD magnetite:

  (A) Partial volume / resolution: ROI-mean iron contrast is invariant (mass
      conservation); the per-voxel peak only approaches the resolved 10x contrast
      as the voxel shrinks toward the 150 um vessel size.
  (B) Beam-hardening nonlinearity: concentrating iron into vessels makes the ray
      integral of a convex exp differ from the homogeneous case (Jensen gap). We
      compute the resulting line-integral / HU error.
  (C) Structural noise: the per-voxel vessel count is binomial, adding texture;
      we give its magnitude per voxel and after ROI averaging.

Outputs results/study_b/{study_b.json, fig_study_b.pdf}.
"""
from __future__ import annotations
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import conrad_backend
import materials
import spectrum as spec
from config import tumor_iron_conc, PHANTOM, SPECTRUM

VESSEL_UM = 150.0
VESSEL_FRAC = 0.10
LOCAL_MULT = 1.0 / VESSEL_FRAC          # 10x local concentration
C_FE_MEAN = tumor_iron_conc(10.0)       # realistic 6 mg dose, mg Fe/ml
L_TUMOR_CM = 2.0 * PHANTOM.tumor_radius_mm / 10.0
L_BODY_CM = PHANTOM.body_diameter_mm / 10.0
N0 = SPECTRUM.photons_per_pixel


def partial_volume(voxel_mm_list, rng):
    """Fine two-compartment tumor texture -> ROI-mean and per-voxel stats vs voxel size.

    Fine grid at 30 um; 150 um vessel blocks, 10% set to 10x mean (mass conserved).
    Values are in units of the homogeneous mean iron concentration (so 1.0 == Study A).
    """
    fine_um = 30.0
    block = int(round(VESSEL_UM / fine_um))          # vessel block = 5 fine px
    nb = 200                                          # 200x200 blocks ~ 30 mm patch
    blocks = (rng.random((nb, nb)) < VESSEL_FRAC).astype(np.float32) * LOCAL_MULT
    fine = np.kron(blocks, np.ones((block, block), np.float32))   # relative iron conc
    out = {}
    for v in voxel_mm_list:
        step = max(1, int(round(v * 1000.0 / fine_um)))
        nx = (fine.shape[0] // step) * step
        f = fine[:nx, :nx]
        vox = f.reshape(f.shape[0] // step, step, f.shape[1] // step, step).mean(axis=(1, 3))
        out[v] = dict(roi_mean=float(vox.mean()),
                      voxel_std=float(vox.std()),
                      voxel_peak=float(vox.max()))
    return out


def beam_hardening_gap():
    """Line-integral error from concentrating iron into vessels (Jensen/BH gap).

    Homogeneous ray: uniform c_Fe over the tumor chord. Vessel ray: the same mean
    iron, but the chord alternates between 10x-iron vessel segments (fraction 0.10)
    and iron-free tissue. Both have identical iron optical-depth MEAN; the
    polychromatic EID line integral differs only through the curvature of
    -log(sum flux e^{-tau}). We report the HU-equivalent difference at the realistic
    dose and at 5x dose.
    """
    E, flux, _ = spec.standard_spectrum()
    s = flux / flux.sum()
    mu_t = materials.linear_attenuation("water", E)                 # 1/cm
    ox = materials.oxide_contrast_massatten(E)                      # per g Fe

    def eid_lineintegral(c_profile_gcm3):
        # c_profile: list of (length_cm, c_Fe_gcm3) tumor segments; body tissue added
        tau = mu_t * L_BODY_CM
        for length_cm, c in c_profile_gcm3:
            tau = tau + (mu_t * length_cm + c * ox * length_cm)     # tissue+iron in segment
        I = np.sum(N0 * s * np.exp(-tau) * E)
        Iair = np.sum(N0 * s * E)
        return -np.log(I / Iair)

    def mu_water_eff():
        # effective linear att for HU scaling (spectrum-weighted, 1/cm)
        return np.sum(s * mu_t)

    res = {}
    for scale, tag in ((1.0, "realistic_6mg"), (5.0, "5x_dose")):
        c_mean = 1e-3 * C_FE_MEAN * scale                          # g/cm^3
        # homogeneous: whole chord at c_mean
        p_hom = eid_lineintegral([(L_TUMOR_CM, c_mean)])
        # vessel: 10% of chord at 10x, 90% at 0 (mass conserved)
        p_ves = eid_lineintegral([(VESSEL_FRAC * L_TUMOR_CM, c_mean * LOCAL_MULT),
                                  ((1 - VESSEL_FRAC) * L_TUMOR_CM, 0.0)])
        # convert line-integral difference to an approximate HU error over the tumor
        d_hu = 1000.0 * (p_ves - p_hom) / (mu_water_eff() * L_TUMOR_CM)
        res[tag] = dict(p_hom=float(p_hom), p_ves=float(p_ves),
                        lineint_gap=float(p_ves - p_hom), hu_gap=float(d_hu))
    return res


def structural_noise(voxel_mm=0.5):
    """Per-voxel binomial vessel-count texture, relative to the homogeneous mean."""
    n_vessels_per_voxel = (voxel_mm * 1000.0 / VESSEL_UM) ** 2      # 2D areal count
    # per-voxel iron (relative to mean) has mean 1, std from Binomial(n, p)*10/n
    rel_std = LOCAL_MULT * np.sqrt(n_vessels_per_voxel * VESSEL_FRAC * (1 - VESSEL_FRAC)) \
        / n_vessels_per_voxel
    return dict(n_vessels_per_voxel=float(n_vessels_per_voxel),
                per_voxel_rel_std=float(rel_std))


def main():
    rng = np.random.default_rng(0)
    voxels = [0.15, 0.25, 0.39, 0.5, 1.0]
    pv = partial_volume(voxels, rng)
    bh = beam_hardening_gap()
    sn = structural_noise(0.5)

    outdir = str(conrad_backend.REPO_ROOT / "results" / "study_b")
    os.makedirs(outdir, exist_ok=True)
    result = dict(partial_volume=pv, beam_hardening=bh, structural_noise=sn,
                  c_fe_mean=C_FE_MEAN, vessel_um=VESSEL_UM, vessel_frac=VESSEL_FRAC)
    with open(os.path.join(outdir, "study_b.json"), "w") as f:
        json.dump(result, f, indent=2)

    # figure: ROI-mean (flat, mass conserved) vs per-voxel peak/std vs voxel size
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    vs = np.array(voxels)
    ax.plot(vs * 1000, [pv[v]["roi_mean"] for v in voxels], "o-", color="#2c3e50",
            label="ROI-mean (mass conserved)")
    ax.plot(vs * 1000, [pv[v]["voxel_peak"] for v in voxels], "s--", color="#c0392b",
            label="per-voxel peak")
    ax.plot(vs * 1000, [pv[v]["voxel_std"] for v in voxels], "^:", color="#e67e22",
            label="per-voxel std (structural)")
    ax.axvline(VESSEL_UM, color="#888", ls=":", lw=1)
    ax.text(VESSEL_UM, ax.get_ylim()[1] * 0.9, " 150 um vessel", fontsize=8, color="#555")
    ax.axvline(500, color="#27ae60", ls=":", lw=1)
    ax.text(500, ax.get_ylim()[1] * 0.7, " 0.5 mm voxel", fontsize=8, color="#27ae60")
    ax.set_xlabel("reconstruction voxel size [um]")
    ax.set_ylabel("iron conc. relative to homogeneous mean")
    ax.set_title("Vascular uptake: partial-volume averaging")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.3)
    figpath = str(conrad_backend.REPO_ROOT / "paper" / "figures" / "fig_study_b.pdf")
    fig.savefig(figpath, bbox_inches="tight")
    dashpng = str(conrad_backend.REPO_ROOT / "docs" / "assets" / "fig_study_b.png")
    fig.savefig(dashpng, dpi=140, bbox_inches="tight"); plt.close(fig)

    print("[Study B] mass conservation: ROI-mean(rel) =",
          {int(v*1000): round(pv[v]["roi_mean"], 3) for v in voxels})
    print("[Study B] per-voxel peak(rel):",
          {int(v*1000): round(pv[v]["voxel_peak"], 2) for v in voxels})
    print(f"[Study B] BH nonlinearity gap: realistic={bh['realistic_6mg']['hu_gap']:.3f} HU, "
          f"5x={bh['5x_dose']['hu_gap']:.3f} HU")
    print(f"[Study B] structural noise @0.5mm: {sn['n_vessels_per_voxel']:.1f} vessels/voxel, "
          f"per-voxel rel-std={sn['per_voxel_rel_std']:.2f}")
    print("[Study B] wrote", outdir, "and paper/figures/fig_study_b.pdf")
    return result


if __name__ == "__main__":
    main()
