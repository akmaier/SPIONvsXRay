"""Spectral optimization for iron (SPION) CT contrast.

Iron has no usable K-edge (7.1 keV), so contrast is pure photoelectric and lives
at low energy. This module quantifies, for the realistic 6 mg dose through the
rabbit, (a) the CNR-optimal monochromatic energy, (b) how added filters / kVp
reshape achievable CNR, and (c) the optimal photon-counting energy thresholds.

Framework: for an energy-resolving detector the maximum detectability (matched
filter / optimal energy weighting, Tapiovaara-Wagner / Cahn) is

    SNR^2_ideal = sum_E N_t(E) * c(E)^2

with background transmitted photons N_t(E) and per-energy line-integral contrast
c(E) = dmu(E) * L_tumor. For a weighting w(E):

    SNR^2(w) = (sum N_t w c)^2 / (sum N_t w^2)

so energy-integrating uses w=E, photon-counting-unweighted uses w=1, and an
M-bin detector with optimal per-bin weights achieves
    SNR^2_bins = sum_bins (sum_bin N_t c)^2 / (sum_bin N_t).
Thresholds are optimized by maximizing that sum.
"""
from __future__ import annotations
import numpy as np

import materials
from config import tumor_iron_conc

# Geometry for the detection task (SPEC §5.1): background ray ~10 cm soft tissue,
# tumor chord = tumor diameter.
L_BODY_CM = 10.0
L_TUMOR_CM = 2.48
C_FE_REALISTIC = 1e-3 * tumor_iron_conc(10.0)   # g/cm^3 at the 6 mg dose
N0_AIR = 70000.0                                # photons/pixel in air (post-filter)


def bremsstrahlung(E, kvp):
    """Kramers thin-target bremsstrahlung shape (unnormalized), 0 above kVp."""
    s = np.where(E < kvp, (kvp - E) / E, 0.0)
    return np.clip(s, 0.0, None)


def filter_transmission(E, filters):
    """Product of exp(-mu(E) t) for a list of (material_name, thickness_mm)."""
    T = np.ones_like(E)
    for name, t_mm in filters:
        mu = materials.linear_attenuation(name, E)   # 1/cm
        T *= np.exp(-mu * (t_mm / 10.0))             # mm -> cm
    return T


def spectrum(E, kvp=120.0, filters=(("aluminium", 2.5),)):
    """Normalized photon energy distribution s(E) (sums to 1)."""
    s = bremsstrahlung(E, kvp) * filter_transmission(E, list(filters))
    tot = np.trapezoid(s, E)
    return s / tot if tot > 0 else s


def _metrics(E, s):
    """Return CNR figures for a given normalized spectrum s(E)."""
    mu_tissue = materials.linear_attenuation("water", E)          # 1/cm
    murho_ox = materials.oxide_contrast_massatten(E)             # magnetite, per g Fe
    Nt = N0_AIR * s * np.exp(-mu_tissue * L_BODY_CM)              # transmitted bkg photons/energy
    c = C_FE_REALISTIC * murho_ox * L_TUMOR_CM                    # per-energy contrast (line integral)

    dE = np.gradient(E)
    Nt_d = Nt * dE                                               # photons per bin
    ideal = np.sum(Nt_d * c**2)
    eid = np.sum(Nt_d * E * c) ** 2 / np.sum(Nt_d * E**2)
    pcd1 = np.sum(Nt_d * c) ** 2 / np.sum(Nt_d)
    return dict(E=E, Nt=Nt_d, c=c, ideal=ideal, eid=eid, pcd1=pcd1)


def snr2_bins(E, Nt_d, c, thresholds):
    """Optimal-weight SNR^2 for bins defined by interior thresholds (keV)."""
    edges = [E[0] - 1] + list(thresholds) + [E[-1] + 1]
    tot = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (E > lo) & (E <= hi)
        V = np.sum(Nt_d[m])
        if V > 0:
            S = np.sum(Nt_d[m] * c[m])
            tot += S * S / V
    return tot


def optimize_thresholds(E, Nt_d, c, n_bins, grid=None):
    """Brute-force optimal interior thresholds maximizing binned SNR^2."""
    if grid is None:
        grid = np.arange(25, 115, 2.5)
    best, best_snr2 = None, -1.0
    if n_bins == 2:
        for t in grid:
            v = snr2_bins(E, Nt_d, c, [t])
            if v > best_snr2:
                best_snr2, best = v, [t]
    elif n_bins == 3:
        for i, t1 in enumerate(grid):
            for t2 in grid[i + 1:]:
                v = snr2_bins(E, Nt_d, c, [t1, t2])
                if v > best_snr2:
                    best_snr2, best = v, [t1, t2]
    return best, best_snr2


def report():
    materials._api()
    E = np.arange(15.0, 120.5, 0.5)

    print("== CNR-optimal monochromatic energy (per-photon detectability) ==")
    mu_tissue = materials.linear_attenuation("water", E)
    murho_ox = materials.oxide_contrast_massatten(E)
    # detectability density for a flat (fixed air-photon) source, mono:
    d = np.exp(-mu_tissue * L_BODY_CM) * (murho_ox ** 2)
    print(f"  argmax over 15-120 keV: E* = {E[np.argmax(d)]:.1f} keV")

    print("\n== Filter / kVp comparison (relative ideal CNR at fixed air flux) ==")
    configs = [
        ("120 kVp, 2.5 mm Al (baseline)", 120, (("aluminium", 2.5),)),
        ("120 kVp, 0.5 mm Al (soft)",     120, (("aluminium", 0.5),)),
        ("80 kVp, 2.5 mm Al",              80, (("aluminium", 2.5),)),
        ("60 kVp, 2.5 mm Al",              60, (("aluminium", 2.5),)),
        ("120 kVp, +0.3 mm Cu (hard)",    120, (("aluminium", 2.5), ("copper", 0.3))),
        ("120 kVp, +0.5 mm Sn (hard)",    120, (("aluminium", 2.5), ("tin", 0.5))),
        ("120 kVp, 0.1 mm Er (quasi-mono)",120,(("aluminium", 1.0), ("erbium", 0.1))),
    ]
    base = None
    for label, kvp, filt in configs:
        s = spectrum(E, kvp, filt)
        mm = _metrics(E, s)
        cnr = np.sqrt(mm["ideal"])
        if base is None:
            base = cnr
        print(f"  {label:34s}  ideal-CNR = {cnr/base:5.2f}x   "
              f"EID/ideal={np.sqrt(mm['eid']/mm['ideal']):.2f}  "
              f"PCD1/ideal={np.sqrt(mm['pcd1']/mm['ideal']):.2f}")

    print("\n== Optimal PCD thresholds (120 kVp, 2.5 mm Al, through rabbit) ==")
    s = spectrum(E, 120, (("aluminium", 2.5),))
    mm = _metrics(E, s)
    for nb in (2, 3):
        thr, snr2 = optimize_thresholds(mm["E"], mm["Nt"], mm["c"], nb)
        gain_vs_eid = np.sqrt(snr2 / mm["eid"])
        gain_vs_pcd1 = np.sqrt(snr2 / mm["pcd1"])
        frac_ideal = np.sqrt(snr2 / mm["ideal"])
        print(f"  {nb}-bin optimal thresholds = {[round(t,1) for t in thr]} keV | "
              f"CNR vs EID = {gain_vs_eid:.2f}x, vs PCD-unweighted = {gain_vs_pcd1:.2f}x, "
              f"= {frac_ideal*100:.0f}% of ideal")
    print(f"  [ref] EID CNR / ideal = {np.sqrt(mm['eid']/mm['ideal']):.2f}, "
          f"optimal-weight ceiling = 1.00")


def figures(outdir: str = None):
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import conrad_backend
    if outdir is None:
        outdir = str(conrad_backend.REPO_ROOT / "results" / "spectral")
    os.makedirs(outdir, exist_ok=True)
    E = np.arange(15.0, 120.5, 0.5)

    # (1) spectra under different filters/kVp
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for label, kvp, filt in [
        ("120 kVp, 2.5 mm Al", 120, (("aluminium", 2.5),)),
        ("60 kVp, 2.5 mm Al", 60, (("aluminium", 2.5),)),
        ("120 kVp +0.5 mm Sn", 120, (("aluminium", 2.5), ("tin", 0.5))),
    ]:
        ax.plot(E, spectrum(E, kvp, filt), lw=1.8, label=label)
    ax.set_xlabel("energy [keV]"); ax.set_ylabel("normalized photon fluence")
    ax.set_title("Beam spectra: filter / kVp shaping"); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(f"{outdir}/spectrum.png", dpi=130); plt.close(fig)

    # (2) per-energy detectability density with optimal PCD thresholds
    s = spectrum(E, 120, (("aluminium", 2.5),)); mm = _metrics(E, s)
    d = mm["Nt"] * mm["c"] ** 2
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.fill_between(E, d, color="#c0392b", alpha=0.5)
    ax.plot(E, d, color="#c0392b", lw=1.5, label=r"iron CNR$^2$ density")
    for t in (35.0, 47.5):
        ax.axvline(t, color="#2c3e50", ls="--", lw=1.2)
    ax.text(24, ax.get_ylim()[1]*0.8, "bin 1\n(high contrast)", fontsize=8, ha="center")
    ax.text(41, ax.get_ylim()[1]*0.8, "bin 2", fontsize=8, ha="center")
    ax.text(80, ax.get_ylim()[1]*0.8, "bin 3\n(low contrast)", fontsize=8, ha="center")
    ax.set_xlabel("energy [keV]"); ax.set_ylabel("detectability contribution")
    ax.set_title("Where iron contrast lives + optimal PCD thresholds (35 / 47.5 keV)")
    ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(f"{outdir}/pcd_bins.png", dpi=130); plt.close(fig)
    print("[ok] wrote", outdir, "-> spectrum.png, pcd_bins.png")


if __name__ == "__main__":
    report()
    figures()
