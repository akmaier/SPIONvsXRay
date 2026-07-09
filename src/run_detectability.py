"""Detectability sweep: EID vs multi-bin PCD across iron concentration.

Base-material sinograms + the polychromatic accumulators depend only on
phantom/geometry/spectrum, so compute them ONCE; each cell is a cheap noise draw
+ CONRAD fan FBP + per-insert measurement. Reports per-insert iron ΔHU + CNR
(over EVAL.noise_realizations noise draws) and Rose detection thresholds per
detector. Water beam-hardening precorrection and the minimal short scan are ALWAYS
on (not factors); with_bone stays a factor -- a hard BH source, present in the
rabbit study (see run() docstring / DEVLOG 2026-07-09).
"""
from __future__ import annotations
import numpy as np

import conrad_backend
import conrad_ct
import conrad_phantom
import conrad_project as cp
import materials
import spectrum as spec
from config import DETECTORS, SPECTRUM, EVAL, PHANTOM

N0 = SPECTRUM.photons_per_pixel


def polychromatic_accumulators(base, kvp=None, filters=()):
    """Noise-free detector accumulators (energy sums), computed once."""
    if kvp is None:
        E, flux, _ = spec.standard_spectrum()
    else:
        E, flux = spec.conrad_spectrum(kvp)
    if filters:
        flux = spec.apply_filters(E, flux, list(filters))
    s = flux / flux.sum()
    names = list(base.keys())
    mu = cp._material_mu(names, E)
    L = {n: base[n] * 0.1 for n in names}          # mm -> cm
    shape = next(iter(base.values())).shape
    edges = np.array(DETECTORS.pcd_bin_edges_kev)
    nb = len(edges) - 1

    S_det = np.zeros(shape); S_det_E2 = np.zeros(shape)
    S_air = float(np.sum(N0 * s * E))
    C_det = [np.zeros(shape) for _ in range(nb)]
    C_air = [0.0] * nb
    for j, Ej in enumerate(E):
        tau = np.zeros(shape)
        for n in names:
            tau += L[n] * mu[n][j]
        nph = N0 * s[j] * np.exp(-tau)
        S_det += nph * Ej
        S_det_E2 += nph * Ej * Ej
        b = int(np.searchsorted(edges, Ej, side="right") - 1)
        if 0 <= b < nb:
            C_det[b] += nph
            C_air[b] += N0 * s[j]
    return dict(S_det=S_det, S_det_E2=S_det_E2, S_air=S_air,
                C_det=C_det, C_air=C_air, E=E, s=s, edges=edges)


def _pcd_weights(acc):
    """CNR-optimal per-bin weights for COUNT-domain combination (matched filter).

    The bins are combined as a weighted sum of photon COUNTS (like an EID, but
    with optimal weights instead of energy E), then a *single* log is taken —
    see `line_integral`. The CNR-optimal weight for that scheme is
        w_b ∝ S_b / V_b   with  S_b = Σ_{E∈bin} N_t(E)·c(E),  V_b = Σ_{E∈bin} N_t(E)
    i.e. the detected-count-weighted mean iron contrast in the bin. This attains
    the binned detectability ceiling ΣS_b²/V_b derived in src/spectral.py.

    N_t = N0·s·exp(−μ_tissue·L_body) is the background photon count reaching the
    detector (spectrum shape × body attenuation); c(E) is the per-energy iron-oxide
    contrast. The low bin (10–37.5 keV) has the highest contrast-per-count so it
    gets the largest weight, but because we combine COUNTS (not per-bin logs) its
    few surviving photons keep its contribution small and, crucially, finite.
    """
    E, edges, s = acc["E"], acc["edges"], acc["s"]
    L_body_cm = PHANTOM.body_diameter_mm / 10.0
    mu_tissue = materials.linear_attenuation("water", E)          # 1/cm
    ox = materials.oxide_contrast_massatten(E)                    # iron contrast shape
    Nt = N0 * s * np.exp(-mu_tissue * L_body_cm)                  # detected bkg photons/energy
    w = []
    for b in range(len(edges) - 1):
        m = (E >= edges[b]) & (E < edges[b + 1])
        S_b = float((Nt[m] * ox[m]).sum())
        V_b = float(Nt[m].sum())
        w.append(S_b / (V_b + 1e-12))                            # count-weighted mean contrast
    w = np.clip(np.array(w), 0.0, None)
    return w / (w.sum() + 1e-12)


def line_integral(acc, detector, seed, bh_polys=None):
    """One noisy line-integral sinogram for a detector model.

    Spectral processing (see README "Spectral processing"):
      EID: energy-integrated signal S = sum_E N(E)*E, Gaussian quantum noise with
           variance sum_E N(E)*E^2, single log; optional water precorrection via a
           single polynomial `bh_polys` calibrated for the EID spectrum.
      PCD: per energy bin b, Poisson counts C_b; each bin is beam-hardening-corrected
           on ITS OWN spectrum (per-bin polynomial `bh_polys[b]`, applied as a
           corrected count), then combined in the count domain M = sum_b w_b*C_b_corr
           (matched-filter weights) and a single log -> one FBP.
    """
    rng = np.random.default_rng(seed)
    eps = 1e-6
    if detector == "EID":
        S = acc["S_det"] + rng.normal(0.0, np.sqrt(np.maximum(acc["S_det_E2"], 1e-30)))
        p = -np.log(np.clip(S, eps, None) / acc["S_air"])
        return np.polyval(bh_polys, p) if bh_polys is not None else p
    # PCD: matched-filter COUNT-domain combination (robust to the photon-starved low
    # bin), with per-bin water precorrection applied as CORRECTED COUNTS
    # C_b_corr = C_air_b*exp(-poly_b(p_b)). Correcting each bin on its own spectrum is
    # the energy-dependent (per-bin) BH approach; recombining in the count domain
    # keeps the noise robustness. With bh_polys=None this is exactly the count-domain
    # combination M = sum_b w_b*C_b (poly_b = identity for the uncorrected case).
    w = _pcd_weights(acc)
    M = np.zeros_like(acc["S_det"]); M_air = 0.0
    for b, (cd, ca) in enumerate(zip(acc["C_det"], acc["C_air"])):
        cm = rng.poisson(np.maximum(cd, 0.0))
        if bh_polys is not None:
            p_b = -np.log(np.maximum(cm, 1.0) / max(ca, eps))
            cm = ca * np.exp(-np.polyval(bh_polys[b], p_b))     # per-bin corrected count
        M += w[b] * cm
        M_air += w[b] * ca
    return -np.log(np.clip(M, eps, None) / max(M_air, eps))


def bh_poly_for(acc, detector):
    """Calibrated water precorrection for a detector, applied in `line_integral`.

    EID: a single polynomial for the energy-weighted (s*E) spectrum.
    PCD: a LIST of per-bin polynomials, each calibrated for that bin's own spectrum.
    """
    E, s, edges = acc["E"], acc["s"], acc["edges"]
    mu_w = materials.linear_attenuation("water", E)     # 1/cm
    if detector == "EID":
        return conrad_ct.water_precorrection_poly(E, s * E, mu_w)
    return [conrad_ct.water_precorrection_poly(E, np.where((E >= edges[b]) & (E < edges[b + 1]), s, 0.0), mu_w)
            for b in range(len(edges) - 1)]


def run(bones=(False, True)):
    """Detectability study across the phantom factors. Each cell = mean iron signal
    / std of the local background over EVAL.noise_realizations noise draws.

    Water beam-hardening precorrection is ALWAYS applied (calibrated per detector /
    per PCD bin), and the acquisition is ALWAYS the minimal short scan (180 deg +
    fan) with Parker redundancy weighting. Neither is a factor anymore.

    Factors (see DEVLOG 2026-07-09):
      - detector: EID vs multi-bin PCD.
      - with_bone: cortical-bone rod absent / present. Bone is a HARD beam-hardening
        source, present in the rabbit study, so it is an explicit FACTOR here, not
        the baseline. Its streak fans across the ring inserts -- a phantom effect,
        for the discussion.
      - c_Fe concentration; EVAL.noise_realizations noise draws per cell.
    """
    rows = []
    inserts_ref = None
    for with_bone in bones:
        scene, inserts = conrad_phantom.build_phantom(with_bone=with_bone)
        inserts_ref = inserts
        geo = conrad_ct.fan_geometry(n_pix=512, short_scan=True)
        base, geo = cp.project_base_materials(inserts, geo)
        acc = polychromatic_accumulators(base)          # once per phantom
        spions = [i for i in inserts if i["c_form"] is not None and i["name"] != "SPION_c0"]
        for detector in DETECTORS.types:                # EID, PCD
            bh_polys = bh_poly_for(acc, detector)       # water precorrection always on
            signal = {i["name"]: [] for i in inserts}
            for seed in range(EVAL.noise_realizations):
                sino = line_integral(acc, detector, seed, bh_polys)
                recon = conrad_ct.fbp(sino, geo)
                meas = cp.measure_inserts(recon, geo, inserts)
                for m in meas:
                    signal[m["name"]].append((m["iron_delta_hu"], m["delta_hu"]))
            for ins in spions:
                arr = np.array(signal[ins["name"]])     # (reps, 2)
                d_hu = float(arr[:, 0].mean())          # c0-corrected iron signal
                noise = float(arr[:, 1].std())          # quantum noise on the insert
                cnr = d_hu / (noise + 1e-9)
                rows.append(dict(detector=detector, with_bone=with_bone,
                                 name=ins["name"], c_fe=ins["c_fe"],
                                 delta_hu=d_hu, noise=noise, cnr=cnr))
            print(f"[{detector} bone={int(with_bone)}] "
                  + "  ".join(f"{r['c_fe']:.2f}:CNR{r['cnr']:.1f}" for r in rows[-len(spions):]))
    return rows, inserts_ref


def thresholds(rows, rose=EVAL.rose_cnr_threshold):
    """Lowest detectable c_Fe (CNR >= Rose) per (detector, with_bone) cell."""
    out = {}
    for det, wb in sorted({(r["detector"], r["with_bone"]) for r in rows}):
        cells = sorted([r for r in rows if r["detector"] == det and r["with_bone"] == wb],
                       key=lambda r: r["c_fe"])
        out[f"{det}_bone{int(wb)}"] = next((r["c_fe"] for r in cells if r["cnr"] >= rose), None)
    return out


def save_results(rows, outdir=None):
    """Persist the detectability sweep to results/detectability/ (CSV + JSON)."""
    import csv, json, os
    if outdir is None:
        outdir = str(conrad_backend.REPO_ROOT / "results" / "detectability")
    os.makedirs(outdir, exist_ok=True)
    cols = ["detector", "with_bone", "name", "c_fe", "delta_hu", "noise", "cnr"]
    with open(os.path.join(outdir, "detectability.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in cols})
    meta = dict(n_realizations=EVAL.noise_realizations,
                photons_per_pixel=SPECTRUM.photons_per_pixel,
                pcd_bin_edges_kev=list(DETECTORS.pcd_bin_edges_kev),
                thresholds=thresholds(rows),
                thresholds_rose3=thresholds(rows, rose=3.0),
                rows=rows)
    with open(os.path.join(outdir, "detectability.json"), "w") as f:
        json.dump(meta, f, indent=2)
    return outdir


if __name__ == "__main__":
    import time
    t = time.time()
    rows, inserts = run()
    dt = time.time() - t
    print(f"\n[detectability] {len(rows)} cells, {EVAL.noise_realizations} reps each, in {dt:.0f}s")
    for rose in (3.0, 5.0):
        thr = thresholds(rows, rose=rose)
        print(f"Detection thresholds (lowest c_Fe with CNR>={rose:.0f}, mg Fe/ml):")
        for det, c in thr.items():
            print(f"  {det}: {c if c is not None else 'none (undetectable)'}")
    outdir = save_results(rows)
    print(f"[ok] wrote {outdir} -> detectability.csv, detectability.json")
