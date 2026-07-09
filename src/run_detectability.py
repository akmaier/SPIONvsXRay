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
from config import (DETECTORS, SPECTRUM, EVAL, PHANTOM, DOSE_LEVELS,
                    CELL_DENSITY_LEVELS)

N0 = SPECTRUM.photons_per_pixel


def load_spectrum(spectrum_name):
    """Resolve a named study spectrum -> (kvp, filters, bin_edges, label).

    Reads the E2 optimum (results/spectral/optimum.json, written by
    run_spectral_sweep.py). spectrum_name is "optimum" (the CNR-optimal kVp/Al +
    re-optimized PCD bins) or "baseline" (90 kVp Al2.5, default PCD bins). Both the
    studies run at these two spectra so the optimum can be compared to the C-arm
    baseline. Falls back to the config standard spectrum if optimum.json is absent.
    """
    import json
    path = conrad_backend.REPO_ROOT / "results" / "spectral" / "optimum.json"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found -- run src/run_spectral_sweep.py (E2) first")
    with open(path) as f:
        opt = json.load(f)
    if spectrum_name == "optimum":
        kvp = float(opt["kvp"])
        filters = tuple(tuple(x) for x in opt["filters"])
        edges = tuple(opt["pcd_bin_edges_kev"])
        label = f"{int(kvp)}kVp_{opt['filter']}"
    elif spectrum_name == "baseline":
        b = opt["baseline"]
        kvp = float(b["kvp"])
        filters = tuple(tuple(x) for x in b["filters"])
        edges = tuple(DETECTORS.pcd_bin_edges_kev)   # default PCD bins for the baseline
        label = f"{int(kvp)}kVp_{b['filter']}"
    else:
        raise ValueError(f"unknown spectrum {spectrum_name!r} (use 'optimum'/'baseline')")
    return dict(name=spectrum_name, kvp=kvp, filters=filters,
                bin_edges=edges, label=label)


def polychromatic_accumulators(base, kvp=None, filters=(), n0=None, bin_edges=None):
    """Noise-free detector accumulators (energy sums), computed once per dose.

    `n0` is the unattenuated I0 [photons/pixel] driving the signal level (and thus
    the Poisson/Gaussian noise floor). Defaults to the documented reference dose
    (SPECTRUM.photons_per_pixel); the dose factor passes DOSE_LEVELS values.

    `bin_edges` overrides the PCD energy-bin edges [keV]; defaults to
    DETECTORS.pcd_bin_edges_kev. The E2 optimum re-optimizes these per spectrum
    (results/spectral/optimum.json) and the studies pass them through here so the
    matched-filter weights / BH precorrection use the winning spectrum's bins.
    """
    if n0 is None:
        n0 = N0
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
    if bin_edges is None:
        bin_edges = DETECTORS.pcd_bin_edges_kev
    edges = np.array(bin_edges, dtype=float)
    nb = len(edges) - 1

    S_det = np.zeros(shape); S_det_E2 = np.zeros(shape)
    S_air = float(np.sum(n0 * s * E))
    C_det = [np.zeros(shape) for _ in range(nb)]
    C_air = [0.0] * nb
    for j, Ej in enumerate(E):
        tau = np.zeros(shape)
        for n in names:
            tau += L[n] * mu[n][j]
        nph = n0 * s[j] * np.exp(-tau)
        S_det += nph * Ej
        S_det_E2 += nph * Ej * Ej
        b = int(np.searchsorted(edges, Ej, side="right") - 1)
        if 0 <= b < nb:
            C_det[b] += nph
            C_air[b] += n0 * s[j]
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
    mu_tissue = materials.linear_attenuation_soft(E)             # ICRU soft tissue, 1/cm
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


def line_integral(acc, detector, seed, bh_polys=None, noiseless=False):
    """One line-integral sinogram for a detector model.

    noiseless=True skips the Poisson/Gaussian draw entirely (mean counts), for a
    TRUE noise-free reconstruction -- so EID and PCD are compared on identical, exactly
    noiseless signals (the PCD path otherwise always draws Poisson counts).

    Spectral processing (see README "Spectral processing"):
      EID: energy-integrated signal S = sum_E N(E)*E, Gaussian quantum noise with
           variance sum_E N(E)*E^2, single log; optional water precorrection via a
           single polynomial `bh_polys` calibrated for the EID spectrum.
      PCD: per energy bin b, Poisson counts C_b are summed in the intensity (I0) domain
           with matched-filter weights (M = sum_b w_b*C_b), a single log gives the combined
           line integral, then ONE water precorrection for the combined effective spectrum
           (bh_poly_for) linearizes it. Correcting AFTER the count sum keeps water flat AND
           weights the photon-starved low bin by its actual count -> CNR ~1.2x EID (a per-bin
           pre-sum correction instead inflates the low bin's noise, giving only CNR ~1.0x).
    """
    rng = np.random.default_rng(seed)
    eps = 1e-6
    if detector == "EID":
        S = acc["S_det"]
        if not noiseless:
            S = S + rng.normal(0.0, np.sqrt(np.maximum(acc["S_det_E2"], 1e-30)))
        p = -np.log(np.clip(S, eps, None) / acc["S_air"])
        return np.polyval(bh_polys, p) if bh_polys is not None else p
    # PCD: sum the RAW bin counts in the INTENSITY (I0) domain with matched-filter weights
    # (M = sum_b w_b*C_b), single log, then ONE water precorrection on the combined line
    # integral (bh_polys = the single combined-spectrum poly). Correcting after the sum
    # keeps water flat while the count-domain sum keeps the starved low bin weighted by its
    # actual count -> PCD CNR ~1.2x EID (per-bin pre-sum correction gives only ~1.0x).
    w = _pcd_weights(acc)
    M = np.zeros_like(acc["S_det"]); M_air = 0.0
    for b, (cd, ca) in enumerate(zip(acc["C_det"], acc["C_air"])):
        cm = cd if noiseless else rng.poisson(np.maximum(cd, 0.0))
        M += w[b] * cm
        M_air += w[b] * ca
    p = -np.log(np.clip(M, eps, None) / max(M_air, eps))
    return np.polyval(bh_polys, p) if bh_polys is not None else p


def bh_poly_for(acc, detector):
    """Calibrated water precorrection for a detector, applied in `line_integral`.

    EID: a single polynomial for the energy-weighted (s*E) spectrum.
    PCD: a SINGLE polynomial for the count-combined effective spectrum
         w_eff(E) = w[bin(E)]*s(E) -- applied AFTER the I0-domain count sum (like the EID),
         so the combined line integral is water-flat. (Correcting per bin BEFORE the sum
         instead inflates the photon-starved low bin's noise; see line_integral / DEVLOG.)
    """
    E, s, edges = acc["E"], acc["s"], acc["edges"]
    mu_w = materials.linear_attenuation("water", E)     # 1/cm
    if detector == "EID":
        return conrad_ct.water_precorrection_poly(E, s * E, mu_w)
    w = _pcd_weights(acc)                                # matched-filter bin weights
    w_eff = np.zeros_like(s)
    for b in range(len(edges) - 1):
        m = (E >= edges[b]) & (E < edges[b + 1])
        w_eff[m] = w[b] * s[m]                           # bin weight x spectrum = combined weighting
    return conrad_ct.water_precorrection_poly(E, w_eff, mu_w)


def _vessel_effects(c_fe, d_hu, quantum_noise, roi_mm=8.0, kvp=None, filters=()):
    """Second-order vessel (Study B) corrections applied to a homogenized insert.

    SPEC §5.9: iron in a fresh injection is confined to 150 um vessels at 10% volume
    fraction (=> 10x LOCAL concentration). The 0.5 mm CT voxel cannot resolve the
    vessels, so the ROI-mean iron (and thus first-order ΔHU) is UNCHANGED (mass
    conserved). This adds the two second-order effects the homogeneous model misses,
    reusing the Study B primitives (src/run_study_b.py):

      (1) Beam-hardening nonlinearity (Jensen gap): rays through the local-conc vessel
          segments harden more than proportionally, so the ROI ΔHU shifts by the
          homogeneous-vs-vessel line-integral gap -- computed at THIS insert's OWN mean
          iron (c_fe) with local = c_fe / VESSEL_VOLUME_FRACTION. No hard-coded
          delivered-mass reference (the old C_FE_MEAN scaling is gone).
      (2) Structural noise: the per-voxel vessel fraction is binomial, adding a
          texture std that survives ROI averaging (reduced by sqrt(n_voxels_ROI)).
          Added in quadrature to the quantum noise on the insert.

    The gap is evaluated on the study's ACTUAL spectrum (kvp/filters), matching the
    reconstruction. Returns (d_hu_vessel, noise_vessel).
    """
    import run_study_b as sb
    from config import VESSEL_VOLUME_FRACTION
    local_mult = 1.0 / VESSEL_VOLUME_FRACTION                    # 10x local conc
    # (1) BH nonlinearity: polychromatic homogeneous-vs-vessel line-integral gap,
    # evaluated at THIS insert's OWN mean iron (no delivered-mass reference).
    c_mean = 1e-3 * c_fe                                          # g/cm^3 (tumor mean)
    if kvp is None:
        E, flux, _ = spec.standard_spectrum()
    else:
        E, flux = spec.conrad_spectrum(kvp)
    if filters:
        flux = spec.apply_filters(E, flux, list(filters))
    s = flux / flux.sum()
    mu_t = materials.linear_attenuation_soft(E)                  # soft-tissue matrix
    ox = materials.oxide_contrast_massatten(E)
    def _p(profile):
        tau = mu_t * sb.L_BODY_CM
        for length_cm, c in profile:
            tau = tau + (mu_t * length_cm + c * ox * length_cm)
        return -np.log(np.sum(N0 * s * np.exp(-tau) * E) / np.sum(N0 * s * E))
    p_hom = _p([(sb.L_TUMOR_CM, c_mean)])
    p_ves = _p([(VESSEL_VOLUME_FRACTION * sb.L_TUMOR_CM, c_mean * local_mult),
                ((1 - VESSEL_VOLUME_FRACTION) * sb.L_TUMOR_CM, 0.0)])
    mu_w_eff = float(np.sum(s * mu_t))
    hu_gap = 1000.0 * (p_ves - p_hom) / (mu_w_eff * sb.L_TUMOR_CM)
    d_hu_vessel = d_hu + hu_gap

    # (2) structural noise: per-voxel binomial texture, reduced by ROI averaging.
    sn = sb.structural_noise(voxel_mm=cp_voxel_mm())
    n_vox_roi = max(1.0, (roi_mm / cp_voxel_mm()) ** 2 * np.pi)   # ~circular ROI voxel count
    struct_hu = abs(d_hu) * sn["per_voxel_rel_std"] / np.sqrt(n_vox_roi)
    noise_vessel = float(np.hypot(quantum_noise, struct_hu))
    return float(d_hu_vessel), noise_vessel


def cp_voxel_mm():
    """Recon voxel size [mm] used by the fan FBP grid (config RECON_VOXEL_MM)."""
    from config import RECON_VOXEL_MM
    return RECON_VOXEL_MM


def run(bones=(False, True), doses=None, tumor_models=("homogeneous", "vessel")):
    """Detectability study across the phantom factors. Each cell = mean iron signal
    / std of the local background over EVAL.noise_realizations noise draws.

    Beam-hardening precorrection is ALWAYS applied (calibrated per detector / per PCD
    bin), and the acquisition is ALWAYS the minimal short scan (180 deg + fan) with
    Parker redundancy weighting. The spectrum is FIXED at the 90 kVp standard. None
    of these is a factor.

    Factors (SPEC §5.6/§5.9):
      - detector: EID vs multi-bin PCD.
      - with_bone: cortical-bone rod absent / present (a hard beam-hardening source).
      - dose: DOSE_LEVELS (low/high) photons/pixel -- parameterizes the Poisson/
        Gaussian noise floor (low dose => higher noise => lower CNR).
      - tumor_model: two experiments -- homogeneous (Study A, cellular uptake:
        SPIONs internalised -> ~uniform tumor iron) vs vessel (Study B, vascular/
        fresh delivery: SPIONs still in 150 um vessels at 10% volume fraction, not
        yet taken up => 10x local conc; mass-conserved mean, plus BH nonlinearity +
        structural noise, reusing src/run_study_b.py).
      - c_Fe concentration; EVAL.noise_realizations noise draws per cell.
    """
    if doses is None:
        doses = DOSE_LEVELS
    rows = []
    inserts_ref = None
    for with_bone in bones:
        scene, inserts = conrad_phantom.build_phantom(with_bone=with_bone)
        inserts_ref = inserts
        geo = conrad_ct.fan_geometry(n_pix=512, short_scan=True)
        base, geo = cp.project_base_materials(inserts, geo)
        spions = [i for i in inserts if i["c_form"] is not None and i["name"] != "SPION_c0"]
        for dose_name, n0 in doses.items():
            acc = polychromatic_accumulators(base, n0=n0)   # once per (phantom, dose)
            for detector in DETECTORS.types:                # EID, PCD
                bh_polys = bh_poly_for(acc, detector)       # water precorrection always on
                signal = {i["name"]: [] for i in inserts}
                for seed in range(EVAL.noise_realizations):
                    sino = line_integral(acc, detector, seed, bh_polys)
                    recon = conrad_ct.fbp(sino, geo)
                    meas = cp.measure_inserts(recon, geo, inserts)
                    for m in meas:
                        signal[m["name"]].append((m["iron_delta_hu"], m["delta_hu"]))
                for tumor_model in tumor_models:
                    for ins in spions:
                        arr = np.array(signal[ins["name"]])     # (reps, 2)
                        d_hu = float(arr[:, 0].mean())          # c0-corrected iron signal
                        noise = float(arr[:, 1].std())          # quantum noise on the insert
                        if tumor_model == "vessel":
                            d_hu, noise = _vessel_effects(ins["c_fe"], d_hu, noise)
                        cnr = d_hu / (noise + 1e-9)
                        rows.append(dict(detector=detector, with_bone=with_bone,
                                         dose=dose_name, tumor_model=tumor_model,
                                         name=ins["name"], c_fe=ins["c_fe"],
                                         delta_hu=d_hu, noise=noise, cnr=cnr))
                    print(f"[{detector} bone={int(with_bone)} dose={dose_name} {tumor_model}] "
                          + "  ".join(f"{r['c_fe']:.2f}:CNR{r['cnr']:.1f}" for r in rows[-len(spions):]))
    return rows, inserts_ref


def run_study_a(spectrum_name="optimum", n_reps=None, densities=None,
                doses=None, bones=(False, True)):
    """STUDY A -- homogeneous cellular uptake (SPIONs internalised -> ~uniform iron).

    Factors: detector {EID,PCD} x cell_density (config.CELL_DENSITY_LEVELS) x dose
    (config.DOSE_LEVELS) x with_bone {False,True}. The 8 loading configs (4 SPION
    formulations x {0h,24h}) are the per-cell inserts, over
    build_phantom(loading=True, density=d, with_bone=b).

    Run at a named spectrum (load_spectrum): "optimum" (E2 winner + re-optimized PCD
    bins) or "baseline" (90 kVp Al2.5). BH water precorrection + minimal short scan
    are ALWAYS on. n_reps overrides EVAL.noise_realizations (smoke test only).

    Returns (rows, spec_info). Per config: iron ΔHU (c0-corrected) + CNR over the
    noise draws; rows carry detector, cell_density, dose, with_bone, config, c_fe,
    delta_hu, noise, cnr.
    """
    if n_reps is None:
        n_reps = EVAL.noise_realizations
    if densities is None:
        densities = CELL_DENSITY_LEVELS
    if doses is None:
        doses = DOSE_LEVELS
    sp = load_spectrum(spectrum_name)
    rows = []
    for dens_name, dens in densities.items():
        for with_bone in bones:
            scene, inserts = conrad_phantom.build_phantom(loading=True, density=dens,
                                                          with_bone=with_bone)
            geo = conrad_ct.fan_geometry(n_pix=512, short_scan=True)
            base, geo = cp.project_base_materials(inserts, geo)
            spions = [i for i in inserts if i["c_form"] is None and i["name"] != "SPION_c0"]
            for dose_name, n0 in doses.items():
                acc = polychromatic_accumulators(base, kvp=sp["kvp"], filters=sp["filters"],
                                                 n0=n0, bin_edges=sp["bin_edges"])
                for detector in DETECTORS.types:
                    bh_polys = bh_poly_for(acc, detector)     # water precorrection always on
                    signal = {i["name"]: [] for i in inserts}
                    for seed in range(n_reps):
                        sino = line_integral(acc, detector, seed, bh_polys)
                        recon = conrad_ct.fbp(sino, geo)
                        for m in cp.measure_inserts(recon, geo, inserts):
                            signal[m["name"]].append((m["iron_delta_hu"], m["delta_hu"]))
                    for ins in spions:
                        arr = np.array(signal[ins["name"]])
                        d_hu = float(arr[:, 0].mean())
                        noise = float(arr[:, 1].std())
                        rows.append(dict(detector=detector, cell_density=dens_name,
                                         dose=dose_name, with_bone=with_bone,
                                         config=ins["name"], c_fe=ins["c_fe"],
                                         delta_hu=d_hu, noise=noise,
                                         cnr=d_hu / (noise + 1e-9)))
                    print(f"[A {sp['label']} {detector} dens={dens_name} "
                          f"dose={dose_name} bone={int(with_bone)}] "
                          + "  ".join(f"{r['config']}:c{r['c_fe']:.2f}:CNR{r['cnr']:.1f}"
                                      for r in rows[-len(spions):]))
    return rows, sp


def run_study_b(spectrum_name="optimum", n_reps=None, doses=None, bones=(False, True)):
    """STUDY B -- fresh vascular injection (SPIONs still in 150 um vessels @10%).

    Factors: detector {EID,PCD} x dose (config.DOSE_LEVELS) x with_bone {False,True},
    over build_phantom(study_b=True, with_bone=b). Each insert is HOMOGENIZED to its
    tumor-mean iron (c_fe = vessel_level * VESSEL_VOLUME_FRACTION); the vessel
    second-order model (BH nonlinearity + structural noise, _vessel_effects) then uses
    each insert's OWN mean c_fe with local = mean / VESSEL_VOLUME_FRACTION (10x). No
    hard-coded delivered-mass reference.

    Same spectrum handling / always-on BH precorrection + short scan as Study A.
    Returns (rows, spec_info); rows carry detector, dose, with_bone, config,
    vessel_level, c_fe, delta_hu, noise, cnr.
    """
    if n_reps is None:
        n_reps = EVAL.noise_realizations
    if doses is None:
        doses = DOSE_LEVELS
    sp = load_spectrum(spectrum_name)
    rows = []
    for with_bone in bones:
        scene, inserts = conrad_phantom.build_phantom(study_b=True, with_bone=with_bone)
        geo = conrad_ct.fan_geometry(n_pix=512, short_scan=True)
        base, geo = cp.project_base_materials(inserts, geo)
        spions = [i for i in inserts if i.get("vessel_level", 0.0) > 0.0]
        for dose_name, n0 in doses.items():
            acc = polychromatic_accumulators(base, kvp=sp["kvp"], filters=sp["filters"],
                                             n0=n0, bin_edges=sp["bin_edges"])
            for detector in DETECTORS.types:
                bh_polys = bh_poly_for(acc, detector)
                signal = {i["name"]: [] for i in inserts}
                for seed in range(n_reps):
                    sino = line_integral(acc, detector, seed, bh_polys)
                    recon = conrad_ct.fbp(sino, geo)
                    for m in cp.measure_inserts(recon, geo, inserts):
                        signal[m["name"]].append((m["iron_delta_hu"], m["delta_hu"]))
                for ins in spions:
                    arr = np.array(signal[ins["name"]])
                    d_hu = float(arr[:, 0].mean())
                    noise = float(arr[:, 1].std())
                    # vessel second-order model on THIS insert's own mean iron
                    d_hu, noise = _vessel_effects(ins["c_fe"], d_hu, noise,
                                                  kvp=sp["kvp"], filters=sp["filters"])
                    rows.append(dict(detector=detector, dose=dose_name,
                                     with_bone=with_bone, config=ins["name"],
                                     vessel_level=ins["vessel_level"], c_fe=ins["c_fe"],
                                     delta_hu=d_hu, noise=noise,
                                     cnr=d_hu / (noise + 1e-9)))
                print(f"[B {sp['label']} {detector} dose={dose_name} bone={int(with_bone)}] "
                      + "  ".join(f"v{r['vessel_level']:.0f}:c{r['c_fe']:.2f}:CNR{r['cnr']:.1f}"
                                  for r in rows[-len(spions):]))
    return rows, sp


def _study_thresholds(rows, cellkeys, rose):
    """Lowest detectable c_Fe (CNR>=rose) per cell (grouped by cellkeys)."""
    out = {}
    keyset = sorted({tuple(r[k] for k in cellkeys) for r in rows})
    for kv in keyset:
        cells = sorted([r for r in rows if all(r[k] == v for k, v in zip(cellkeys, kv))],
                       key=lambda r: r["c_fe"])
        label = "_".join(f"{v}" for v in kv)
        out[label] = next((r["c_fe"] for r in cells if r["cnr"] >= rose), None)
    return out


def save_study(rows, sp, study, cellkeys, outdir=None):
    """Persist a study sweep to results/detectability/study_<x>_<spectrum>.json (+csv)."""
    import csv, json, os
    if outdir is None:
        outdir = str(conrad_backend.REPO_ROOT / "results" / "detectability")
    os.makedirs(outdir, exist_ok=True)
    stem = f"study_{study}_{sp['name']}"
    cols = list(rows[0].keys()) if rows else []
    with open(os.path.join(outdir, stem + ".csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for r in rows:
            w.writerow(r)
    meta = dict(study=study, spectrum=sp, n_realizations=EVAL.noise_realizations,
                dose_levels=dict(DOSE_LEVELS),
                thresholds_rose5=_study_thresholds(rows, cellkeys, 5.0),
                thresholds_rose3=_study_thresholds(rows, cellkeys, 3.0),
                rows=rows)
    with open(os.path.join(outdir, stem + ".json"), "w") as f:
        json.dump(meta, f, indent=2)
    return os.path.join(outdir, stem + ".json")


def thresholds(rows, rose=EVAL.rose_cnr_threshold):
    """Lowest detectable c_Fe (CNR >= Rose) per (detector, bone, dose, tumor_model) cell."""
    out = {}
    keys = sorted({(r["detector"], r["with_bone"], r["dose"], r["tumor_model"]) for r in rows})
    for det, wb, dose, tm in keys:
        cells = sorted([r for r in rows if r["detector"] == det and r["with_bone"] == wb
                        and r["dose"] == dose and r["tumor_model"] == tm],
                       key=lambda r: r["c_fe"])
        out[f"{det}_bone{int(wb)}_{dose}_{tm}"] = next(
            (r["c_fe"] for r in cells if r["cnr"] >= rose), None)
    return out


def save_results(rows, outdir=None):
    """Persist the detectability sweep to results/detectability/ (CSV + JSON)."""
    import csv, json, os
    if outdir is None:
        outdir = str(conrad_backend.REPO_ROOT / "results" / "detectability")
    os.makedirs(outdir, exist_ok=True)
    cols = ["detector", "with_bone", "dose", "tumor_model", "name", "c_fe",
            "delta_hu", "noise", "cnr"]
    with open(os.path.join(outdir, "detectability.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in cols})
    meta = dict(n_realizations=EVAL.noise_realizations,
                photons_per_pixel=SPECTRUM.photons_per_pixel,
                dose_levels=dict(DOSE_LEVELS),
                pcd_bin_edges_kev=list(DETECTORS.pcd_bin_edges_kev),
                thresholds=thresholds(rows),
                thresholds_rose3=thresholds(rows, rose=3.0),
                rows=rows)
    with open(os.path.join(outdir, "detectability.json"), "w") as f:
        json.dump(meta, f, indent=2)
    return outdir


def main(spectra=("optimum", "baseline")):
    """Run Study A and Study B at each named spectrum and persist per-study results.

    This is the FULL detectability run (EVAL.noise_realizations reps). For a quick
    end-to-end check use run_study_a / run_study_b with restricted factors + n_reps.
    """
    import time
    t = time.time()
    for spectrum in spectra:
        rowsA, spA = run_study_a(spectrum)
        pA = save_study(rowsA, spA, "a",
                        ["detector", "cell_density", "dose", "with_bone"])
        rowsB, spB = run_study_b(spectrum)
        pB = save_study(rowsB, spB, "b", ["detector", "dose", "with_bone"])
        print(f"[ok] {spectrum}: wrote {pA} ({len(rowsA)} rows), "
              f"{pB} ({len(rowsB)} rows)")
    print(f"[detectability] done in {time.time() - t:.0f}s")


if __name__ == "__main__":
    main()
