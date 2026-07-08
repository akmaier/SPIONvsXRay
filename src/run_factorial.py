"""M6b — full factorial: EID vs multi-bin PCD x beam-hardening off/on x 10 noise.

Base-material sinograms + the polychromatic accumulators depend only on
phantom/geometry/spectrum, so compute them ONCE; each of the 40 cells is a cheap
noise draw + CONRAD fan FBP + per-insert measurement. Reports per-insert iron
ΔHU + CNR (from the 10 noise realizations) and Rose detection thresholds.
"""
from __future__ import annotations
import numpy as np

import conrad_backend
import conrad_ct
import conrad_phantom
import conrad_project as cp
import materials
import spectrum as spec
from config import DETECTORS, SPECTRUM, EVAL

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
    """Optimal-ish per-bin weights ~ mean iron (oxide) contrast in each bin."""
    ox = materials.oxide_contrast_massatten(acc["E"])
    edges = acc["edges"]; s = acc["s"]
    w = []
    for b in range(len(edges) - 1):
        m = (acc["E"] >= edges[b]) & (acc["E"] < edges[b + 1])
        w.append(float((ox[m] * s[m]).sum() / (s[m].sum() + 1e-12)))
    w = np.array(w)
    return w / w.sum()


def line_integral(acc, detector, seed):
    """One noisy line-integral sinogram for a detector model."""
    rng = np.random.default_rng(seed)
    eps = 1e-6
    if detector == "EID":
        S = acc["S_det"] + rng.normal(0.0, np.sqrt(np.maximum(acc["S_det_E2"], 1e-30)))
        return -np.log(np.clip(S, eps, None) / acc["S_air"])
    # PCD: per-bin Poisson counts -> per-bin line integrals -> optimal weighting
    w = _pcd_weights(acc)
    p = np.zeros_like(acc["S_det"])
    for b, (cd, ca) in enumerate(zip(acc["C_det"], acc["C_air"])):
        cm = rng.poisson(np.maximum(cd, 0.0))
        p += w[b] * (-np.log(np.clip(cm, eps, None) / max(ca, eps)))
    return p


def run():
    scene, inserts = conrad_phantom.build_phantom()
    geo = conrad_ct.fan_geometry(n_pix=512)
    base, geo = cp.project_base_materials(scene, geo)
    acc = polychromatic_accumulators(base)              # once (90 kVp standard)
    spions = [i for i in inserts if i["c_form"] is not None and i["name"] != "SPION_c0"]

    rows = []
    for detector in DETECTORS.types:                    # EID, PCD
        for bh in (False, True):
            signal = {i["name"]: [] for i in inserts}
            for seed in range(EVAL.noise_realizations):
                sino = line_integral(acc, detector, seed)
                recon = conrad_ct.fbp(sino, geo, bh_correction=bh)
                meas = cp.measure_inserts(recon, geo, inserts)
                for m in meas:
                    signal[m["name"]].append((m["iron_delta_hu"], m["delta_hu"]))
            for ins in spions:
                arr = np.array(signal[ins["name"]])         # (reps, 2)
                d_hu = float(arr[:, 0].mean())              # c0-corrected signal
                noise = float(arr[:, 1].std())              # quantum noise on the insert
                cnr = d_hu / (noise + 1e-9)
                rows.append(dict(detector=detector, bh=bh, name=ins["name"],
                                 c_fe=ins["c_fe"], delta_hu=d_hu, noise=noise, cnr=cnr))
            print(f"[{detector} BH={int(bh)}] "
                  + "  ".join(f"{r['c_fe']:.2f}:CNR{r['cnr']:.1f}" for r in rows[-len(spions):]))
    return rows, inserts


def thresholds(rows, rose=EVAL.rose_cnr_threshold):
    """Lowest detectable c_Fe (CNR >= Rose) per (detector, BH)."""
    out = {}
    for det in DETECTORS.types:
        for bh in (False, True):
            cells = sorted([r for r in rows if r["detector"] == det and r["bh"] == bh],
                           key=lambda r: r["c_fe"])
            det_c = next((r["c_fe"] for r in cells if r["cnr"] >= rose), None)
            out[(det, bh)] = det_c
    return out


if __name__ == "__main__":
    import time
    t = time.time()
    rows, inserts = run()
    thr = thresholds(rows)
    print(f"\n[factorial] {len(rows)} cells in {time.time()-t:.0f}s")
    print("Detection thresholds (lowest c_Fe with CNR>=5, mg Fe/ml):")
    for (det, bh), c in thr.items():
        print(f"  {det} BH={int(bh)}: {c if c is not None else 'none (undetectable)'}")
