"""E2 — spectrum / detector optimization for the SPION detectability study.

Sweeps the feasible C-arm tube voltage (config.KVP_LEVELS, 70-120 kVp) crossed with
inherent aluminium filtration (config.FILTER_CONFIGS, Al 1.0/2.5/5.0/8.0 mm) for both
detectors (EID and multi-bin PCD), scoring the reconstructed iron CNR at a
representative iron concentration through the rabbit body -- end-to-end through the
CONRAD pipeline (base-material sinograms are geometry-only, computed once; each
spectrum rebuilds the polychromatic accumulators, draws noise, reconstructs, and
measures the per-insert iron dHU / CNR). Beam-hardening water precorrection and the
minimal short scan are ALWAYS on, exactly as in the study.

E2 then (a) picks the CNR-optimal (kVp, Al filter), (b) re-optimizes the PCD 3-bin
energy thresholds for that winning spectrum via the analytic ideal-observer framework
(src/spectral.py), and (c) writes results/spectral/optimum.json (consumed by
run_detectability.run_study_a / run_study_b) + a summary figure e2_sweep.png.

Outputs:
  results/spectral/optimum.json   {kvp, filter, filters, pcd_bin_edges_kev, baseline, cnr_table}
  results/spectral/e2_sweep.png   CNR vs kVp per Al filter, EID vs PCD
  results/spectral_sweep/{sweep.csv,sweep.json}   full per-insert rows
"""
from __future__ import annotations
import json
import os

import numpy as np

import conrad_backend
import conrad_ct
import conrad_phantom
import conrad_project as cp
import run_detectability as rf
import spectral
import spectrum as spec
import materials
from config import DETECTORS, EVAL, KVP_LEVELS, FILTER_CONFIGS, DOSE_LEVELS

# Representative iron concentration for the score [mg Fe/ml]: mid-range of the study
# tasks (Study A ~0.8-2.5, Study B homogenized mean ~0.5-1.5), so the ranking reflects
# the concentrations we actually care about detecting.
REPRESENTATIVE_CFE = 1.0
# Score at the realistic top-end dose so the CNR ranking is at the operating point.
SCORE_N0 = DOSE_LEVELS["high"]


def _representative_insert(spions):
    """The insert whose iron c_fe is closest to REPRESENTATIVE_CFE."""
    return min(spions, key=lambda i: abs(i["c_fe"] - REPRESENTATIVE_CFE))


def sweep(n_reps=15):
    """End-to-end CNR sweep over KVP_LEVELS x FILTER_CONFIGS x {EID,PCD}.

    Returns (rows, spions, rep_name): rows carry per-insert iron dHU / noise / CNR;
    rep_name is the representative insert used to rank spectra.
    """
    scene, inserts = conrad_phantom.build_phantom(with_bone=False)
    geo = conrad_ct.fan_geometry(n_pix=512, short_scan=True)
    base, geo = cp.project_base_materials(inserts, geo)          # geometry-only, once
    spions = [i for i in inserts if i["c_form"] is not None and i["name"] != "SPION_c0"]
    rep = _representative_insert(spions)

    rows = []
    for fname, filt in FILTER_CONFIGS.items():
        for kvp in KVP_LEVELS:
            acc = rf.polychromatic_accumulators(base, kvp=kvp, filters=tuple(filt),
                                                n0=SCORE_N0)
            for detector in DETECTORS.types:
                bh_polys = rf.bh_poly_for(acc, detector)        # water precorrection ALWAYS on
                signal = {i["name"]: [] for i in inserts}
                for seed in range(n_reps):
                    sino = rf.line_integral(acc, detector, seed, bh_polys)
                    recon = conrad_ct.fbp(sino, geo)             # short scan + distance-weighted
                    for m in cp.measure_inserts(recon, geo, inserts):
                        signal[m["name"]].append((m["iron_delta_hu"], m["delta_hu"]))
                for ins in spions:
                    arr = np.array(signal[ins["name"]])
                    d_hu = float(arr[:, 0].mean())
                    noise = float(arr[:, 1].std())
                    rows.append(dict(filter=fname, kvp=float(kvp), detector=detector,
                                     name=ins["name"], c_fe=ins["c_fe"],
                                     delta_hu=d_hu, noise=noise,
                                     cnr=d_hu / (noise + 1e-9)))
            rep_cnr = {d: next(r["cnr"] for r in rows[::-1]
                               if r["filter"] == fname and r["kvp"] == float(kvp)
                               and r["detector"] == d and r["name"] == rep["name"])
                       for d in DETECTORS.types}
            print(f"[{fname:6s} {int(kvp):3d}kVp] "
                  + "  ".join(f"{d}:CNR(c{rep['c_fe']:.2f})={rep_cnr[d]:.2f}"
                              for d in DETECTORS.types))
    return rows, spions, rep["name"]


def cnr_table(rows, rep_name):
    """CNR at the representative insert, keyed [filter][kvp][detector]."""
    tab = {}
    for r in rows:
        if r["name"] != rep_name:
            continue
        tab.setdefault(r["filter"], {}).setdefault(f"{int(r['kvp'])}", {})[r["detector"]] = r["cnr"]
    return tab


def pick_optimum(rows, rep_name):
    """CNR-optimal (kVp, filter, detector) at the representative insert."""
    cells = [r for r in rows if r["name"] == rep_name]
    best = max(cells, key=lambda r: r["cnr"])
    return best


def optimize_pcd_bins(kvp, filters):
    """Re-optimize the PCD 3-bin thresholds for a spectrum (analytic ideal observer).

    Uses src/spectral.py: transmitted background photons N_t(E) through the rabbit and
    per-energy iron line-integral contrast c(E), then brute-force maximizes the binned
    matched-filter SNR^2. Returns the 4 bin EDGES [10, t1, t2, kvp] (keV), matching
    DETECTORS.pcd_bin_edges_kev / acc["edges"] layout.
    """
    E, s = spectral.real_spectrum(kvp, tuple(filters))
    mm = spectral._metrics(E, s)
    thr3, snr2 = spectral.optimize_thresholds(mm["E"], mm["Nt"], mm["c"], 3)
    lo = float(E[0])
    hi = float(kvp)
    edges = [round(lo, 1)] + [round(float(t), 1) for t in thr3] + [round(hi, 1)]
    return edges, dict(snr2_bins=float(snr2), eid=float(mm["eid"]),
                       ideal=float(mm["ideal"]),
                       cnr_gain_vs_eid=float(np.sqrt(snr2 / mm["eid"])),
                       frac_ideal=float(np.sqrt(snr2 / mm["ideal"])))


def save_optimum(rows, rep_name, spions):
    best = pick_optimum(rows, rep_name)
    kvp, fname = best["kvp"], best["filter"]
    filters = FILTER_CONFIGS[fname]
    edges, bin_info = optimize_pcd_bins(kvp, filters)
    baseline = dict(kvp=90, filter="Al2.5", filters=[["aluminium", 2.5]])
    optimum = dict(
        kvp=float(kvp), filter=fname, filters=[[m, t] for (m, t) in filters],
        detector=best["detector"], pcd_bin_edges_kev=edges,
        representative_cfe=REPRESENTATIVE_CFE, representative_insert=rep_name,
        score_n0=SCORE_N0, cnr_optimum=float(best["cnr"]),
        pcd_bin_optimization=bin_info, baseline=baseline,
        cnr_table=cnr_table(rows, rep_name))
    outdir = str(conrad_backend.REPO_ROOT / "results" / "spectral")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "optimum.json"), "w") as f:
        json.dump(optimum, f, indent=2)
    return optimum, outdir


def figure(rows, rep_name, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    tab = cnr_table(rows, rep_name)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
    for ax, det in zip(axes, DETECTORS.types):
        for fname in FILTER_CONFIGS:
            kvps = sorted(int(k) for k in tab[fname])
            cnrs = [tab[fname][str(k)][det] for k in kvps]
            ax.plot(kvps, cnrs, "o-", lw=1.8, label=fname)
        ax.set_title(f"{det} — iron CNR vs tube voltage")
        ax.set_xlabel("tube voltage [kVp]")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, title="Al filter")
    axes[0].set_ylabel(f"iron CNR @ c_Fe={REPRESENTATIVE_CFE:g} mg/ml (rabbit body)")
    fig.suptitle("E2 spectrum / detector optimization (short scan + BH precorrection)")
    fig.tight_layout()
    figpath = os.path.join(outdir, "e2_sweep.png")
    fig.savefig(figpath, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return figpath


def save_rows(rows):
    import csv
    outdir = str(conrad_backend.REPO_ROOT / "results" / "spectral_sweep")
    os.makedirs(outdir, exist_ok=True)
    cols = ["filter", "kvp", "detector", "name", "c_fe", "delta_hu", "noise", "cnr"]
    with open(os.path.join(outdir, "sweep.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in cols})
    with open(os.path.join(outdir, "sweep.json"), "w") as f:
        json.dump(dict(rows=rows), f, indent=2)
    return outdir


if __name__ == "__main__":
    import time
    t = time.time()
    rows, spions, rep_name = sweep()
    optimum, outdir = save_optimum(rows, rep_name, spions)
    figpath = figure(rows, rep_name, outdir)
    rowdir = save_rows(rows)
    print(f"\n[E2] {len(rows)} rows in {time.time()-t:.0f}s")
    print(f"[E2] OPTIMUM: {optimum['kvp']:.0f} kVp {optimum['filter']} "
          f"({optimum['detector']}), CNR={optimum['cnr_optimum']:.2f} "
          f"at c_Fe={optimum['representative_cfe']:g} mg/ml")
    print(f"[E2] PCD 3-bin edges (re-optimized): {optimum['pcd_bin_edges_kev']} keV "
          f"(CNR {optimum['pcd_bin_optimization']['cnr_gain_vs_eid']:.2f}x EID, "
          f"{optimum['pcd_bin_optimization']['frac_ideal']*100:.0f}% of ideal)")
    print(f"[E2] baseline: {optimum['baseline']}")
    print(f"[E2] wrote {outdir}/optimum.json, {figpath}, {rowdir}/sweep.*")
