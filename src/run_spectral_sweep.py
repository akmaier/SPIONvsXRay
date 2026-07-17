"""E2 -- spectrum / detector optimization via the ANALYTIC ideal-observer CNR.

For each feasible C-arm setting -- tube voltage (config.KVP_LEVELS, 70-120 kVp) x
ADDED aluminum filtration (config.FILTER_CONFIGS; the CONRAD spectrum is already
calibrated on a real C-arm with its inherent filtration, so these are *added* mm Al)
-- we score iron detectability with the ideal-observer / matched-filter model
(src/spectral.py, Tapiovaara-Wagner / Cahn):

    EID CNR  = sqrt(SNR^2_EID)          energy weighting w(E)=E
    PCD CNR  = sqrt(SNR^2_bins)         3 bins with CNR-optimal per-bin weights

This score is deterministic and bounded by the ideal observer, so it gives the
correct EID-vs-PCD relationship (~1.2x, capped at the ~1.22x ideal ceiling) and a
physically sensible optimum (the softest feasible beam -- iron contrast is
photoelectric). It replaces the earlier recon-measured per-realization CNR, which was
noisy, over-credited PCD (a ratio above the ideal ceiling), and mis-picked the filter.
The realistic recon-measured CNRs are reported by the detectability studies (E3/E4),
which run at the optimum this module writes, and at the 90 kVp baseline.

Outputs:
  results/spectral/optimum.json   {kvp, filter, filters, pcd_bin_edges_kev, baseline,
                                   cnr_table, cnr_optimum_*, pcd_bin_optimization}
  results/spectral/e2_sweep.png   ideal-observer iron CNR vs kVp per Al filter, EID vs PCD
  results/spectral_sweep/{sweep.csv,sweep.json}   full per-setting rows
"""
from __future__ import annotations
import csv
import json
import os

import numpy as np

import conrad_backend
import spectral
from config import KVP_LEVELS, FILTER_CONFIGS


def sweep():
    """Analytic ideal-observer CNR over KVP_LEVELS x FILTER_CONFIGS for EID and PCD.

    Returns rows with per-setting EID/PCD CNR, their ratio, and the re-optimized PCD
    3-bin edges + SNR^2 bookkeeping for that spectrum.
    """
    rows = []
    for fname, filt in FILTER_CONFIGS.items():
        for kvp in KVP_LEVELS:
            E, s = spectral.real_spectrum(kvp, tuple(filt))
            mm = spectral._metrics(E, s)
            eid = float(np.sqrt(mm["eid"]))
            thr, snr2 = spectral.optimize_thresholds(mm["E"], mm["Nt"], mm["c"], 3)
            pcd = float(np.sqrt(snr2))
            edges = [round(float(E[0]), 1)] + [round(float(t), 1) for t in thr] + [round(float(kvp), 1)]
            rows.append(dict(filter=fname, kvp=float(kvp), EID=eid, PCD=pcd,
                             ratio=pcd / eid, pcd_edges=edges,
                             snr2_bins=float(snr2), eid_snr2=float(mm["eid"]),
                             ideal_snr2=float(mm["ideal"])))
            print(f"[{fname:6s} {int(kvp):3d}kVp] EID={eid:.4f}  PCD={pcd:.4f}  ({pcd/eid:.2f}x)")
    return rows


def pick_optimum(rows):
    """CNR-optimal setting = the one maximizing the (matched-filter) PCD CNR."""
    return max(rows, key=lambda r: r["PCD"])


def cnr_table(rows):
    """Ideal-observer CNR keyed [filter][kvp] -> {EID, PCD}."""
    tab = {}
    for r in rows:
        tab.setdefault(r["filter"], {})[f"{int(r['kvp'])}"] = {"EID": r["EID"], "PCD": r["PCD"]}
    return tab


def save_optimum(rows):
    best = pick_optimum(rows)
    kvp, fname = best["kvp"], best["filter"]
    filters = FILTER_CONFIGS[fname]
    bin_info = dict(snr2_bins=best["snr2_bins"], eid=best["eid_snr2"], ideal=best["ideal_snr2"],
                    cnr_gain_vs_eid=best["PCD"] / best["EID"],
                    frac_ideal=float(np.sqrt(best["snr2_bins"] / best["ideal_snr2"])))
    baseline = dict(kvp=90, filter="Al2.5", filters=[["aluminium", 2.5]])
    optimum = dict(
        kvp=float(kvp), filter=fname, filters=[[m, t] for (m, t) in filters],
        detector="PCD", pcd_bin_edges_kev=best["pcd_edges"],
        metric="ideal-observer iron CNR = sqrt(SNR^2) (spectral.py)",
        cnr_optimum_eid=best["EID"], cnr_optimum_pcd=best["PCD"],
        cnr_gain_pcd_vs_eid=best["PCD"] / best["EID"],
        pcd_bin_optimization=bin_info, baseline=baseline,
        cnr_table=cnr_table(rows))
    outdir = str(conrad_backend.REPO_ROOT / "results" / "spectral")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "optimum.json"), "w") as f:
        json.dump(optimum, f, indent=2)
    return optimum, outdir


def figure(rows, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    tab = cnr_table(rows)
    best = pick_optimum(rows)
    ymax = max(r["PCD"] for r in rows) * 1.12
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
    for ax, det in zip(axes, ("EID", "PCD")):
        for fname in FILTER_CONFIGS:
            kvps = sorted(int(k) for k in tab[fname])
            cnrs = [tab[fname][str(k)][det] for k in kvps]
            ax.plot(kvps, cnrs, "o-", lw=1.8, label=fname)
        ax.set_title(f"{det} -- iron CNR vs tube voltage")
        ax.set_xlabel("tube voltage [kVp]"); ax.set_ylim(0, ymax); ax.grid(alpha=0.3)
        ax.legend(fontsize=8, title="added Al filter")
    # mark the optimum on the PCD panel
    axes[1].plot([best["kvp"]], [best["PCD"]], "o", ms=16, mfc="none", mec="#16a085", mew=2.2)
    axes[1].annotate(f"optimum {int(best['kvp'])} kVp / {best['filter']}",
                     xy=(best["kvp"], best["PCD"]),
                     xytext=(best["kvp"] + 8, best["PCD"] + 0.10 * ymax),
                     color="#16a085", fontsize=9,
                     arrowprops=dict(arrowstyle="->", color="#16a085"))
    axes[0].set_ylabel("iron CNR (ideal observer, per ray)")
    fig.suptitle("E2 spectrum / detector optimization -- ideal-observer detectability")
    fig.tight_layout()
    figpath = os.path.join(outdir, "e2_sweep.png")
    fig.savefig(figpath, dpi=140, bbox_inches="tight"); plt.close(fig)
    return figpath


def save_rows(rows):
    outdir = str(conrad_backend.REPO_ROOT / "results" / "spectral_sweep")
    os.makedirs(outdir, exist_ok=True)
    cols = ["filter", "kvp", "EID", "PCD", "ratio"]
    with open(os.path.join(outdir, "sweep.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore"); w.writeheader()
        for r in rows:
            w.writerow(r)
    with open(os.path.join(outdir, "sweep.json"), "w") as f:
        json.dump(dict(rows=rows), f, indent=2)
    return outdir


if __name__ == "__main__":
    conrad_backend.setup()
    rows = sweep()
    optimum, outdir = save_optimum(rows)
    figpath = figure(rows, outdir)
    rowdir = save_rows(rows)
    b = optimum
    print(f"\n[E2] OPTIMUM: {b['kvp']:.0f} kVp {b['filter']} (PCD) -- "
          f"EID CNR {b['cnr_optimum_eid']:.4f}, PCD CNR {b['cnr_optimum_pcd']:.4f} "
          f"({b['cnr_gain_pcd_vs_eid']:.2f}x EID, "
          f"{b['pcd_bin_optimization']['frac_ideal']*100:.0f}% of ideal)")
    print(f"[E2] PCD 3-bin edges: {b['pcd_bin_edges_kev']} keV")
    print(f"[E2] baseline: {b['baseline']['kvp']:.0f} kVp {b['baseline']['filter']}")
    print(f"[E2] wrote {outdir}/optimum.json, {figpath}, {rowdir}/sweep.*")
