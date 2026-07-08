"""Filter / tube-voltage sweep through the full CONRAD pipeline.

The M6 factorial fixes the spectrum at 90 kVp. Iron contrast is photoelectric, so
it should improve at lower kVp and degrade under hardening filters (spectral.py
predicts this for the ideal observer). Here we confirm it end-to-end: for each
spectrum we rebuild the polychromatic accumulators (base-material sinograms are
geometry-only, computed once), draw noise, reconstruct, and measure the per-insert
iron CNR + detection threshold for EID and PCD (beam-hardening off).

Outputs results/spectral_sweep/{sweep.csv,sweep.json}.
"""
from __future__ import annotations
import numpy as np

import conrad_backend
import conrad_ct
import conrad_phantom
import conrad_project as cp
import run_factorial as rf
from config import DETECTORS, EVAL, KVP_LEVELS, FILTER_CONFIGS


def sweep(n_reps=15):
    scene, inserts = conrad_phantom.build_phantom()
    geo = conrad_ct.fan_geometry(n_pix=512)
    base, geo = cp.project_base_materials(scene, geo)          # geometry-only, once
    spions = [i for i in inserts if i["c_form"] is not None and i["name"] != "SPION_c0"]

    # spectra: kVp variants (no added filter) + filter configs at the 90 kVp base
    configs = [(f"{int(k)}kVp", k, ()) for k in KVP_LEVELS]
    configs += [(name, None, tuple(filt)) for name, filt in FILTER_CONFIGS.items()]

    rows = []
    for label, kvp, filters in configs:
        acc = rf.polychromatic_accumulators(base, kvp=kvp, filters=filters)
        for detector in DETECTORS.types:
            signal = {i["name"]: [] for i in inserts}
            for seed in range(n_reps):
                sino = rf.line_integral(acc, detector, seed)
                recon = conrad_ct.fbp(sino, geo, bh_correction=False)
                for m in cp.measure_inserts(recon, geo, inserts):
                    signal[m["name"]].append((m["iron_delta_hu"], m["delta_hu"]))
            for ins in spions:
                arr = np.array(signal[ins["name"]])
                d_hu = float(arr[:, 0].mean())
                noise = float(arr[:, 1].std())
                rows.append(dict(spectrum=label, kvp=(kvp or 90.0), detector=detector,
                                 name=ins["name"], c_fe=ins["c_fe"],
                                 delta_hu=d_hu, noise=noise, cnr=d_hu / (noise + 1e-9)))
        # threshold at the realistic top-end dose readout
        cnr_top = {d: next(r["cnr"] for r in rows[::-1]
                           if r["spectrum"] == label and r["detector"] == d) for d in DETECTORS.types}
        print(f"[{label:16s}] " + "  ".join(f"{d}:CNR(c1.09)={cnr_top[d]:.1f}" for d in DETECTORS.types))
    return rows, spions


def thresholds(rows, rose=EVAL.rose_cnr_threshold):
    out = {}
    labels = sorted({r["spectrum"] for r in rows})
    for lab in labels:
        for det in DETECTORS.types:
            cells = sorted([r for r in rows if r["spectrum"] == lab and r["detector"] == det],
                           key=lambda r: r["c_fe"])
            out[f"{lab}_{det}"] = next((r["c_fe"] for r in cells if r["cnr"] >= rose), None)
    return out


if __name__ == "__main__":
    import csv, json, os, time
    t = time.time()
    rows, spions = sweep()
    thr = thresholds(rows)
    outdir = str(conrad_backend.REPO_ROOT / "results" / "spectral_sweep")
    os.makedirs(outdir, exist_ok=True)
    cols = ["spectrum", "kvp", "detector", "name", "c_fe", "delta_hu", "noise", "cnr"]
    with open(os.path.join(outdir, "sweep.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in cols})
    with open(os.path.join(outdir, "sweep.json"), "w") as f:
        json.dump(dict(rows=rows, thresholds_rose5=thr,
                       thresholds_rose3=thresholds(rows, rose=3.0)), f, indent=2)
    print(f"[sweep] {len(rows)} rows in {time.time()-t:.0f}s -> {outdir}")
    print("Detection thresholds (lowest c_Fe with CNR>=5):")
    for k, v in thr.items():
        print(f"  {k}: {v if v is not None else 'none'}")
