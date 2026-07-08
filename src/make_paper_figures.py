"""Generate publication figures for the SPIE manuscript from the persisted results.

Reads results/factorial/factorial.json (+ results/spectral_sweep/sweep.json if
present) and writes vector PDFs into paper/figures/. Kept separate from the
simulation so figures can be regenerated without re-running CONRAD.
"""
from __future__ import annotations
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import conrad_backend
from config import tumor_iron_conc, DELIVERED_AT_10_MG

REPO = conrad_backend.REPO_ROOT
FIGDIR = str(REPO / "paper" / "figures")
os.makedirs(FIGDIR, exist_ok=True)

C_REALISTIC = tumor_iron_conc(10.0)   # mg Fe/ml at the 6 mg delivered dose

plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3,
                     "figure.dpi": 150, "savefig.bbox": "tight"})
EID_C, PCD_C = "#2c3e50", "#c0392b"


def _load(path):
    with open(path) as f:
        return json.load(f)


def fig_dhu(fac):
    rows = [r for r in fac["rows"] if not r["bh"]]
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    for det, c, mk in (("EID", EID_C, "o"), ("PCD", PCD_C, "s")):
        rr = sorted([r for r in rows if r["detector"] == det], key=lambda r: r["c_fe"])
        ax.plot([r["c_fe"] for r in rr], [r["delta_hu"] for r in rr], mk + "-",
                color=c, label=det, ms=5)
    ax.axvline(C_REALISTIC, color="#888", ls=":", lw=1)
    ax.text(C_REALISTIC, ax.get_ylim()[1] * 0.96, " 6 mg dose", color="#555",
            fontsize=8, va="top")
    ax.set_xlabel("tumor iron concentration $c_{\\mathrm{Fe}}$ [mg/ml]")
    ax.set_ylabel("iron contrast $\\Delta$HU")
    ax.set_title("Iron contrast vs. dose (beam-hardening off)")
    ax.legend(frameon=False)
    fig.savefig(f"{FIGDIR}/fig_dhu_vs_conc.pdf"); plt.close(fig)


def fig_cnr(fac):
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    styles = {("EID", False): (EID_C, "o", "-", "EID, BH off"),
              ("EID", True): (EID_C, "o", "--", "EID, BH on"),
              ("PCD", False): (PCD_C, "s", "-", "PCD, BH off"),
              ("PCD", True): (PCD_C, "s", "--", "PCD, BH on")}
    for (det, bh), (c, mk, ls, lab) in styles.items():
        rr = sorted([r for r in fac["rows"] if r["detector"] == det and r["bh"] == bh],
                    key=lambda r: r["c_fe"])
        ax.plot([r["c_fe"] for r in rr], [r["cnr"] for r in rr], mk + ls,
                color=c, label=lab, ms=4, lw=1.3)
    ax.axhline(5, color="#27ae60", lw=1); ax.text(ax.get_xlim()[1], 5, " Rose 5",
                                                  color="#27ae60", fontsize=8, va="bottom", ha="right")
    ax.axhline(3, color="#27ae60", lw=1, ls=":"); ax.text(ax.get_xlim()[1], 3, " Rose 3",
                                                          color="#27ae60", fontsize=8, va="bottom", ha="right")
    ax.axvline(C_REALISTIC, color="#888", ls=":", lw=1)
    ax.set_xlabel("tumor iron concentration $c_{\\mathrm{Fe}}$ [mg/ml]")
    ax.set_ylabel("CNR")
    ax.set_title("Detectability vs. dose")
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(f"{FIGDIR}/fig_cnr_vs_conc.pdf"); plt.close(fig)


def fig_sweep(sw):
    rows = sw["rows"]
    top = max(r["c_fe"] for r in rows)
    labels = []
    seen = set()
    for r in rows:                       # preserve encounter order
        if r["spectrum"] not in seen:
            seen.add(r["spectrum"]); labels.append(r["spectrum"])
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    x = np.arange(len(labels)); wdt = 0.38
    for i, (det, c) in enumerate((("EID", EID_C), ("PCD", PCD_C))):
        vals = []
        for lab in labels:
            m = [r["cnr"] for r in rows if r["spectrum"] == lab and r["detector"] == det
                 and abs(r["c_fe"] - top) < 1e-6]
            vals.append(m[0] if m else 0.0)
        ax.bar(x + (i - 0.5) * wdt, vals, wdt, color=c, label=det)
    ax.axhline(5, color="#27ae60", lw=1, ls="--")
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel(f"CNR at $c_{{\\mathrm{{Fe}}}}$={top:.2f} mg/ml")
    ax.set_title("Spectral shaping: tube voltage and filtration")
    ax.legend(frameon=False)
    fig.savefig(f"{FIGDIR}/fig_spectral_sweep.pdf"); plt.close(fig)


def main():
    fac = _load(str(REPO / "results" / "factorial" / "factorial.json"))
    fig_dhu(fac); fig_cnr(fac)
    print("[figs] fig_dhu_vs_conc.pdf, fig_cnr_vs_conc.pdf")
    sw_path = REPO / "results" / "spectral_sweep" / "sweep.json"
    if sw_path.exists():
        fig_sweep(_load(str(sw_path)))
        print("[figs] fig_spectral_sweep.pdf")
    print("[figs] wrote", FIGDIR)


if __name__ == "__main__":
    main()
