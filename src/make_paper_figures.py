"""Generate publication figures + LaTeX tables for the SPION-vs-X-ray manuscript.

Reads the COMPLETED results on disk and regenerates every paper figure and table:

  Figures  (paper/figures/*.pdf  [vector]  +  docs/assets/*.png  [dashboard])
    fig_physics          mu(E) iron-oxide vs soft tissue + 70/90 kVp spectra; mono dHU
    fig_phantom_recon    Study A phantom recon, EID vs PCD, noise-free + 1 realization
    fig_studyA_cnr       CNR vs tumor c_Fe, EID vs PCD (optimum), density range, Rose 5
    fig_density_montage  loading x density {1e8,3e8,1e9} recon montage, EID + PCD
    fig_studyB           CNR vs vessel level, EID vs PCD, optimum vs baseline, high dose
    fig_spectral         analytic ideal-observer iron CNR vs kVp per Al filter, EID vs
                         PCD (shared y), optimum marked on PCD (read from optimum.json)
    fig_eid_vs_pcd       PCD/EID CNR ratio across configs, optimum vs baseline

  Tables  (paper/tables/*.tex, booktabs)
    tab_composition      magnetite Fe3O4 + PAA, phi per formulation, c_Fe->c_NP example
    tab_loading_density  8 loading configs x {1e8,3e8,1e9} c_Fe grid
    tab_thresholds       Rose CNR>=5 lowest detectable c_Fe per detector/density/dose
    tab_cnr              PCD/EID CNR per config, optimum & baseline

IMPORTANT: only bone=False (with_bone False) cells are used for detectability -- the
bone rod corrupts the c0 reference (CNR 40-600) and is a known limitation, not a result.

The physics/recon figures (fig_physics, fig_phantom_recon, fig_density_montage) require
a live CONRAD (pyconrad) session and MUST run from the src/ directory. The data-driven
figures and all tables read only the persisted JSON and always regenerate.
"""
from __future__ import annotations
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np

import conrad_backend
from config import (tumor_iron_conc, FE_FRACTION, PAA_MASS_FRAC_BY_FORMULATION,
                    CELLULAR_LOADING, CELL_DENSITY_LEVELS, cfe_from_loading,
                    nanoparticle_conc_from_fe, STUDY_B_VESSEL_LEVELS,
                    VESSEL_VOLUME_FRACTION, DOSE_LEVELS)

REPO = conrad_backend.REPO_ROOT
FIGDIR = str(REPO / "paper" / "figures")
DASHDIR = str(REPO / "docs" / "assets")
TABDIR = str(REPO / "paper" / "tables")
DET = str(REPO / "results" / "detectability")
for d in (FIGDIR, DASHDIR, TABDIR):
    os.makedirs(d, exist_ok=True)

# ----------------------------------------------------------------------------
# Global style: clean, legible, colorblind-safe (Okabe-Ito derived).
# ----------------------------------------------------------------------------
plt.rcParams.update({
    # Base sizes raised ~35% for in-figure legibility (reviewer request).
    "font.size": 13, "font.family": "sans-serif",
    "axes.labelsize": 14, "axes.titlesize": 14,
    "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 12,
    "axes.grid": True, "grid.alpha": 0.3, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150, "savefig.bbox": "tight", "savefig.dpi": 150,
    "legend.frameon": False, "axes.linewidth": 0.8,
})
EID_C = "#0072B2"      # blue   (Okabe-Ito) - energy-integrating
PCD_C = "#D55E00"      # vermilion            - photon-counting
OPT_C = "#009E73"      # green  - optimum spectrum
BASE_C = "#999999"     # grey   - baseline spectrum
ROSE_C = "#CC79A7"     # reddish-purple - Rose threshold


def _save(fig, name, png=True):
    fig.savefig(f"{FIGDIR}/{name}.pdf")
    if png:
        fig.savefig(f"{DASHDIR}/{name}.png", dpi=140)
    plt.close(fig)
    print(f"[fig] {name}.pdf" + (f" + {name}.png" if png else ""))


def _load(path):
    with open(path) as f:
        return json.load(f)


def _bone_off(rows):
    """Only bone=False cells (bone=True is degenerate, see module docstring)."""
    return [r for r in rows if not r["with_bone"]]


def _spec_label(kvp, filters):
    """Human-readable spectrum label, e.g. '70 kVp / Al 1.0'."""
    mat, mm = filters[0]
    sym = {"aluminium": "Al", "aluminum": "Al", "copper": "Cu"}.get(mat, mat)
    return f"{float(kvp):.0f} kVp / {sym} {float(mm):.1f}"


# ============================================================================
# FIGURE 1 -- physics: mu(E) + spectra, and mono dHU vs c_Fe
# ============================================================================
def fig_physics():
    import materials
    import spectrum as spec

    E = np.arange(20.0, 120.01, 0.5)
    mu_soft = materials.linear_attenuation_soft(E)
    ox = materials.oxide_contrast_massatten(E)                 # cm^2/g per g Fe
    # tumor mu at a representative iron load (1 mg Fe/ml == 1e-3 g/cm^3)
    mu_tumor1 = mu_soft + 1e-3 * ox
    mu_tumor10 = mu_soft + 10e-3 * ox

    opt = _load(str(REPO / "results" / "spectral" / "optimum.json"))
    opt_kvp = float(opt["kvp"]); opt_filters = [tuple(x) for x in opt["filters"]]
    bas = opt["baseline"]
    bas_kvp = float(bas["kvp"]); bas_filters = [tuple(x) for x in bas["filters"]]
    opt_lab = _spec_label(opt_kvp, opt_filters) + " (optimum)"
    bas_lab = _spec_label(bas_kvp, bas_filters) + " (baseline)"

    E70, f70 = spec.conrad_spectrum(opt_kvp)
    f70 = spec.apply_filters(E70, f70, opt_filters)            # current optimum
    E90, f90 = spec.conrad_spectrum(bas_kvp)
    f90 = spec.apply_filters(E90, f90, bas_filters)            # baseline
    f70n = f70 / f70.max()
    f90n = f90 / f90.max()

    fig = plt.figure(figsize=(10.4, 3.7))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.25, 1.0], wspace=0.52)

    # (a) attenuation + spectra
    ax = fig.add_subplot(gs[0])
    ax.plot(E, mu_soft, color="#444", lw=2, label="soft tissue")
    ax.plot(E, mu_tumor1, color=PCD_C, lw=1.8, label=r"tumor +1 mg Fe/ml")
    ax.plot(E, mu_tumor10, color=PCD_C, lw=1.8, ls="--", label=r"tumor +10 mg Fe/ml")
    ax.set_xlabel("energy [keV]")
    ax.set_ylabel(r"linear attenuation $\mu$ [cm$^{-1}$]")
    ax.set_yscale("log")
    ax.set_xlim(20, 120)
    ax.legend(loc="upper right", fontsize=11)
    ax.set_title("(a) Iron-oxide contrast vs. spectra", fontsize=13, loc="left")

    axr = ax.twinx()
    axr.grid(False)
    axr.fill_between(E70, 0, f70n, color=OPT_C, alpha=0.18)
    axr.plot(E70, f70n, color=OPT_C, lw=1.2, label=opt_lab)
    axr.plot(E90, f90n, color=BASE_C, lw=1.2, ls="-", label=bas_lab)
    axr.set_ylabel("normalized photon flux", fontsize=12)
    axr.set_ylim(0, 1.05)
    axr.spines["top"].set_visible(False)
    axr.legend(loc="center right", fontsize=10, bbox_to_anchor=(1.0, 0.62))

    # (b) mono dHU vs c_Fe at 3 mono energies
    axb = fig.add_subplot(gs[1])
    cfe = np.linspace(0, 10, 60)
    for Emono, col in ((40.0, "#56B4E9"), (60.0, EID_C), (80.0, "#004C6D")):
        mu_w = float(materials.linear_attenuation_soft(np.array([Emono]))[0])
        oxm = float(materials.oxide_contrast_massatten(np.array([Emono]))[0])
        dHU = 1000.0 * (1e-3 * cfe * oxm) / mu_w
        axb.plot(cfe, dHU, color=col, lw=1.8, label=f"{Emono:.0f} keV")
    axb.set_xlabel(r"tumor iron $c_{\mathrm{Fe}}$ [mg/ml]")
    axb.set_ylabel(r"mono. iron contrast $\Delta$HU")
    axb.legend(title="mono energy", fontsize=11, title_fontsize=11)
    axb.set_title("(b) Monoenergetic iron contrast", fontsize=13, loc="left")

    _save(fig, "fig_physics")
    # headline: dHU at 60 keV for 1 mg Fe/ml
    mu_w60 = float(materials.linear_attenuation_soft(np.array([60.0]))[0])
    ox60 = float(materials.oxide_contrast_massatten(np.array([60.0]))[0])
    return 1000.0 * (1e-3 * 1.0 * ox60) / mu_w60


# ============================================================================
# Recon helpers (live CONRAD)
# ============================================================================
def _noise_free_sino(acc, det, rf):
    if det == "EID":
        return -np.log(np.clip(acc["S_det"], 1e-6, None) / acc["S_air"])
    w = rf._pcd_weights(acc)
    M = np.zeros_like(acc["S_det"]); Ma = 0.0
    for b, (cd, ca) in enumerate(zip(acc["C_det"], acc["C_air"])):
        M += w[b] * cd; Ma += w[b] * ca
    return -np.log(np.clip(M, 1e-6, None) / max(Ma, 1e-6))


def _to_hu(img, sp, roi_mm=8.0):
    N = img.shape[0]
    yy, xx = np.mgrid[0:N, 0:N]
    c = N / 2.0
    bg = ((xx - c) ** 2 + (yy - c) ** 2) < (roi_mm / sp) ** 2
    mu_w = float(np.median(img[bg]))
    return 1000.0 * (img - mu_w) / mu_w, c


# ============================================================================
# FIGURE 2 -- Study A phantom recon: EID vs PCD, noise-free + one realization
# ============================================================================
def fig_phantom_recon():
    import conrad_ct
    import conrad_phantom
    import conrad_project as cp
    import run_detectability as rf

    opt = _load(str(REPO / "results" / "spectral" / "optimum.json"))
    kvp = float(opt["kvp"]); filters = [tuple(x) for x in opt["filters"]]
    edges = tuple(opt["pcd_bin_edges_kev"])

    # Study A loading phantom at the mid density, high dose, bone omitted for display.
    scene, inserts = conrad_phantom.build_phantom(loading=True,
                                                  density=CELL_DENSITY_LEVELS["1e9"],
                                                  with_bone=False)
    geo = conrad_ct.fan_geometry(n_pix=512, short_scan=True)
    base, geo = cp.project_base_materials(inserts, geo)
    sp = geo["voxel_mm"]
    n0 = DOSE_LEVELS["high"]
    acc = rf.polychromatic_accumulators(base, kvp=kvp, filters=filters, n0=n0,
                                        bin_edges=edges)

    panels = {}
    for det in ("EID", "PCD"):
        bh = rf.bh_poly_for(acc, det)
        nf = _noise_free_sino(acc, det, rf)
        nf = np.polyval(bh, nf)
        panels[(det, "noise-free")] = conrad_ct.fbp(nf, geo)
        panels[(det, "1 realization")] = conrad_ct.fbp(
            rf.line_integral(acc, det, seed=7, bh_polys=bh), geo)

    WIN = 20
    fig, ax = plt.subplots(2, 2, figsize=(7.4, 7.7))
    order = [("EID", "noise-free"), ("PCD", "noise-free"),
             ("EID", "1 realization"), ("PCD", "1 realization")]
    for a, key in zip(ax.ravel(), order):
        hu, c = _to_hu(panels[key], sp)
        im = a.imshow(hu, cmap="gray", vmin=-WIN, vmax=WIN)
        a.set_title(f"{key[0]} - {key[1]}", fontsize=14)
        a.axis("off")
    # label the highest-iron insert on the top-left panel
    hu, c = _to_hu(panels[("EID", "noise-free")], sp)
    for ins in inserts:
        if ins.get("c_fe", 0) and ins["c_fe"] > 0:
            cx, cy = ins["center_mm"]
            ax[0, 0].text(cx / sp + c, cy / sp + c + 16, f"{ins['c_fe']:.1f}",
                          color="#ffd23b", fontsize=9, ha="center", va="center")
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02, shrink=0.85)
    cbar.set_label("HU (iron window)", fontsize=12)
    fig.suptitle("Study A phantom (cellular loading, density $10^9$/cm$^3$, high dose)\n"
                 f"{_spec_label(kvp, filters)} optimum  ·  window [$-${WIN}, {WIN}] HU  ·  bone omitted for display",
                 fontsize=13, y=0.99)
    _save(fig, "fig_phantom_recon")


# ============================================================================
# FIGURE 3 -- Study A: CNR vs tumor c_Fe, EID vs PCD (optimum), density range
# ============================================================================
def fig_studyA_cnr():
    a = _load(f"{DET}/study_a_optimum.json")
    rows = _bone_off(a["rows"])
    sp = a["spectrum"]; opt_lab = _spec_label(sp["kvp"], sp["filters"])
    fig, ax = plt.subplots(figsize=(5.4, 3.9))

    for det, col, mk in (("EID", EID_C, "o"), ("PCD", PCD_C, "s")):
        rr = sorted([r for r in rows if r["detector"] == det and r["dose"] == "high"],
                    key=lambda r: r["c_fe"])
        ax.plot([r["c_fe"] for r in rr], [r["cnr"] for r in rr], mk,
                color=col, ms=5, alpha=0.85, label=f"{det} (high dose)")
        rrl = sorted([r for r in rows if r["detector"] == det and r["dose"] == "low"],
                     key=lambda r: r["c_fe"])
        ax.plot([r["c_fe"] for r in rrl], [r["cnr"] for r in rrl], mk,
                color=col, ms=4, alpha=0.35, mfc="none", label=f"{det} (low dose)")

    ax.axhline(5, color=ROSE_C, lw=1.4, ls="--")
    ax.text(ax.get_xlim()[1], 5, " Rose CNR = 5", color=ROSE_C, fontsize=11,
            va="bottom", ha="right")

    # SPION I fresh (I_113_0h) marker: its c_Fe across the density range
    fresh = [c for (lab, form, pg0, pg24) in CELLULAR_LOADING if lab == "I_113"
             for c in [pg0]][0]
    cfe_lo = cfe_from_loading(fresh, CELL_DENSITY_LEVELS["1e8"])
    cfe_hi = cfe_from_loading(fresh, CELL_DENSITY_LEVELS["1e9"])
    ax.axvspan(cfe_lo, cfe_hi, color="#f0e442", alpha=0.18, lw=0,
               label="SPION I fresh, $10^{8}$-$10^{9}$/cm$^3$")
    ax.axvline(cfe_from_loading(fresh, CELL_DENSITY_LEVELS["3e8"]), color="#B8860B",
               lw=1.0, ls=":")

    ax.set_xlabel(r"tumor iron $c_{\mathrm{Fe}}$ [mg/ml]")
    ax.set_ylabel("CNR")
    ax.set_title(f"Study A - cellular loading ({opt_lab} optimum)", fontsize=13)
    ax.legend(fontsize=10, loc="upper left")
    ax.set_ylim(bottom=0)
    _save(fig, "fig_studyA_cnr")

    # headline: max PCD CNR at high dose
    pcdhi = [r["cnr"] for r in rows if r["detector"] == "PCD" and r["dose"] == "high"]
    return max(pcdhi)


# ============================================================================
# FIGURE 4 -- density montage: loading x density recon, EID + PCD, noise-free
# ============================================================================
def fig_density_montage():
    import conrad_ct
    import conrad_phantom
    import conrad_project as cp
    import run_detectability as rf

    opt = _load(str(REPO / "results" / "spectral" / "optimum.json"))
    kvp = float(opt["kvp"]); filters = [tuple(x) for x in opt["filters"]]
    edges = tuple(opt["pcd_bin_edges_kev"])
    dens_keys = ["1e8", "3e8", "1e9"]
    n0 = DOSE_LEVELS["high"]

    grid = {}
    sp = None
    inserts_ref = None
    for dk in dens_keys:
        scene, inserts = conrad_phantom.build_phantom(loading=True,
                                                      density=CELL_DENSITY_LEVELS[dk],
                                                      with_bone=False)
        geo = conrad_ct.fan_geometry(n_pix=512, short_scan=True)
        base, geo = cp.project_base_materials(inserts, geo)
        sp = geo["voxel_mm"]; inserts_ref = inserts
        acc = rf.polychromatic_accumulators(base, kvp=kvp, filters=filters, n0=n0,
                                            bin_edges=edges)
        for det in ("EID", "PCD"):
            bh = rf.bh_poly_for(acc, det)
            nf = np.polyval(bh, _noise_free_sino(acc, det, rf))
            grid[(det, dk)] = conrad_ct.fbp(nf, geo)

    WIN = 20
    fig, ax = plt.subplots(2, 3, figsize=(9.6, 6.6))
    dens_lab = {"1e8": r"$10^{8}$/cm$^3$", "3e8": r"$3{\times}10^{8}$/cm$^3$",
                "1e9": r"$10^{9}$/cm$^3$"}
    for i, det in enumerate(("EID", "PCD")):
        for j, dk in enumerate(dens_keys):
            hu, c = _to_hu(grid[(det, dk)], sp)
            im = ax[i, j].imshow(hu, cmap="gray", vmin=-WIN, vmax=WIN)
            ax[i, j].axis("off")
            if i == 0:
                ax[i, j].set_title(dens_lab[dk], fontsize=14)
            if j == 0:
                ax[i, j].text(-0.06, 0.5, det, transform=ax[i, j].transAxes,
                              rotation=90, va="center", ha="center", fontsize=15,
                              fontweight="bold")
    cbar = fig.colorbar(im, ax=ax, fraction=0.022, pad=0.02, shrink=0.8)
    cbar.set_label("HU (iron window)", fontsize=12)
    fig.suptitle(f"Study A recon vs. tumor cell density (noise-free, {_spec_label(kvp, filters)})\n"
                 f"iron scales with density  ·  window [$-${WIN}, {WIN}] HU  ·  bone omitted for display",
                 fontsize=13, y=0.99)
    _save(fig, "fig_density_montage")


# ============================================================================
# FIGURE 5 -- Study B: CNR vs vessel level, EID vs PCD, optimum vs baseline
# ============================================================================
def fig_studyB():
    opt = _load(f"{DET}/study_b_optimum.json")
    base = _load(f"{DET}/study_b_baseline.json")
    fig, ax = plt.subplots(figsize=(5.4, 3.9))

    for data, sname, ls, alpha in ((opt, "optimum", "-", 1.0),
                                    (base, "baseline", "--", 0.6)):
        rows = _bone_off(data["rows"])
        for det, col, mk in (("EID", EID_C, "o"), ("PCD", PCD_C, "s")):
            rr = sorted([r for r in rows if r["detector"] == det and r["dose"] == "high"],
                        key=lambda r: r["vessel_level"])
            ax.plot([r["vessel_level"] for r in rr], [r["cnr"] for r in rr],
                    mk + ls, color=col, ms=5, alpha=alpha,
                    label=f"{det} ({sname})")

    ax.axhline(5, color=ROSE_C, lw=1.4, ls="--")
    ax.text(ax.get_xlim()[1], 5, " Rose CNR = 5", color=ROSE_C, fontsize=11,
            va="bottom", ha="right")
    ax.set_xlabel("vessel-local iron concentration [mg Fe/ml]")
    ax.set_ylabel("CNR")
    ax.set_title("Study B - fresh vascular injection (high dose)", fontsize=13)
    ax.legend(fontsize=10, ncol=2, loc="upper left")
    ax.set_ylim(bottom=0)
    _save(fig, "fig_studyB")

    pcd = [r["cnr"] for r in _bone_off(opt["rows"])
           if r["detector"] == "PCD" and r["dose"] == "high"]
    return max(pcd)


# ============================================================================
# FIGURE 6 -- spectral sweep: CNR vs kVp per Al filter, EID vs PCD
# ============================================================================
def fig_spectral():
    """Analytic ideal-observer iron CNR, EID vs PCD, per added-Al filter vs kVp.

    Reads the corrected ideal-observer CNR table (optimum.json / sweep.csv), NOT
    the old recon-measured CNR. Two panels share a y-axis so the ~1.2x PCD gain
    at the 70 kVp / Al 1.0 optimum is read directly off the plot.
    """
    opt = _load(str(REPO / "results" / "spectral" / "optimum.json"))
    table = opt["cnr_table"]
    opt_kvp = float(opt["kvp"]); opt_filter = opt["filter"]      # e.g. 70, "Al1.0"
    opt_lab = _spec_label(opt_kvp, opt["filters"])              # "70 kVp / Al 1.0"

    filters = ["Al1.0", "Al2.5", "Al5.0", "Al8.0"]
    kvps = sorted(float(k) for k in table[filters[0]].keys())

    def cnr(fl, kvp, det):
        return table[fl][str(int(kvp))][det]

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.9), sharey=True)
    cmap = plt.cm.viridis(np.linspace(0.15, 0.85, len(filters)))
    allv = [cnr(fl, k, det) for fl in filters for k in kvps for det in ("EID", "PCD")]
    ymax = max(allv) * 1.15
    for ax, det in zip(axes, ("EID", "PCD")):
        for fl, col in zip(filters, cmap):
            ax.plot(kvps, [cnr(fl, k, det) for k in kvps], "o-",
                    color=col, ms=5, lw=1.6, label=fl.replace("Al", "Al "))
        ax.set_xlabel("tube voltage [kVp]")
        ax.set_title(f"{det} - iron CNR vs tube voltage", fontsize=14)
        ax.set_xticks(kvps)
        ax.set_ylim(0, ymax)
    axes[0].set_ylabel("iron CNR (ideal observer)")
    axes[0].legend(title="added Al filter", fontsize=11, title_fontsize=11,
                   loc="upper right")

    # mark the current optimum (70 kVp / Al 1.0) on the PCD panel
    cval = cnr(opt_filter, opt_kvp, "PCD")
    axes[1].scatter([opt_kvp], [cval], s=150, facecolors="none",
                    edgecolors=OPT_C, linewidths=2.2, zorder=5)
    axes[1].annotate(f"optimum\n{opt_lab}", xy=(opt_kvp, cval),
                     xytext=(108, cval - 0.02 * ymax), fontsize=11, color=OPT_C,
                     ha="center", va="center",
                     arrowprops=dict(arrowstyle="->", color=OPT_C, lw=1.2))
    fig.suptitle("Spectral shaping sweep: tube voltage x Al filtration (ideal observer)",
                 fontsize=13, y=1.0)
    fig.tight_layout()
    _save(fig, "fig_spectral")
    return cval


# ============================================================================
# FIGURE 7 -- PCD/EID CNR ratio across configs, optimum vs baseline
# ============================================================================
def _ratio_by_config(data, cfgkey, dose="high"):
    rows = _bone_off(data["rows"])
    out = {}
    for r in rows:
        if r["dose"] != dose:
            continue
        # aggregate across densities (Study A has 3): average CNR per config/detector
        out.setdefault((r[cfgkey], r["detector"]), []).append(r["cnr"])
    cfgs = sorted({k[0] for k in out})
    ratios = {}
    for cfg in cfgs:
        e = out.get((cfg, "EID")); p = out.get((cfg, "PCD"))
        if e and p:
            ratios[cfg] = np.mean(p) / np.mean(e)
    return ratios


def fig_eid_vs_pcd():
    aopt = _load(f"{DET}/study_a_optimum.json")
    abas = _load(f"{DET}/study_a_baseline.json")
    bopt = _load(f"{DET}/study_b_optimum.json")
    bbas = _load(f"{DET}/study_b_baseline.json")
    opt_lab = _spec_label(aopt["spectrum"]["kvp"], aopt["spectrum"]["filters"])
    bas_lab = _spec_label(abas["spectrum"]["kvp"], abas["spectrum"]["filters"])

    ra_o = _ratio_by_config(aopt, "config")
    ra_b = _ratio_by_config(abas, "config")
    rb_o = _ratio_by_config(bopt, "config")
    rb_b = _ratio_by_config(bbas, "config")

    cfgs = list(ra_o.keys()) + list(rb_o.keys())
    x = np.arange(len(cfgs))
    opt_vals = [ra_o.get(c, rb_o.get(c)) for c in cfgs]
    bas_vals = [ra_b.get(c, rb_b.get(c)) for c in cfgs]

    fig, ax = plt.subplots(figsize=(8.4, 3.9))
    wd = 0.4
    ax.bar(x - wd / 2, opt_vals, wd, color=OPT_C, label=f"optimum ({opt_lab})")
    ax.bar(x + wd / 2, bas_vals, wd, color=BASE_C, label=f"baseline ({bas_lab})")
    ax.axhline(1.0, color="#333", lw=1.0, ls="--")
    ax.text(len(cfgs) - 0.5, 1.0, " PCD = EID", fontsize=11, va="bottom", ha="right",
            color="#333")

    # Headroom so the section labels sit in a clear band above every bar and the
    # (now larger) legend drops just below that band without colliding.
    top = max(v for v in opt_vals + bas_vals if v is not None) * 1.20
    ax.set_ylim(top=top)

    nA = len(ra_o)
    ax.axvline(nA - 0.5, color="#ccc", lw=1.0)
    ax.text((nA - 1) / 2, top * 0.98, "Study A", ha="center",
            va="top", fontsize=12, color="#555")
    ax.text((nA + len(cfgs) - 1) / 2, top * 0.98, "Study B", ha="center",
            va="top", fontsize=12, color="#555")

    ax.set_xticks(x)
    ax.set_xticklabels([c.replace("_", " ") for c in cfgs], rotation=40, ha="right",
                       fontsize=10)
    ax.set_ylabel("PCD / EID CNR ratio")
    ax.set_title("Photon-counting CNR gain over energy-integrating (bone off, high dose)",
                 fontsize=13)
    ax.legend(fontsize=11, loc="upper right", bbox_to_anchor=(1.0, 0.88))
    _save(fig, "fig_eid_vs_pcd")

    allr = [v for v in opt_vals if v is not None]
    return float(np.mean(allr))


# ============================================================================
# TABLES
# ============================================================================
def _write_tex(name, body):
    with open(f"{TABDIR}/{name}.tex", "w") as f:
        f.write(body)
    print(f"[tab] {name}.tex")


def tab_composition():
    phiI = PAA_MASS_FRAC_BY_FORMULATION["SPION_I"]
    phiII = PAA_MASS_FRAC_BY_FORMULATION["SPION_II"]
    # example: c_Fe = 1 mg/ml -> c_NP
    cnpI = nanoparticle_conc_from_fe(1.0, phiI)
    cnpII = nanoparticle_conc_from_fe(1.0, phiII)
    body = r"""\begin{tabular}{@{}llc@{}}
\toprule
Component & Description & Value \\
\midrule
Core & Magnetite Fe$_3$O$_4$ & iron fraction $F_{\mathrm{Fe}}=%.3f$ \\
Coating & Polyacrylic acid (PAA), (C$_3$H$_4$O$_2$)$_n$ & low-$Z$, tissue-equivalent \\
$\varphi$ (SPION~I) & PAA mass fraction (TGA, 12\,nm core) & %.2f \\
$\varphi$ (SPION~II) & PAA mass fraction (TGA, 8\,nm core) & %.2f \\
\midrule
\multicolumn{3}{@{}l}{Example conversion $c_{\mathrm{Fe}}=1.0$\,mg\,Fe/ml $\rightarrow$ whole-particle $c_{\mathrm{NP}}$:} \\
SPION~I  & $c_{\mathrm{NP}}=c_{\mathrm{Fe}}/[F_{\mathrm{Fe}}(1-\varphi)]$ & %.2f\,mg\,SPION/ml \\
SPION~II & & %.2f\,mg\,SPION/ml \\
\bottomrule
\end{tabular}""" % (FE_FRACTION, phiI, phiII, cnpI, cnpII)
    _write_tex("tab_composition", body)


def tab_loading_density():
    dens = [("1e8", r"$10^{8}$"), ("3e8", r"$3{\times}10^{8}$"), ("1e9", r"$10^{9}$")]
    lines = []
    for lab, form, pg0, pg24 in CELLULAR_LOADING:
        fshort = "I" if form == "SPION_I" else "II"
        cfg = lab.replace("_", r"\_")
        for tp, pg in ((r"0\,h", pg0), (r"24\,h", pg24)):
            cells = " & ".join(f"{cfe_from_loading(pg, CELL_DENSITY_LEVELS[dk]):.2f}"
                               for dk, _ in dens)
            lines.append(rf"{cfg} ({tp}) & {fshort} & {pg:.2f} & {cells} \\")
    header = r"""\begin{tabular}{@{}llccccc@{}}
\toprule
 & & & \multicolumn{3}{c}{$c_{\mathrm{Fe}}$ [mg/ml] at cell density} \\
\cmidrule(l){4-6}
Configuration & Formulation & pg\,Fe/cell & """ + " & ".join(d[1] for d in dens) + r""" \\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}"""
    _write_tex("tab_loading_density", header)


def _thr_fmt(v):
    return "--" if v is None else f"{v:.2f}"


def tab_thresholds():
    aopt_j = _load(f"{DET}/study_a_optimum.json")
    abas_j = _load(f"{DET}/study_a_baseline.json")
    aopt = aopt_j["thresholds_rose5"]
    abas = abas_j["thresholds_rose5"]
    opt_lab = _spec_label(aopt_j["spectrum"]["kvp"], aopt_j["spectrum"]["filters"]).replace(" ", "")
    bas_lab = _spec_label(abas_j["spectrum"]["kvp"], abas_j["spectrum"]["filters"]).replace(" ", "")
    dens = [("1e8", r"$10^{8}$"), ("3e8", r"$3{\times}10^{8}$"), ("1e9", r"$10^{9}$")]
    lines = []
    for det in ("EID", "PCD"):
        for dk, dlab in dens:
            row = [dlab]
            for spname, table in (("opt", aopt), ("base", abas)):
                for dose in ("low", "high"):
                    key = f"{det}_{dk}_{dose}_False"
                    row.append(_thr_fmt(table.get(key)))
            lines.append(f"{det} & " + " & ".join(row) + r" \\")
        if det == "EID":
            lines.append(r"\midrule")
    body = r"""\begin{tabular}{@{}llcccc@{}}
\toprule
 & & \multicolumn{2}{c}{Optimum (""" + opt_lab + r""")} & \multicolumn{2}{c}{Baseline (""" + bas_lab + r""")} \\
\cmidrule(lr){3-4}\cmidrule(l){5-6}
Detector & Cell density & low dose & high dose & low dose & high dose \\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
% Lowest detectable $c_{\mathrm{Fe}}$ [mg/ml] at Rose CNR $\geq 5$; bone off. -- = not reached in swept range."""
    _write_tex("tab_thresholds", body)


def tab_cnr():
    aopt_j = _load(f"{DET}/study_a_optimum.json")
    abas_j = _load(f"{DET}/study_a_baseline.json")
    aopt = _bone_off(aopt_j["rows"])
    abas = _bone_off(abas_j["rows"])
    opt_lab = _spec_label(aopt_j["spectrum"]["kvp"], aopt_j["spectrum"]["filters"]).replace(" ", "")
    bas_lab = _spec_label(abas_j["spectrum"]["kvp"], abas_j["spectrum"]["filters"]).replace(" ", "")

    def cnr_for(rows, cfg, det):
        vals = [r["cnr"] for r in rows if r["config"] == cfg and r["detector"] == det
                and r["dose"] == "high" and r["cell_density"] == "1e9"]
        return np.mean(vals) if vals else float("nan")

    configs = []
    for lab, form, pg0, pg24 in CELLULAR_LOADING:
        for tp in ("0h", "24h"):
            configs.append(f"{lab}_{tp}")

    lines = []
    for cfg in configs:
        eo = cnr_for(aopt, cfg, "EID"); po = cnr_for(aopt, cfg, "PCD")
        eb = cnr_for(abas, cfg, "EID"); pb = cnr_for(abas, cfg, "PCD")
        ratio = po / eo if eo else float("nan")
        lines.append(rf"{cfg.replace('_',' ')} & {eo:.1f} & {po:.1f} & {ratio:.2f} & "
                     rf"{eb:.1f} & {pb:.1f} \\")
    body = r"""\begin{tabular}{@{}lccccc@{}}
\toprule
 & \multicolumn{3}{c}{Optimum (""" + opt_lab + r""")} & \multicolumn{2}{c}{Baseline (""" + bas_lab + r""")} \\
\cmidrule(lr){2-4}\cmidrule(l){5-6}
Configuration & EID & PCD & PCD/EID & EID & PCD \\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
% Study A CNR at cell density $10^9$/cm$^3$, high dose, bone off."""
    _write_tex("tab_cnr", body)


# ============================================================================
def main():
    print("=== data-driven figures + tables (JSON only) ===")
    h1 = None
    h3 = fig_studyA_cnr()
    h5 = fig_studyB()
    h6 = fig_spectral()
    h7 = fig_eid_vs_pcd()
    tab_composition()
    tab_loading_density()
    tab_thresholds()
    tab_cnr()

    print("=== live-CONRAD figures (physics + recon) ===")
    try:
        h1 = fig_physics()
        fig_phantom_recon()
        fig_density_montage()
    except Exception as e:
        print(f"[warn] live-CONRAD figures skipped: {e!r}")

    print("\n=== HEADLINE NUMBERS ===")
    if h1: print(f"  fig_physics:    1 mg Fe/ml -> {h1:.1f} dHU @ 60 keV mono")
    print(f"  fig_studyA_cnr: peak PCD CNR (high dose) = {h3:.1f}")
    print(f"  fig_studyB:     peak PCD CNR (optimum, high dose) = {h5:.1f}")
    if h6: print(f"  fig_spectral:   optimum PCD ideal-observer CNR = {h6:.3f}")
    print(f"  fig_eid_vs_pcd: mean PCD/EID ratio (optimum) = {h7:.2f}x")


if __name__ == "__main__":
    main()
