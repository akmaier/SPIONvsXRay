"""Single source of truth for all experiment parameters (see SPEC.md).

All lengths in mm unless noted. Concentrations in mg/ml. Energies in keV.
"""
from __future__ import annotations
from dataclasses import dataclass, field

# --------------------------------------------------------------------------
# Dose model (SPEC §5.2) — delivered mass, anchored at 6 mg for 10 mg/ml.
# --------------------------------------------------------------------------
FE_FRACTION = 0.724          # iron mass fraction of magnetite Fe3O4
DELIVERED_AT_10_MG = 6.0     # mg SPION delivered into the tumor at c_form = 10 mg/ml
TUMOR_VOLUME_CM3 = 8.0       # cm^3

# --------------------------------------------------------------------------
# Full-nanoparticle model (magnetite core + PAA coating).
# --------------------------------------------------------------------------
# The nanoparticle is a PAA-coated magnetite cluster (Heinen et al.): a magnetite
# (Fe3O4) core + a polyacrylic-acid (PAA, monomer (C3H4O2)n) coating. The reported
# whole-nanoparticle mass (mg SPION/ml) therefore = magnetite mass + coating mass.
#
# PAA_MASS_FRAC (phi) is the coating's fraction of the whole-particle mass. The
# article's supplementary TGA (mmc1.pdf, Table A.1 / Fig A.1 -- now in the repo root)
# fixes it PER FORMULATION: the iron-oxide residual is 87.7% (SPION I) / 66.4%
# (SPION II); on heating, magnetite -> hematite (+3.4%, article Eq A.1), so magnetite
# = R_m/1.034 and the PAA complement is phi ~= 0.15 (SPION I, 12 nm cores) and ~= 0.36
# (SPION II, 8 nm cores). TODO: this module still uses a SINGLE central phi = 0.15
# (exact for SPION I); per-formulation phi enters with the Study A per-configuration
# sampling (see HANDOFF.md). phi sets ONLY the reported whole-particle concentration
# (mg SPION/ml); it is NEGLIGIBLE for the X-ray mu (PAA is low-Z C/H/O, tissue/water-
# equivalent) and does NOT move the iron contrast. Independent variable: tumor c_Fe.
PAA_MASS_FRAC = 0.15          # phi (SPION I): default single-formulation value
PAA_MASS_FRAC_RANGE = (0.15, 0.36)   # SPION I .. SPION II (TGA, supplement Table A.1)

# Per-formulation coating (supplement TGA, Table A.1): SPION I (12 nm) phi = 0.15,
# SPION II (8 nm) phi = 0.36. build_phantom(formulation=...) picks one; the generic
# concentration-sweep phantom defaults to SPION I (the realistic high-loading particle).
PAA_MASS_FRAC_BY_FORMULATION = {"SPION_I": 0.15, "SPION_II": 0.36}
DEFAULT_FORMULATION = "SPION_I"


def coating_frac(formulation: str = DEFAULT_FORMULATION) -> float:
    """PAA coating mass fraction phi for a formulation (TGA, supplement Table A.1)."""
    return PAA_MASS_FRAC_BY_FORMULATION[formulation]


# --------------------------------------------------------------------------
# Study A concentrations: MEASURED cellular loading (Heinen et al. Fig 5, B16-F10
# melanoma). Each configuration = (formulation, hydrodynamic size) with its 0 h and
# 24 h iron load in pg Fe/cell. Converted to tumor mg Fe/ml via the tumor cell
# density, and simulated with that formulation's coating phi (SPION I 0.15 / II 0.36).
# --------------------------------------------------------------------------
CELL_DENSITY_PER_CM3 = 3.0e8   # default tumor cellularity (mid level below).
# Tumor cellularity is uncertain (solid tumors ~1e8-1e9 /cm^3), so it is a FACTOR:
# the Study A iron spans this range (c_Fe = pg Fe/cell x density x 1e-9). E.g. SPION I
# fresh 8.23 pg/cell -> 0.82 / 2.47 / 8.23 mg Fe/ml at 1e8 / 3e8 / 1e9.
CELL_DENSITY_LEVELS = {"1e8": 1.0e8, "3e8": 3.0e8, "1e9": 1.0e9}
CELLULAR_LOADING = [           # label, formulation, pg Fe/cell @ 0 h, @ 24 h
    ("I_113",  "SPION_I",  8.23, 3.78),
    ("II_115", "SPION_II", 3.86, 1.52),
    ("II_98",  "SPION_II", 3.60, 1.37),
    ("II_76",  "SPION_II", 2.51, 0.85),
]


def cfe_from_loading(pg_fe_per_cell: float, density: float = CELL_DENSITY_PER_CM3) -> float:
    """Tumor iron conc [mg Fe/ml] from cellular loading [pg Fe/cell] x density [cells/cm^3].

    c_Fe [mg/ml] = pg_fe * density * 1e-12 g/pg * 1e3 mg/g = pg_fe * density * 1e-9.
    """
    return pg_fe_per_cell * density * 1e-9


def study_a_inserts(density: float = CELL_DENSITY_PER_CM3):
    """Study A insert specs: each cellular-loading configuration x timepoint.

    Returns dicts {name, formulation, phi, pg_fe, timepoint, c_fe} -- the tumor iron
    grounded in the measured loading, each with its formulation's coating phi.
    """
    out = []
    for label, form, pg0, pg24 in CELLULAR_LOADING:
        for tp, pg in (("0h", pg0), ("24h", pg24)):
            out.append(dict(name=f"{label}_{tp}", formulation=form,
                            phi=coating_frac(form), pg_fe=pg, timepoint=tp,
                            c_fe=cfe_from_loading(pg, density)))
    return out

# Independent variable: article formulation concentration [mg SPION/ml]
C_FORM_LEVELS = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]

# Body / tumor-matrix material (SPEC §5.1/§5.2): iron is loaded into a SOFT-TISSUE
# matrix, NOT water. CONRAD ships no literal "soft_tissue"/ICRU-44 XML; its
# ICRU-family soft-tissue phantom proxy is the "body" Mixture (rho=1.0 g/cm^3,
# H/O, water-equivalent density -- matches SPEC "ICRU soft-tissue (water proxy)").
# This single name links the phantom body cylinder, the base-material sinogram key,
# and every mu / precorrection lookup, so a zero-iron insert == background exactly.
BODY_MATERIAL = "body"


def delivered_spion_mg(c_form: float) -> float:
    """Delivered SPION particle mass [mg] for a formulation concentration."""
    return DELIVERED_AT_10_MG * (c_form / 10.0)


def tumor_iron_conc(c_form: float) -> float:
    """Tumor iron concentration [mg Fe/ml] from formulation conc [mg SPION/ml].

    c_Fe = FE_FRACTION * (delivered mass / tumor volume) = 0.0543 * c_form.
    """
    c_spion = delivered_spion_mg(c_form) / TUMOR_VOLUME_CM3
    return FE_FRACTION * c_spion


def nanoparticle_conc_from_fe(c_fe: float, phi: float = PAA_MASS_FRAC) -> float:
    """Whole-nanoparticle concentration [mg SPION/ml] from tumor iron conc [mg Fe/ml].

    A whole particle = magnetite core + PAA coating. Iron is FE_FRACTION (0.724) of
    the magnetite, and magnetite is the (1 - phi) non-coating fraction of the whole
    particle, so the whole-particle iron fraction is FE_FRACTION*(1 - phi) and

        c_NP = c_Fe / (FE_FRACTION * (1 - phi)).

    This is the REPORTED mg SPION/ml (particle basis); phi is an estimate (see
    PAA_MASS_FRAC) and is negligible for the X-ray mu.
    """
    return c_fe / (FE_FRACTION * (1.0 - phi))


def tumor_paa_conc(c_form: float, phi: float = PAA_MASS_FRAC) -> float:
    """PAA-coating mass concentration in the tumor [mg PAA/ml] = phi * c_NP."""
    return phi * nanoparticle_conc_from_fe(tumor_iron_conc(c_form), phi)


# --------------------------------------------------------------------------
# Phantom (SPEC §5.1) — rabbit-scale soft tissue + bone + iron-loaded tumor.
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Phantom:
    body_diameter_mm: float = 110.0     # rabbit trunk (~11 cm) inside 20 cm FOV
    body_height_mm: float = 140.0
    body_material: str = "body"  # CONRAD ICRU-family soft-tissue proxy (see BODY_MATERIAL)

    # 8 cm^3 tumor sphere -> radius = (3V/4pi)^(1/3)
    tumor_volume_cm3: float = TUMOR_VOLUME_CM3
    tumor_center_mm: tuple = (25.0, 0.0, 0.0)   # offset from iso-center

    # cortical-bone insert (beam-hardening source), a rod
    bone_radius_mm: float = 8.0
    bone_center_mm: tuple = (-30.0, 0.0, 0.0)
    bone_material: str = "cortical_bone"

    @property
    def tumor_radius_mm(self) -> float:
        return (3.0 * self.tumor_volume_cm3 * 1000.0 / (4.0 * 3.141592653589793)) ** (1.0 / 3.0)


# --------------------------------------------------------------------------
# Acquisition geometry (SPEC §5.5) — standard C-arm cone beam.
# Defaults below are recorded from CONRAD's standard config at M4 and may be
# refined; FOV and projection count are fixed by the study.
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Geometry:
    num_projections: int = 500
    fov_mm: float = 200.0               # 20 cm reconstruction FOV
    angular_range_deg: float = 200.0    # short scan (180 + fan); confirm at M4
    source_iso_mm: float = 750.0        # SID  (CONRAD default; confirm)
    source_detector_mm: float = 1200.0  # SDD  (CONRAD default; confirm)
    detector_cols: int = 620
    detector_rows: int = 480
    detector_pixel_mm: float = 0.616    # confirm from CONRAD default


# --------------------------------------------------------------------------
# Reconstruction (SPEC §5.5)
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Recon:
    voxels: tuple = (512, 512, 384)
    voxel_mm: float = 200.0 / 512.0     # ~0.39 mm isotropic over 20 cm FOV
    bh_correction_states: tuple = (False, True)   # beam-hardening correction off/on


# --------------------------------------------------------------------------
# Spectrum & dose (SPEC §5.3) — CONRAD standard polychromatic spectrum.
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Spectrum:
    # CONRAD standard spectrum (confirmed at M3): 90 kVp, 10-150 keV, 0.5 keV,
    # mean 55.4 keV, W anode with characteristic lines (peak flux at 59.5 keV).
    kvp: float = 90.0
    e_min_kev: float = 10.0
    e_max_kev: float = 150.0
    e_delta_kev: float = 0.5
    photons_per_pixel: float = 70000.0   # documented reference dose (SPEC §5.3)


# --------------------------------------------------------------------------
# Dose as a factor (SPEC §5.3 fixed 70k -> now swept low/high per the C-arm dose
# research). photons/pixel = unattenuated I0 driving the Poisson/Gaussian noise.
# 70000 above stays the documented reference; the factorial sweeps these two.
# --------------------------------------------------------------------------
DOSE_LEVELS = {"low": 20000.0, "high": 100000.0}


# --------------------------------------------------------------------------
# Detectors (SPEC §5.4)
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Detectors:
    types: tuple = ("EID", "PCD")        # energy-integrating, photon-counting
    pcd_num_bins: int = 3
    # Iron has NO usable K-edge (7.1 keV), so thresholds are NOT placed at a
    # K-edge. They separate the photoelectric-rich low band (high iron contrast)
    # from the Compton high band. Optimized on the REAL CONRAD standard (90 kVp)
    # spectrum through the rabbit via src/spectral.py (M3):
    #   3-bin edges 10/37.5/50/90 keV -> 1.35x CNR vs EID (97% of ideal),
    #   2-bin split at 47.5 keV        -> 1.31x CNR vs EID.
    pcd_bin_edges_kev: tuple = (10.0, 37.5, 50.0, 90.0)
    pcd_bin_edges_2_kev: tuple = (10.0, 47.5, 90.0)


# --------------------------------------------------------------------------
# Spectral optimization sweep (filters + tube voltage) — see src/spectral.py.
# Iron contrast is photoelectric -> favors LOW energy. Hardening filters
# (Cu/Sn) HURT iron contrast; lower kVp and optimal PCD weighting help most.
# --------------------------------------------------------------------------
KVP_LEVELS = [60.0, 80.0, 100.0, 120.0]
FILTER_CONFIGS = {
    "Al2.5_baseline": [("aluminium", 2.5)],
    "Al0.5_soft":     [("aluminium", 0.5)],
    "Cu0.3_hard":     [("aluminium", 2.5), ("copper", 0.3)],
    "Sn0.5_hard":     [("aluminium", 2.5), ("tin", 0.5)],
    "Er0.1_quasimono":[("aluminium", 1.0), ("erbium", 0.1)],
}


# --------------------------------------------------------------------------
# Evaluation (SPEC §5.6/§5.7)
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Evaluation:
    noise_realizations: int = 30         # per cell; σ(CNR) estimate ~1/√(2·(n−1))
    rose_cnr_threshold: float = 5.0      # detectable if CNR >= 3..5 (report both)


# Phantom background and tumor distribution models (SPEC §5.1 / §5.9).
# Decision 2026-07-08: single round geometric phantom (no digital rabbit exists;
# ROBY=rat, MOBY=mouse). Realistic-anatomy arm dropped.
PHANTOM_BACKGROUNDS = ["round"]              # geometric rabbit-scale cylinder
# Two experiments (one comparable factor):
#   Study A homogeneous = cellular uptake (SPIONs internalised -> ~uniform tumor iron)
#   Study B vessel      = vascular/fresh delivery (SPIONs still in 150 µm vessels @10%
#                         volume, not yet taken up -> 10x local conc, heterogeneous)
TUMOR_MODELS = ["homogeneous", "vessel"]

# Study B vessel model (SPEC §5.9)
VESSEL_DIAMETER_UM = 150.0
VESSEL_VOLUME_FRACTION = 0.10        # vessels occupy 10% of tumor volume
# mass-conserved: local vessel iron concentration = tumor mean / volume fraction (10x)

# Reconstruction voxel size (isotropic). The fixed CL fan backprojector
# (conrad_ext) honors this via setSpacing; 512 px * 0.5 mm = 256 mm recon FOV.
RECON_VOXEL_MM = 0.5

PHANTOM = Phantom()
GEOMETRY = Geometry()
RECON = Recon()
SPECTRUM = Spectrum()
DETECTORS = Detectors()
EVAL = Evaluation()


def summary() -> str:
    lines = ["SPIONvsXRay configuration", "=" * 40]
    lines.append(f"Formulation levels (mg SPION/ml): {C_FORM_LEVELS}")
    lines.append(f"Full nanoparticle: magnetite core + PAA coating (phi={PAA_MASS_FRAC:g},"
                 f" range {PAA_MASS_FRAC_RANGE[0]:g}-{PAA_MASS_FRAC_RANGE[1]:g}; estimate)")
    lines.append("Tumor loading per level (iron mg Fe/ml -> whole-particle mg SPION/ml, PAA mg/ml):")
    for c in C_FORM_LEVELS:
        cfe = tumor_iron_conc(c)
        lines.append(f"  c_form={c:5.1f} -> delivered={delivered_spion_mg(c):5.2f} mg"
                     f" -> {cfe:.4f} mg Fe/ml -> {nanoparticle_conc_from_fe(cfe):.4f} mg SPION/ml"
                     f" (PAA {tumor_paa_conc(c):.4f} mg/ml)")
    lines.append(f"Tumor radius: {PHANTOM.tumor_radius_mm:.2f} mm (8 cm^3 sphere)")
    lines.append(f"Projections: {GEOMETRY.num_projections}, FOV: {GEOMETRY.fov_mm} mm")
    lines.append(f"Photons/pixel: {SPECTRUM.photons_per_pixel:.0f}")
    lines.append(f"Detectors: {DETECTORS.types}, PCD bins: {DETECTORS.pcd_num_bins}")
    n_cells = (len(C_FORM_LEVELS) * len(DETECTORS.types) * 2 * len(DOSE_LEVELS)
               * len(TUMOR_MODELS) * EVAL.noise_realizations)
    lines.append(f"Detectability cells: conc {len(C_FORM_LEVELS)} x detector "
                 f"{len(DETECTORS.types)} x bone 2 x dose {len(DOSE_LEVELS)} x "
                 f"distribution {len(TUMOR_MODELS)} x {EVAL.noise_realizations} reps "
                 f"= {n_cells}  (BH + short scan always on, not factors)")
    return "\n".join(lines)


if __name__ == "__main__":
    print(summary())
