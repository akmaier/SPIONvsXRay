"""M2 (CONRAD-native) — ED-phantom as a real CONRAD AnalyticPhantom.

Builds the phantom at runtime via pyconrad from CONRAD classes: an ICRU
SOFT-TISSUE (config.BODY_MATERIAL) body Cylinder + 7 SPION concentration insert
Cylinders on a circle at equal radius + a cortical-bone insert, each a
PhysicalObject with a Material. The SPION inserts use registered magnetite-oxide
Mixture materials (custom-materials pattern from
edu.stanford.rsl.tutorial.physics.CreateCustomMaterial) built by ADDING magnetite
(Fe3O4) to the SAME soft-tissue matrix (NOT water; SPEC §5.2), so CONRAD's own
attenuation model is the source of truth and the zero-iron insert (SPION_c0)
equals plain soft tissue -- ΔHU is pure iron.

The body cylinder is assigned the SPION_c0 material (soft-tissue matrix, zero
iron), so background == zero-iron insert exactly and ΔHU(c_Fe=0) == 0.

Inserts are added AFTER the body so the PriorityRayTracer lets them override the
body material where they overlap.
"""
from __future__ import annotations
import math
import jpype

import conrad_backend
from config import (C_FORM_LEVELS, tumor_iron_conc, tumor_paa_conc,
                    nanoparticle_conc_from_fe, PAA_MASS_FRAC, BODY_MATERIAL,
                    coating_frac, DEFAULT_FORMULATION,
                    study_a_inserts, study_b_inserts, CELL_DENSITY_PER_CM3,
                    VESSEL_VOLUME_FRACTION)

CG = conrad_backend.class_getter

# geometry [mm]
BODY_RADIUS_MM = 80.0
INSERT_CIRCLE_MM = 50.0
INSERT_RADIUS_MM = 12.5
BONE_RADIUS_MM = 12.5
CYL_HEIGHT_MM = 200.0

FE_FRACTION = 0.724                       # magnetite Fe mass fraction
M_FE3O4 = 3 * 55.845 + 4 * 15.999         # g/mol

_registered_formulation = None


def spion_name(c_form: float) -> str:
    return f"SPION_c{c_form:g}"


# magnetite element mass fractions (Fe3O4 = 3 Fe + 4 O)
M_FE = 55.845
M_O = 15.999
M_C = 12.011
M_H = 1.008
FE_MASS_FRAC_FE3O4 = 3 * M_FE / M_FE3O4
O_MASS_FRAC_FE3O4 = 4 * M_O / M_FE3O4

# PAA coating: polyacrylic acid, monomer (C3H4O2)n. Low-Z (C/H/O), tissue/water-
# equivalent -> negligible for the X-ray mu, but included so the FULL nanoparticle
# (magnetite core + PAA coating) is simulated, not just the iron.
M_PAA_MONOMER = 3 * M_C + 4 * M_H + 2 * M_O          # (C3H4O2) g/mol = 72.06
C_MASS_FRAC_PAA = 3 * M_C / M_PAA_MONOMER
H_MASS_FRAC_PAA = 4 * M_H / M_PAA_MONOMER
O_MASS_FRAC_PAA = 2 * M_O / M_PAA_MONOMER


def _spion_wac(grams_magnetite, grams_paa=0.0):
    """Build a WAC for 1 g of ICRU soft-tissue matrix + grams_magnetite g of Fe3O4
    + grams_paa g of PAA coating (polyacrylic acid, monomer C3H4O2).

    CRITICAL: WeightedAtomicComposition.getCompositionTable() stores per-element
    *mass* (moles*atomicMass), and WAC.add(element, moles) multiplies its argument
    by the atomic mass. The 4f11336 regression copied the soft-tissue table back in
    via add(elem, src.get(elem)) -- feeding already-mass-weighted values through the
    *mass again, inflating the matrix ~16x (O jumped 14.2 -> 227) so the small Fe3O4
    add washed out (Fe mass fraction ~22x too small; mu even dropped below body).

    Instead we work purely in the mass basis: scale the soft-tissue composition to
    1 gram total, then add magnetite as elemental *masses* (Fe + O) via a TreeMap
    passed to setCompositionTable(). This mirrors the known-good pre-4f11336 water
    host (1 g H2O + grams_magnetite g Fe3O4) but with the soft-tissue matrix as host.
    grams_magnetite = 0 reproduces plain soft tissue exactly (ΔHU == 0).
    """
    utils = CG("edu.stanford.rsl.conrad.physics.materials.utils")
    WAC = utils.WeightedAtomicComposition
    MU = utils.MaterialUtils
    DB = CG("edu.stanford.rsl.conrad.physics.materials.database").MaterialsDB
    TreeMap = jpype.JClass("java.util.TreeMap")

    body_wac = DB.getMaterialWithName(BODY_MATERIAL).getWeightedAtomicComposition()
    src = body_wac.getCompositionTable()               # element -> mass (moles*A)
    mm_soft = float(MU.computeMolarMass(body_wac))      # total mass of one formula unit

    tbl = TreeMap()
    for k in src.keySet():
        tbl.put(str(k), float(src.get(k)) / mm_soft)   # scale to 1 g soft tissue
    def _add(elem, mass):
        cur = float(tbl.get(elem)) if tbl.containsKey(elem) else 0.0
        tbl.put(elem, cur + mass)
    if grams_magnetite > 0:
        _add("Fe", grams_magnetite * FE_MASS_FRAC_FE3O4)
        _add("O", grams_magnetite * O_MASS_FRAC_FE3O4)
    if grams_paa > 0:                                  # PAA coating: C/H/O (low-Z)
        _add("C", grams_paa * C_MASS_FRAC_PAA)
        _add("H", grams_paa * H_MASS_FRAC_PAA)
        _add("O", grams_paa * O_MASS_FRAC_PAA)
    w = WAC()
    w.setCompositionTable(tbl)
    return w


def register_spion_materials(formulation=None):
    """Register a FULL-nanoparticle-in-soft-tissue Mixture per concentration.

    The whole nanoparticle (magnetite Fe3O4 core + PAA coating, monomer C3H4O2) is
    ADDED to the ICRU soft-tissue matrix (config.BODY_MATERIAL), NOT diluted in water
    (SPEC §5.2). The delivered iron mass is exactly c_Fe [mg Fe/ml]
    (c_Fe = 0.0543*c_form, FE_FRACTION magnetite); the magnetite mass = c_Fe/0.724 and
    the PAA-coating mass = phi*c_NP (config.tumor_paa_conc). phi is the formulation's
    PAA coating fraction (supplement TGA, Table A.1): SPION I = 0.15, SPION II = 0.36.
    The PAA coating is low-Z (C/H/O), tissue/water-equivalent -> small for mu but
    included so the full particle is simulated. SPION_c0 carries zero particle mass, so
    it equals plain soft tissue and ΔHU(c_Fe=0) == 0. Re-registers when the formulation
    changes (materials keep the same SPION_c{c} names; one formulation at a time).
    """
    global _registered_formulation
    conrad_backend.setup()
    if formulation is None:
        formulation = DEFAULT_FORMULATION
    if _registered_formulation == formulation:
        return
    phi = coating_frac(formulation)                    # per-formulation PAA coating
    matpkg = CG("edu.stanford.rsl.conrad.physics.materials")
    DB = CG("edu.stanford.rsl.conrad.physics.materials.database").MaterialsDB
    for c in C_FORM_LEVELS:
        c_fe = tumor_iron_conc(c)                      # mg Fe/ml
        grams_magnetite = (1e-3 * c_fe) / FE_FRACTION  # g magnetite per g matrix
        grams_paa = 1e-3 * tumor_paa_conc(c, phi)      # g PAA coating per g matrix (phi)
        wac = _spion_wac(grams_magnetite, grams_paa)
        m = matpkg.Mixture()
        m.setDensity(1.0 + grams_magnetite + grams_paa)  # rho raised by particle mass
        m.setName(spion_name(c))
        m.setWeightedAtomicComposition(wac)
        DB.put(m)
    _registered_formulation = formulation


def register_loading_materials(specs):
    """Register a full-particle material per Study A cellular-loading spec.

    Unlike register_spion_materials (one formulation, concentration sweep), each spec
    carries its OWN c_fe and coating phi, so SPION I (phi=0.15) and SPION II (phi=0.36)
    inserts coexist in one phantom. Always (re)registers the zero-iron SPION_c0 body
    material too. Materials are keyed by the spec name (e.g. I_113_0h).
    """
    global _registered_formulation
    conrad_backend.setup()
    matpkg = CG("edu.stanford.rsl.conrad.physics.materials")
    DB = CG("edu.stanford.rsl.conrad.physics.materials.database").MaterialsDB
    allspecs = [dict(name=spion_name(0.0), c_fe=0.0, phi=PAA_MASS_FRAC)] + list(specs)
    for s in allspecs:
        c_fe, phi = s["c_fe"], s["phi"]
        grams_magnetite = (1e-3 * c_fe) / FE_FRACTION           # g magnetite / g matrix
        c_np = nanoparticle_conc_from_fe(c_fe, phi) if c_fe > 0 else 0.0
        grams_paa = 1e-3 * phi * c_np                           # PAA coating mass = phi*c_NP
        wac = _spion_wac(grams_magnetite, grams_paa)
        m = matpkg.Mixture()
        m.setDensity(1.0 + grams_magnetite + grams_paa)
        m.setName(s["name"])
        m.setWeightedAtomicComposition(wac)
        DB.put(m)
    _registered_formulation = None    # invalidate the sweep cache (materials replaced)


def register_vessel_materials(specs):
    """Register a full-particle material per Study B vessel insert (HOMOGENIZED).

    Same material machinery as register_loading_materials (each insert its OWN mean
    iron c_fe + coating phi, water body + soft-tissue matrix + magnetite + PAA), but
    keyed on the Study B insert names (VES_<level>). The iron is the tumor-MEAN
    (= vessel_level * VESSEL_VOLUME_FRACTION), since the 0.5 mm voxel cannot resolve
    the 150 um vessels; the vessel heterogeneity re-enters as a second-order model
    downstream, not in the material. Always (re)registers the zero-iron SPION_c0 body
    material too.
    """
    global _registered_formulation
    conrad_backend.setup()
    matpkg = CG("edu.stanford.rsl.conrad.physics.materials")
    DB = CG("edu.stanford.rsl.conrad.physics.materials.database").MaterialsDB
    allspecs = [dict(name=spion_name(0.0), c_fe=0.0, phi=PAA_MASS_FRAC)] + list(specs)
    for s in allspecs:
        c_fe, phi = s["c_fe"], s["phi"]
        grams_magnetite = (1e-3 * c_fe) / FE_FRACTION           # g magnetite / g matrix
        c_np = nanoparticle_conc_from_fe(c_fe, phi) if c_fe > 0 else 0.0
        grams_paa = 1e-3 * phi * c_np                           # PAA coating mass = phi*c_NP
        wac = _spion_wac(grams_magnetite, grams_paa)
        m = matpkg.Mixture()
        m.setDensity(1.0 + grams_magnetite + grams_paa)
        m.setName(s["name"])
        m.setWeightedAtomicComposition(wac)
        DB.put(m)
    _registered_formulation = None    # invalidate the sweep cache (materials replaced)


def _new_scene():
    try:
        return CG("edu.stanford.rsl.conrad.rendering").PrioritizableScene()
    except Exception:
        s = CG("edu.stanford.rsl.conrad.phantom").SpherePhantom()
        s.clear()
        return s


def _cyl(radius, cx, cy):
    shapes = CG("edu.stanford.rsl.conrad.geometry.shapes.simple")
    cyl = shapes.Cylinder(radius, radius, CYL_HEIGHT_MM)
    if cx != 0.0 or cy != 0.0:
        Tr = CG("edu.stanford.rsl.conrad.geometry.transforms").Translation
        pt = shapes.PointND(jpype.JArray(jpype.JDouble)([cx, cy, 0.0]))
        cyl.applyTransform(Tr(pt.getAbstractVector()))
    return cyl


def build_phantom(with_bone=True, formulation=None, loading=False, density=None,
                  study_b=False):
    """Return (scene, inserts) where inserts = list of dicts with layout + material.

    with_bone=False omits the cortical-bone rod (used only for clean display
    galleries; the quantitative study always keeps the bone beam-hardening source).

    Three concentration models (loading and study_b are mutually exclusive):
      loading=False, study_b=False (default): the concentration SWEEP (C_FORM_LEVELS)
        of one `formulation` (phi from config: SPION_I=0.15, SPION_II=0.36).
      loading=True: STUDY A -- the measured cellular loading (config.study_a_inserts):
        one insert per (formulation, timepoint) configuration, iron = pg Fe/cell x
        `density` (config.CELL_DENSITY_PER_CM3), each with its formulation's phi. The
        zero-iron SPION_c0 is the body + a reference insert.
      study_b=True: STUDY B -- fresh vascular injection (config.study_b_inserts): one
        insert per vessel-LOCAL level, HOMOGENIZED to its tumor-mean iron
        (c_fe = vessel_level * VESSEL_VOLUME_FRACTION, ~0.5..1.5 mg Fe/ml). Fresh
        SPION I injection (phi=0.15). Each insert is tagged tumor_model="vessel",
        vessel_level, c_fe. Includes the zero-iron SPION_c0 insert; bone stays a factor.
    """
    if loading and study_b:
        raise ValueError("loading and study_b are mutually exclusive")
    conrad_backend.setup()
    if study_b:
        specs = study_b_inserts()
        register_vessel_materials(specs)
        insert_specs = [dict(name=spion_name(0.0), c_form=0.0, c_fe=0.0,
                             tumor_model="vessel", vessel_level=0.0, phi=PAA_MASS_FRAC)]
        insert_specs += [dict(name=s["name"], c_form=None, c_fe=s["c_fe"],
                              tumor_model="vessel", vessel_level=s["vessel_level"],
                              phi=s["phi"]) for s in specs]
    elif loading:
        if density is None:
            density = CELL_DENSITY_PER_CM3
        specs = study_a_inserts(density)
        register_loading_materials(specs)
        insert_specs = [dict(name=spion_name(0.0), c_form=0.0, c_fe=0.0,
                             formulation=DEFAULT_FORMULATION, timepoint="-", pg_fe=0.0)]
        insert_specs += [dict(name=s["name"], c_form=None, c_fe=s["c_fe"],
                              formulation=s["formulation"], timepoint=s["timepoint"],
                              pg_fe=s["pg_fe"]) for s in specs]
    else:
        register_spion_materials(formulation)
        insert_specs = [dict(name=spion_name(c), c_form=float(c), c_fe=tumor_iron_conc(c))
                        for c in C_FORM_LEVELS]

    phys = CG("edu.stanford.rsl.conrad.physics")
    DB = CG("edu.stanford.rsl.conrad.physics.materials.database").MaterialsDB
    scene = _new_scene()

    def add(shape, material):
        po = phys.PhysicalObject()
        po.setMaterial(material)
        po.setShape(shape)
        scene.add(po)

    # body/background = WATER; the inserts are iron-loaded SOFT TISSUE embedded in it
    # (the water-calibrated BH precorrection then flattens the body exactly). The
    # zero-iron SPION_c0 insert (soft tissue) carries the soft-tissue-vs-water offset,
    # which measure_inserts' c0 subtraction removes -> ΔHU stays pure iron.
    add(_cyl(BODY_RADIUS_MM, 0.0, 0.0), DB.getMaterialWithName("water"))  # water background

    inserts = []
    n_slots = len(insert_specs) + (1 if with_bone else 0)
    for k, sp in enumerate(insert_specs):
        theta = 2 * math.pi * k / n_slots + math.pi / 2
        cx, cy = INSERT_CIRCLE_MM * math.cos(theta), INSERT_CIRCLE_MM * math.sin(theta)
        add(_cyl(INSERT_RADIUS_MM, cx, cy), DB.getMaterialWithName(sp["name"]))
        inserts.append({**sp, "center_mm": (cx, cy), "radius_mm": INSERT_RADIUS_MM})
    if with_bone:
        theta = 2 * math.pi * len(insert_specs) / n_slots + math.pi / 2
        bx, by = INSERT_CIRCLE_MM * math.cos(theta), INSERT_CIRCLE_MM * math.sin(theta)
        add(_cyl(BONE_RADIUS_MM, bx, by), DB.getMaterialWithName("bone"))
        inserts.append(dict(name="bone", c_form=None, c_fe=0.0, center_mm=(bx, by), radius_mm=BONE_RADIUS_MM))
    return scene, inserts


if __name__ == "__main__":
    scene, inserts = build_phantom()
    DB = CG("edu.stanford.rsl.conrad.physics.materials.database").MaterialsDB
    print("scene objects:", sum(1 for _ in scene))
    print("SPION material attenuation @60 keV (should rise with c):")
    for c in C_FORM_LEVELS:
        mat = DB.getMaterialWithName(spion_name(c))
        AT = CG("edu.stanford.rsl.conrad.physics.materials.utils").AttenuationType
        mu = float(mat.getAttenuation(60.0, AT.TOTAL_WITH_COHERENT_ATTENUATION))
        print(f"  c_form={c:5.1f} (c_Fe={tumor_iron_conc(c):.3f}): density={mat.getDensity():.5f}  mu60={mu:.5f} /cm")
