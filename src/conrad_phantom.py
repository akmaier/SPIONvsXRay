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
from config import C_FORM_LEVELS, tumor_iron_conc, BODY_MATERIAL

CG = conrad_backend.class_getter

# geometry [mm]
BODY_RADIUS_MM = 80.0
INSERT_CIRCLE_MM = 50.0
INSERT_RADIUS_MM = 12.5
BONE_RADIUS_MM = 12.5
CYL_HEIGHT_MM = 200.0

FE_FRACTION = 0.724                       # magnetite Fe mass fraction
M_FE3O4 = 3 * 55.845 + 4 * 15.999         # g/mol

_registered = False


def spion_name(c_form: float) -> str:
    return f"SPION_c{c_form:g}"


def _soft_tissue_wac():
    """A FRESH WeightedAtomicComposition copy of the ICRU soft-tissue matrix.

    getWeightedAtomicComposition() returns the DB material's LIVE object, so
    mutating it (adding magnetite) would pollute the shared material and
    accumulate across concentrations. We copy its composition table into a new
    WAC per call; the SPION_c0 build below uses the identical copy, so background
    (SPION_c0) == zero-iron insert exactly.
    """
    utils = CG("edu.stanford.rsl.conrad.physics.materials.utils")
    WAC = utils.WeightedAtomicComposition
    DB = CG("edu.stanford.rsl.conrad.physics.materials.database").MaterialsDB
    src = DB.getMaterialWithName(BODY_MATERIAL).getWeightedAtomicComposition().getCompositionTable()
    w = WAC()
    for k in src.keySet():
        w.add(str(k), float(src.get(k)))
    return w


def register_spion_materials():
    """Register a magnetite-oxide-in-soft-tissue Mixture per concentration.

    Iron oxide (Fe3O4) is ADDED to the ICRU soft-tissue matrix (config.BODY_MATERIAL),
    NOT diluted in water (SPEC §5.2). The magnetite mole count is referenced to the
    soft-tissue matrix's own molar basis so the delivered iron mass is exactly
    c_Fe [mg Fe/ml] (c_Fe = 0.0543*c_form, FE_FRACTION magnetite). SPION_c0 carries
    zero iron, so it equals plain soft tissue and ΔHU(c_Fe=0) == 0.
    """
    global _registered
    conrad_backend.setup()
    if _registered:
        return
    matpkg = CG("edu.stanford.rsl.conrad.physics.materials")
    utils = CG("edu.stanford.rsl.conrad.physics.materials.utils")
    DB = CG("edu.stanford.rsl.conrad.physics.materials.database").MaterialsDB
    MU = utils.MaterialUtils
    # grams of soft-tissue matrix per one WAC formula-unit (its own molar mass)
    mm_soft = float(MU.computeMolarMass(DB.getMaterialWithName(BODY_MATERIAL)
                                        .getWeightedAtomicComposition()))
    for c in C_FORM_LEVELS:
        c_fe = tumor_iron_conc(c)                      # mg Fe/ml
        grams_magnetite = (1e-3 * c_fe) / FE_FRACTION  # g magnetite per g matrix
        wac = _soft_tissue_wac()
        if grams_magnetite > 0:
            # magnetite moles per one soft-tissue formula-unit (mole-consistent add)
            wac.add("Fe3O4", grams_magnetite * mm_soft / M_FE3O4)
        m = matpkg.Mixture()
        m.setDensity(1.0 + grams_magnetite)            # rho raised by ~0.001*c_Fe (negligible)
        m.setName(spion_name(c))
        m.setWeightedAtomicComposition(wac)
        DB.put(m)
    _registered = True


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


def build_phantom(with_bone=True):
    """Return (scene, inserts) where inserts = list of dicts with layout + material.

    with_bone=False omits the cortical-bone rod (used only for clean display
    galleries; the quantitative study always keeps the bone beam-hardening source).
    """
    conrad_backend.setup()
    register_spion_materials()
    phys = CG("edu.stanford.rsl.conrad.physics")
    DB = CG("edu.stanford.rsl.conrad.physics.materials.database").MaterialsDB
    scene = _new_scene()

    def add(shape, material):
        po = phys.PhysicalObject()
        po.setMaterial(material)
        po.setShape(shape)
        scene.add(po)

    # body = zero-iron soft tissue (SPION_c0), so background == zero-iron insert
    add(_cyl(BODY_RADIUS_MM, 0.0, 0.0), DB.getMaterialWithName(spion_name(0.0)))  # body first

    inserts = []
    n_slots = len(C_FORM_LEVELS) + 1
    for k, c in enumerate(C_FORM_LEVELS):
        theta = 2 * math.pi * k / n_slots + math.pi / 2
        cx, cy = INSERT_CIRCLE_MM * math.cos(theta), INSERT_CIRCLE_MM * math.sin(theta)
        add(_cyl(INSERT_RADIUS_MM, cx, cy), DB.getMaterialWithName(spion_name(c)))
        inserts.append(dict(name=spion_name(c), c_form=float(c), c_fe=tumor_iron_conc(c),
                            center_mm=(cx, cy), radius_mm=INSERT_RADIUS_MM))
    if with_bone:
        theta = 2 * math.pi * len(C_FORM_LEVELS) / n_slots + math.pi / 2
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
