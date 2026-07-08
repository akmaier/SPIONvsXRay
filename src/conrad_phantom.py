"""M2 (CONRAD-native) — ED-phantom as a real CONRAD AnalyticPhantom.

Builds the phantom at runtime via pyconrad from CONRAD classes: a soft-tissue
(water) body Cylinder + 7 SPION concentration insert Cylinders on a circle at
equal radius + a cortical-bone insert, each a PhysicalObject with a Material.
The SPION inserts use registered magnetite (Fe3O4-in-water) Mixture materials
(custom-materials pattern from edu.stanford.rsl.tutorial.physics.CreateCustomMaterial),
so CONRAD's own attenuation model is the source of truth.

Inserts are added AFTER the body so the PriorityRayTracer lets them override the
body material where they overlap (verified: water 135 mm + iron 25 mm = 160 mm).
"""
from __future__ import annotations
import math
import jpype

import conrad_backend
from config import C_FORM_LEVELS, tumor_iron_conc

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


def register_spion_materials():
    """Register a magnetite (Fe3O4-in-water) Mixture per concentration."""
    global _registered
    conrad_backend.setup()
    if _registered:
        return
    matpkg = CG("edu.stanford.rsl.conrad.physics.materials")
    utils = CG("edu.stanford.rsl.conrad.physics.materials.utils")
    DB = CG("edu.stanford.rsl.conrad.physics.materials.database").MaterialsDB
    WAC = utils.WeightedAtomicComposition
    MU = utils.MaterialUtils
    water_particles_1g = 1.0 / MU.computeMolarMass(WAC("H2O"))
    for c in C_FORM_LEVELS:
        c_fe = tumor_iron_conc(c)                      # mg Fe/ml
        grams_magnetite = (1e-3 * c_fe) / FE_FRACTION  # g magnetite per ml (~per g water)
        wac = WAC("H2O", water_particles_1g)
        if grams_magnetite > 0:
            wac.add("Fe3O4", grams_magnetite / M_FE3O4)
        m = matpkg.Mixture()
        m.setDensity(1.0 + grams_magnetite)
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


def build_phantom():
    """Return (scene, inserts) where inserts = list of dicts with layout + material."""
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

    add(_cyl(BODY_RADIUS_MM, 0.0, 0.0), DB.getMaterialWithName("water"))   # body first

    inserts = []
    n_slots = len(C_FORM_LEVELS) + 1
    for k, c in enumerate(C_FORM_LEVELS):
        theta = 2 * math.pi * k / n_slots + math.pi / 2
        cx, cy = INSERT_CIRCLE_MM * math.cos(theta), INSERT_CIRCLE_MM * math.sin(theta)
        add(_cyl(INSERT_RADIUS_MM, cx, cy), DB.getMaterialWithName(spion_name(c)))
        inserts.append(dict(name=spion_name(c), c_form=float(c), c_fe=tumor_iron_conc(c),
                            center_mm=(cx, cy), radius_mm=INSERT_RADIUS_MM))
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
