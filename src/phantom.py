"""M2 (rev.) — ED-phantom-style module: concentration inserts on a circle.

Like a CT calibration (electron-density) phantom: a soft-tissue cylinder with the
SPION concentration inserts arranged on a circle at EQUAL radius, plus a
cortical-bone insert. All concentrations are imaged in one scan and sit at the
same radial position, so beam-hardening cupping affects them equally and they are
directly comparable. Each insert is a 2.5 cm disk (~8 cm^3-equivalent tumor).

Returns 2D component maps (soft fraction, bone fraction, iron mg/ml) for fan-beam
projection, and the insert table for ROI measurement.

Tumor models (per insert):
  'homogeneous' — uniform iron; 'vessel' — 150 um vessels @10% (mass conserved,
  10x local conc, sub-resolution -> per-voxel binomial partial-volume variance).
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np

import conrad_backend
from config import (C_FORM_LEVELS, tumor_iron_conc,
                    VESSEL_DIAMETER_UM, VESSEL_VOLUME_FRACTION)

_RNG = np.random.default_rng(20260708)

BODY_RADIUS_MM = 80.0        # fits inside the 20 cm FOV
INSERT_CIRCLE_MM = 50.0      # radius of the ring of inserts
INSERT_RADIUS_MM = 12.5      # 2.5 cm dia ~ 8 cm^3-equivalent
BONE_RADIUS_MM = 12.5


@dataclass
class EDPhantom:
    soft: np.ndarray
    bone: np.ndarray
    iron: np.ndarray            # mg Fe/ml per pixel
    voxel_cm: float
    inserts: list = field(default_factory=list)   # dicts: name, center_mm, radius_mm, c_form, c_fe


def _disk(X, Y, cx, cy, r):
    return (X - cx) ** 2 + (Y - cy) ** 2 <= r ** 2


def build_ed_phantom(model="homogeneous", N=512, fov_mm=200.0) -> EDPhantom:
    voxel_mm = fov_mm / N
    ax = (np.arange(N) - (N - 1) / 2.0) * voxel_mm
    X, Y = np.meshgrid(ax, ax, indexing="ij")   # X=axis0(rows), Y=axis1(cols)

    soft = _disk(X, Y, 0, 0, BODY_RADIUS_MM).astype(np.float32)
    bone = np.zeros((N, N), np.float32)
    iron = np.zeros((N, N), np.float32)

    inserts = []
    n_conc = len(C_FORM_LEVELS)
    n_slots = n_conc + 1                       # + bone
    for k, c_form in enumerate(C_FORM_LEVELS):
        theta = 2 * np.pi * k / n_slots + np.pi / 2   # start at top
        cx = INSERT_CIRCLE_MM * np.cos(theta)
        cy = INSERT_CIRCLE_MM * np.sin(theta)
        mask = _disk(X, Y, cx, cy, INSERT_RADIUS_MM)
        c_fe = tumor_iron_conc(c_form)
        _fill_iron(iron, mask, c_fe, model, voxel_mm)
        inserts.append(dict(name=f"c{c_form:g}", center_mm=(float(cx), float(cy)),
                            radius_mm=INSERT_RADIUS_MM, c_form=float(c_form), c_fe=float(c_fe)))

    # bone insert in the last slot
    theta = 2 * np.pi * n_conc / n_slots + np.pi / 2
    bx, by = INSERT_CIRCLE_MM * np.cos(theta), INSERT_CIRCLE_MM * np.sin(theta)
    bmask = _disk(X, Y, bx, by, BONE_RADIUS_MM)
    bone[bmask] = 1.0
    soft[bmask] = 0.0
    inserts.append(dict(name="bone", center_mm=(float(bx), float(by)),
                        radius_mm=BONE_RADIUS_MM, c_form=None, c_fe=0.0))

    return EDPhantom(soft, bone, iron, voxel_mm / 10.0, inserts)


def _fill_iron(iron, mask, c_fe, model, voxel_mm):
    if c_fe <= 0:
        return
    if model == "homogeneous":
        iron[mask] = c_fe
    elif model == "vessel":
        vessel_mm = VESSEL_DIAMETER_UM / 1000.0
        m = int(np.clip(round((voxel_mm / vessel_mm) ** 3), 1, 4000))
        nvox = int(mask.sum())
        frac = _RNG.binomial(m, VESSEL_VOLUME_FRACTION, size=nvox) / m
        iron[mask] = (frac * (c_fe / VESSEL_VOLUME_FRACTION)).astype(np.float32)
    else:
        raise ValueError(f"unknown tumor model {model!r}")


def visualize(outdir=None):
    import os, matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if outdir is None:
        outdir = str(conrad_backend.REPO_ROOT / "results" / "phantom")
    os.makedirs(outdir, exist_ok=True)
    ph = build_ed_phantom("homogeneous", N=512)
    comp = ph.soft * 0.3 + ph.bone * 1.0
    comp[ph.iron > 0] = 0.6
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    ax[0].imshow(comp.T, cmap="bone", origin="lower")
    ax[0].set_title("ED phantom: soft + bone + concentration inserts")
    im = ax[1].imshow(ph.iron.T, cmap="inferno", origin="lower")
    ax[1].set_title("iron [mg Fe/ml] per insert")
    plt.colorbar(im, ax=ax[1], fraction=0.046)
    for a in ax:
        a.axis("off")
    # annotate concentrations
    for ins in ph.inserts:
        cx, cy = ins["center_mm"]
        px = cx / (ph.voxel_cm * 10) + 256
        py = cy / (ph.voxel_cm * 10) + 256
        lbl = "bone" if ins["name"] == "bone" else f"{ins['c_form']:g}"
        ax[0].text(px, py, lbl, color="cyan", ha="center", va="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{outdir}/phantom_axial.png", dpi=130)
    plt.close(fig)
    print("[ok] wrote", outdir, "-> phantom_axial.png; inserts:",
          [(i["name"], round(i["c_fe"], 3)) for i in ph.inserts])


if __name__ == "__main__":
    visualize()
