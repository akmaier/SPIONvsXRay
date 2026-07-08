"""M2 — round geometric rabbit-scale phantom.

Analytic geometry (soft-tissue cylinder + cortical-bone rod + iron-loaded tumor
sphere) plus a voxelizer that returns three component volumes so polychromatic
projection (M4) can integrate each material's path length independently:

    mu(E; ray) = L_soft·mu_soft(E) + L_bone·mu_bone(E) + (∫c_Fe dl)·oxide(E)

Two tumor distributions (SPEC §5.1/§5.9):
  * 'homogeneous' — uniform iron through the sphere (Study A).
  * 'vessel'      — iron confined to 150 µm vessels at 10% volume, mass conserved
                    -> 10x local conc. Vessels are sub-resolution, so per-voxel
                    effective conc has mean = c_Fe with binomial partial-volume
                    variance (structural noise). (Study B)
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

import conrad_backend
from config import PHANTOM, tumor_iron_conc, VESSEL_DIAMETER_UM, VESSEL_VOLUME_FRACTION

_RNG = np.random.default_rng(20260708)   # fixed seed for reproducibility


@dataclass
class ComponentVolumes:
    """Per-voxel material composition (all arrays share the same grid)."""
    soft: np.ndarray      # soft-tissue fraction in [0,1]
    bone: np.ndarray      # cortical-bone fraction in [0,1]
    iron_mgml: np.ndarray # iron concentration [mg Fe/ml] (tumor voxels only)
    voxel_mm: float
    origin_mm: tuple      # (x0,y0,z0) of voxel-index (0,0,0) center


def _grid(dim, voxel_mm):
    nx, ny, nz = dim
    # centered coordinates [mm]
    x = (np.arange(nx) - (nx - 1) / 2.0) * voxel_mm
    y = (np.arange(ny) - (ny - 1) / 2.0) * voxel_mm
    z = (np.arange(nz) - (nz - 1) / 2.0) * voxel_mm
    return np.meshgrid(x, y, z, indexing="ij")


def build_components(c_form: float, model: str = "homogeneous",
                     dim=(256, 256, 192), voxel_mm: float = 0.78) -> ComponentVolumes:
    """Voxelize the phantom for one formulation concentration and tumor model."""
    X, Y, Z = _grid(dim, voxel_mm)
    ph = PHANTOM

    # body cylinder (axis along z)
    body_r = ph.body_diameter_mm / 2.0
    in_body = (X**2 + Y**2 <= body_r**2) & (np.abs(Z) <= ph.body_height_mm / 2.0)

    # cortical-bone rod (cylinder along z)
    bx, by, _ = ph.bone_center_mm
    in_bone = in_body & ((X - bx)**2 + (Y - by)**2 <= ph.bone_radius_mm**2)

    # tumor sphere
    tx, ty, tz = ph.tumor_center_mm
    tr = ph.tumor_radius_mm
    in_tumor = (X - tx)**2 + (Y - ty)**2 + (Z - tz)**2 <= tr**2

    soft = np.where(in_body & ~in_bone, 1.0, 0.0).astype(np.float32)
    bone = np.where(in_bone, 1.0, 0.0).astype(np.float32)

    c_fe = tumor_iron_conc(c_form)   # mean tumor iron [mg Fe/ml]
    iron = np.zeros(dim, dtype=np.float32)
    if c_fe > 0:
        if model == "homogeneous":
            iron[in_tumor] = c_fe
        elif model == "vessel":
            # per-voxel effective conc: 10% vessels at 10x conc, mass conserved.
            # sub-voxels per CT voxel ~ (voxel / vessel)^3, capped for stability.
            vessel_mm = VESSEL_DIAMETER_UM / 1000.0
            m = max(1, int(round((voxel_mm / vessel_mm) ** 3)))
            m = min(m, 4000)
            n_tumor = int(in_tumor.sum())
            vessel_frac = _RNG.binomial(m, VESSEL_VOLUME_FRACTION, size=n_tumor) / m
            local = vessel_frac * (c_fe / VESSEL_VOLUME_FRACTION)   # 10x conc where vessels
            iron[in_tumor] = local.astype(np.float32)
        else:
            raise ValueError(f"unknown tumor model {model!r}")

    origin = (-(dim[0] - 1) / 2.0 * voxel_mm,
              -(dim[1] - 1) / 2.0 * voxel_mm,
              -(dim[2] - 1) / 2.0 * voxel_mm)
    return ComponentVolumes(soft, bone, iron, voxel_mm, origin)


def summary_stats(cv: ComponentVolumes) -> dict:
    tumor = cv.iron_mgml[cv.iron_mgml > 0]
    return dict(
        soft_voxels=int(cv.soft.sum()),
        bone_voxels=int(cv.bone.sum()),
        tumor_voxels=int((cv.iron_mgml > 0).sum()),
        iron_mean=float(cv.iron_mgml[cv.iron_mgml >= 0].mean()) if cv.iron_mgml.size else 0.0,
        tumor_iron_mean=float(tumor.mean()) if tumor.size else 0.0,
        tumor_iron_std=float(tumor.std()) if tumor.size else 0.0,
        tumor_iron_max=float(tumor.max()) if tumor.size else 0.0,
    )


def visualize(outdir: str = None, c_form: float = 20.0):
    """Render axial slices for homogeneous vs vessel tumor -> dashboard figures."""
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if outdir is None:
        outdir = str(conrad_backend.REPO_ROOT / "results" / "phantom")
    os.makedirs(outdir, exist_ok=True)

    dim = (256, 256, 64)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.3))
    cv = build_components(c_form, "homogeneous", dim=dim)
    kz = dim[2] // 2

    # composite material map: soft=0.3, bone=1.0, tumor outlined by iron>0
    comp = cv.soft[:, :, kz] * 0.3 + cv.bone[:, :, kz] * 1.0
    comp[cv.iron_mgml[:, :, kz] > 0] = 0.6
    axes[0].imshow(comp.T, cmap="bone", origin="lower")
    axes[0].set_title("Phantom (soft + bone + tumor)"); axes[0].axis("off")

    # tumor iron map: homogeneous vs vessel (zoom on tumor)
    tx = int(round(PHANTOM.tumor_center_mm[0] / cv.voxel_mm + (dim[0] - 1) / 2))
    ty = int(round(PHANTOM.tumor_center_mm[1] / cv.voxel_mm + (dim[1] - 1) / 2))
    r = int(round(PHANTOM.tumor_radius_mm / cv.voxel_mm)) + 3
    for ax, model in zip(axes[1:], ["homogeneous", "vessel"]):
        cvi = build_components(c_form, model, dim=dim)
        sl = cvi.iron_mgml[tx - r:tx + r, ty - r:ty + r, kz]
        im = ax.imshow(sl.T, cmap="inferno", origin="lower", vmin=0)
        ax.set_title(f"tumor iron [mg Fe/ml] — {model}"); ax.axis("off")
        plt.colorbar(im, ax=ax, fraction=0.046)
    fig.suptitle(f"Round geometric phantom (c_form={c_form:g} mg/ml)")
    fig.tight_layout()
    fig.savefig(f"{outdir}/phantom_axial.png", dpi=130)
    plt.close(fig)
    print("[ok] wrote", outdir, "-> phantom_axial.png")


if __name__ == "__main__":
    for model in ["homogeneous", "vessel"]:
        cv = build_components(20.0, model, dim=(256, 256, 96))
        print(f"[{model}]", summary_stats(cv))
    visualize()
