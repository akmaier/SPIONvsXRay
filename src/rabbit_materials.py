"""Re-materialize the single-energy RabbitCT reconstruction into a polychromatic
phantom: HU calibration -> (material, density) segmentation, so each voxel gets an
energy-dependent attenuation for poly-energy forward projection.

Calibration anchors from the reference_256.vol histogram: air = raw 0 -> -1000 HU,
soft-tissue peak = raw 1023 -> 0 HU. HU->material follows the standard
Schneider-style piecewise map (air / soft tissue with density / soft-bone blend).

Returns per-voxel density and bone-fraction maps; the poly attenuation is then
  mu(E) = rho * [(1-f) (mu/rho)_soft(E) + f (mu/rho)_bone(E)]
which feeds a per-material cone-beam projection (src/rabbit3d.py) and the CONRAD
GPU spectral detectors.
"""
from __future__ import annotations
import numpy as np

import rabbitct
import conrad_backend

RAW_AIR = 0.0
RAW_WATER = 1023.0


def raw_to_hu(vol):
    return (vol - RAW_WATER) * 1000.0 / (RAW_WATER - RAW_AIR)


# HU class boundaries (from the reference_256 histogram: bimodal soft tissue with
# an adipose peak ~ -210 HU and a muscle peak ~ +10 HU, valley ~ -80 HU).
HU_AIR = -350.0
HU_FAT_MUSCLE = -80.0
HU_BONE = 120.0
HU_CORTICAL = 1500.0

# CONRAD material for each class (verified present, energy-dependent coefficients).
CLASS_MATERIAL = {0: None, 1: "adipose", 2: "muscle", 3: "bone"}
CLASS_DENSITY = {1: 0.92, 2: 1.05, 3: 1.92}   # CONRAD densities


def segment(vol):
    """HU -> (density, bone_fraction, label) with adipose/muscle separated.

    label: 0 air, 1 adipose (fat), 2 muscle/soft tissue, 3 muscle/bone blend.
    The bone class blends muscle<->cortical bone by `bone_frac` for partial volume.
    """
    hu = raw_to_hu(vol)
    dens = np.zeros_like(hu, np.float32)
    f = np.zeros_like(hu, np.float32)
    label = np.zeros(hu.shape, np.uint8)

    fat = (hu >= HU_AIR) & (hu < HU_FAT_MUSCLE)
    muscle = (hu >= HU_FAT_MUSCLE) & (hu < HU_BONE)
    bone = hu >= HU_BONE

    label[fat] = 1
    dens[fat] = CLASS_DENSITY[1]
    label[muscle] = 2
    dens[muscle] = CLASS_DENSITY[2]
    fb = np.clip((hu[bone] - HU_BONE) / (HU_CORTICAL - HU_BONE), 0.0, 1.0)
    f[bone] = fb
    label[bone] = 3
    dens[bone] = np.clip(1.05 + fb * (1.92 - 1.05), 1.05, 1.92)
    return dict(hu=hu, density=dens, bone_frac=f, label=label)


def material_massatten(energies_kev):
    """CONRAD mass attenuation (mu/rho) [cm^2/g] for the segmentation materials --
    the energy-dependent (spectral) coefficients used for poly-energy projection."""
    import materials
    return {name: materials.mass_attenuation(name, energies_kev)
            for name in ("adipose", "muscle", "bone")}


def main():
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    root = conrad_backend.REPO_ROOT
    vol = rabbitct.load_reference_volume(str(root / "data" / "rabbitct" / "reference_256.vol"), 256)
    seg = segment(vol)
    label, dens = seg["label"], seg["density"]

    frac = {n: float((label == i).mean()) for i, n in enumerate(["air", "adipose", "muscle", "bone"])}
    print("volume fractions:", {k: round(v, 4) for k, v in frac.items()})
    body = label > 0
    print(f"soft-tissue split: adipose {100*(label==1).sum()/body.sum():.0f}% / "
          f"muscle {100*(label==2).sum()/body.sum():.0f}% / bone {100*(label==3).sum()/body.sum():.0f}% of body")
    print(f"density range (non-air): {dens[body].min():.2f}-{dens[body].max():.2f} g/cm^3")
    ma = material_massatten(np.array([30.0, 60.0, 90.0]))
    print("CONRAD (mu/rho) [cm^2/g] @30/60/90 keV used per class:")
    for k, v in ma.items():
        print(f"  {k:8s} {v[0]:.4f} {v[1]:.4f} {v[2]:.4f}")

    # overlay: 4-class material map, translucent color on the grayscale original
    mat_cmap = ListedColormap([(0, 0, 0, 0),            # air  -> transparent
                               (1.0, 0.85, 0.1, 0.5),   # adipose -> yellow
                               (0.2, 0.8, 1.0, 0.45),   # muscle -> cyan
                               (1.0, 0.3, 0.2, 0.6)])   # bone -> red
    slices = [("axial z=150", vol[150], label[150], dens[150]),
              ("coronal y=150", vol[:, 150], label[:, 150], dens[:, 150]),
              ("sagittal x=128", vol[:, :, 128], label[:, :, 128], dens[:, :, 128])]
    fig, ax = plt.subplots(3, 3, figsize=(11, 11))
    for r, (name, sl, lab, dn) in enumerate(slices):
        lo, hi = np.percentile(sl, [40, 99.5])
        ax[r, 0].imshow(sl, cmap="gray", vmin=lo, vmax=hi); ax[r, 0].set_title(f"original — {name}")
        ax[r, 1].imshow(sl, cmap="gray", vmin=lo, vmax=hi)
        ax[r, 1].imshow(lab, cmap=mat_cmap, vmin=0, vmax=3, interpolation="nearest")
        ax[r, 1].set_title("materials: adipose=yellow, muscle=cyan, bone=red")
        im = ax[r, 2].imshow(dn, cmap="viridis", vmin=0, vmax=1.92); ax[r, 2].set_title("density [g/cm³]")
        for c in range(3):
            ax[r, c].axis("off")
    fig.colorbar(im, ax=ax[:, 2], shrink=0.5, label="g/cm³")
    fig.suptitle("RabbitCT re-materialization: HU→(material, density) for poly-energy projection",
                 fontsize=13)
    outdir = str(root / "results" / "rabbit"); os.makedirs(outdir, exist_ok=True)
    fig.savefig(f"{outdir}/rabbit_materials.png", dpi=120, bbox_inches="tight"); plt.close(fig)
    print("wrote", f"{outdir}/rabbit_materials.png")


if __name__ == "__main__":
    main()
