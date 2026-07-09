"""HARD GATE: one cone-beam FDK reconstruction of a uniform value-1 sphere.

Replicates ConeBeamReconstructionExample.java geometry/order but builds the
CircularTrajectory programmatically (headless, no XML). Instruments every stage:
forward-projection sum, filtered sum, backprojection sum, interior DC.
"""
from __future__ import annotations
import os, sys
import numpy as np
import jpype

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import conrad_backend as cb

OUTDIR = "/private/tmp/claude-501/-Users-maier-Documents-SPIONvsXRay/f77329d1-6870-4df8-b0b7-dd7a968fa1bb/scratchpad"

# fixed geometry for the gate
SID_MM = 750.0
SDD_MM = 1200.0
N_VIEWS = 360
DET_LEN_MM = 400.0
FOV_MM = 128.0
VOXEL_MM = 1.0
DET_PIX_MM = 1.0
SPHERE_R_MM = 40.0


def _cls(pkg, name):
    return cb.class_getter(pkg).__getattr__(name)


def np3d_to_grid3d(vol):
    vol = np.ascontiguousarray(vol, dtype=np.float32)
    Z, Y, X = vol.shape
    G2 = _cls("edu.stanford.rsl.conrad.data.numeric", "Grid2D")
    G3 = _cls("edu.stanford.rsl.conrad.data.numeric", "Grid3D")
    g3 = G3(int(X), int(Y), int(Z))
    for z in range(Z):
        jbuf = jpype.JArray(jpype.JFloat)(vol[z].ravel(order="C"))
        g3.setSubGrid(z, G2(jbuf, int(X), int(Y)))
    return g3


def grid3d_to_np(g3):
    sz = g3.getSize()
    X, Y, Z = int(sz[0]), int(sz[1]), int(sz[2])
    out = np.empty((Z, Y, X), dtype=np.float32)
    for z in range(Z):
        g2 = g3.getSubGrid(z)
        out[z] = np.array(g2.getBuffer()[:], dtype=np.float32).reshape(Y, X)
    return out


def setup_geometry(voxel_mm, det_pix_mm):
    Configuration = _cls("edu.stanford.rsl.conrad.utils", "Configuration")
    CircularTrajectory = _cls("edu.stanford.rsl.conrad.geometry.trajectories", "CircularTrajectory")
    Projection = _cls("edu.stanford.rsl.conrad.geometry", "Projection")
    CameraAxisDirection = Projection.CameraAxisDirection
    SimpleVector = _cls("edu.stanford.rsl.conrad.numerics", "SimpleVector")

    Configuration.loadConfiguration()
    config = Configuration.getGlobalConfiguration()

    n_det = int(round(DET_LEN_MM / det_pix_mm))
    n_vox = int(round(FOV_MM / voxel_mm))

    traj = CircularTrajectory(config.getGeometry())
    rotationAxis = SimpleVector(0.0, 0.0, 1.0)
    traj.setDetectorWidth(int(n_det))
    traj.setDetectorHeight(int(n_det))
    traj.setSourceToAxisDistance(float(SID_MM))
    traj.setSourceToDetectorDistance(float(SDD_MM))
    traj.setReconDimensions(int(n_vox), int(n_vox), int(n_vox))
    traj.setReconDimensionX(int(n_vox))
    traj.setReconDimensionY(int(n_vox))
    traj.setReconDimensionZ(int(n_vox))
    traj.setOriginInPixelsX((n_vox - 1) / 2.0)
    traj.setOriginInPixelsY((n_vox - 1) / 2.0)
    traj.setOriginInPixelsZ((n_vox - 1) / 2.0)
    traj.setDetectorOffsetU(0.0)
    traj.setDetectorOffsetV(0.0)
    traj.setPixelDimensionX(float(det_pix_mm))
    traj.setPixelDimensionY(float(det_pix_mm))
    traj.setVoxelSpacingX(float(voxel_mm))
    traj.setVoxelSpacingY(float(voxel_mm))
    traj.setVoxelSpacingZ(float(voxel_mm))
    ang_inc = 360.0 / N_VIEWS
    traj.setAverageAngularIncrement(float(ang_inc))
    traj.setProjectionStackSize(int(N_VIEWS))
    traj.setDetectorUDirection(CameraAxisDirection.DETECTORMOTION_PLUS)
    traj.setDetectorVDirection(CameraAxisDirection.ROTATIONAXIS_PLUS)
    traj.setTrajectory(int(N_VIEWS), float(SID_MM), float(ang_inc), 0.0, 0.0,
                       CameraAxisDirection.DETECTORMOTION_PLUS,
                       CameraAxisDirection.ROTATIONAXIS_PLUS, rotationAxis)
    config.setGeometry(traj)
    Configuration.setGlobalConfiguration(config)

    # sanity: projection matrices actually populated?
    pm = traj.getProjectionMatrices()
    print(f"  n_det={n_det} n_vox={n_vox} projMats={len(pm) if pm is not None else None} "
          f"stackSize={traj.getProjectionStackSize()}")
    print(f"  P[0]=\n{np.array(pm[0].computeP().copyAsDoubleArray())}")
    return dict(n_det=n_det, n_vox=n_vox, det_pix_mm=det_pix_mm, voxel_mm=voxel_mm,
                focal=SDD_MM, maxU=n_det * det_pix_mm, maxV=n_det * det_pix_mm,
                deltaU=det_pix_mm, deltaV=det_pix_mm)


def make_sphere(n_vox, voxel_mm, radius_mm=SPHERE_R_MM):
    c = (n_vox - 1) / 2.0
    ax = (np.arange(n_vox) - c) * voxel_mm
    Z, Y, X = np.meshgrid(ax, ax, ax, indexing="ij")
    r2 = X * X + Y * Y + Z * Z
    return (r2 <= radius_mm * radius_mm).astype(np.float32)


def main():
    cb.setup(max_ram="8G")
    assert cb.opencl_available(), "OpenCL not available"
    print("=== GATE: single cone-beam FDK, value-1 sphere ===")
    geo = setup_geometry(VOXEL_MM, DET_PIX_MM)
    vol = make_sphere(geo["n_vox"], VOXEL_MM)
    print(f"  phantom voxels>0: {int((vol>0).sum())}  sum={vol.sum():.1f}")

    CBP = _cls("edu.stanford.rsl.tutorial.cone", "ConeBeamProjector")
    CBBP = _cls("edu.stanford.rsl.tutorial.cone", "ConeBeamBackprojector")
    CosF = _cls("edu.stanford.rsl.tutorial.cone", "ConeBeamCosineFilter")
    Ramp = _cls("edu.stanford.rsl.tutorial.filters", "RamLakKernel")

    OG3 = _cls("edu.stanford.rsl.conrad.data.numeric.opencl", "OpenCLGrid3D")
    G3 = _cls("edu.stanford.rsl.conrad.data.numeric", "Grid3D")

    grid = np3d_to_grid3d(vol)
    grid.setSpacing(VOXEL_MM, VOXEL_MM, VOXEL_MM)

    cbp = CBP()
    # IMPORTANT: the convenience projectRayDrivenCL(Grid3D) leaves the volume
    # un-uploaded on this Apple GPU -> all-zero sinogram (the prior-harness bug).
    # Use the fast OpenCLGrid3D path which forces host->device upload.
    nd = int(geo["n_det"])
    gridCL = OG3(grid); gridCL.getDelegate().prepareForDeviceOperation()
    sino = OG3(G3(nd, nd, N_VIEWS)); sino.getDelegate().prepareForDeviceOperation()
    cbp.fastProjectRayDrivenCL(sino, gridCL)
    sino.getDelegate().prepareForHostOperation()
    sino_np = grid3d_to_np(sino)  # (nProj, v, u)
    print(f"[stage1 forward]  sino shape={sino_np.shape} sum={sino_np.sum():.3e} "
          f"max={sino_np.max():.3e}")
    np.save(os.path.join(OUTDIR, "gate_sino.npy"), sino_np[N_VIEWS // 2])

    cbFilter = CosF(float(geo["focal"]), float(geo["maxU"]), float(geo["maxV"]),
                    float(geo["deltaU"]), float(geo["deltaV"]))
    ramK = Ramp(int(geo["n_det"]), float(geo["deltaU"]))
    n_proj = int(sino.getSize()[2])
    maxV = int(geo["n_det"])
    for i in range(n_proj):
        sub = sino.getSubGrid(i)
        cbFilter.applyToGrid(sub)
        for j in range(maxV):
            ramK.applyToGrid(sub.getSubGrid(j))
    filt_np = grid3d_to_np(sino)
    print(f"[stage2 filtered] sum={filt_np.sum():.3e} absmax={np.abs(filt_np).max():.3e}")

    cbbp = CBBP()
    rec = cbbp.backprojectPixelDrivenCL(sino)
    rec_np = grid3d_to_np(rec)  # (Z,Y,X)
    print(f"[stage3 backproj] rec shape={rec_np.shape} sum={rec_np.sum():.3e} "
          f"min={rec_np.min():.3e} max={rec_np.max():.3e}")

    # interior DC on central slice
    n = geo["n_vox"]; zc = n // 2
    sl = rec_np[zc]
    c = (n - 1) / 2.0
    yy, xx = np.mgrid[0:n, 0:n]
    r_mm = np.sqrt((xx - c) ** 2 + (yy - c) ** 2) * VOXEL_MM
    interior = r_mm < 0.6 * SPHERE_R_MM
    dc = float(sl[interior].mean())
    print(f"[GATE RESULT] central-slice interior DC (r<0.6R) = {dc:.5f}  (target ~1.0)")

    np.save(os.path.join(OUTDIR, "gate_recon_central.npy"), sl)
    np.save(os.path.join(OUTDIR, "gate_recon_zprofile.npy"),
            np.array([rec_np[z][interior].mean() for z in range(n)]))

    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    im0 = ax[0].imshow(sino_np[N_VIEWS // 2], cmap="gray"); ax[0].set_title("forward sino (view 180)")
    plt.colorbar(im0, ax=ax[0], fraction=0.046)
    im1 = ax[1].imshow(sl, cmap="gray"); ax[1].set_title(f"recon central slice\nDC={dc:.4f}")
    plt.colorbar(im1, ax=ax[1], fraction=0.046)
    ax[2].plot(sl[n // 2]); ax[2].axhline(1.0, color="r", ls=":"); ax[2].set_title("central row profile")
    ax[2].grid(alpha=0.3)
    fig.suptitle(f"GATE: cone-beam FDK, value-1 sphere R={SPHERE_R_MM}mm, "
                 f"{n}^3 @ {VOXEL_MM}mm, DC={dc:.4f}")
    fig.tight_layout()
    p = os.path.join(OUTDIR, "gate_recon.png")
    fig.savefig(p, dpi=110); plt.close(fig)
    print("wrote", p)

    if dc > 0.1:
        print("\n>>> GATE PASSED: reconstruction is non-empty and near unity.\n")
    else:
        print("\n>>> GATE FAILED: reconstruction empty/too low.\n")


if __name__ == "__main__":
    main()
