"""3D cone-beam FDK resolution study on CONRAD's OpenCL projector/backprojector.

CHARACTERIZATION ONLY. No source is modified. This harness measures how an FDK
reconstruction of a uniform sphere depends on (recon VOXEL size) x (detector
PIXEL size) with the physical object and physical detector held FIXED, using
CONRAD's *cone-beam* tutorial classes exclusively:

  forward     edu.stanford.rsl.tutorial.cone.ConeBeamProjector.projectRayDrivenCL(Grid3D)
  cosine      edu.stanford.rsl.tutorial.cone.ConeBeamCosineFilter (FDK cos weight)
  ramp        edu.stanford.rsl.tutorial.filters.{RamLakKernel,SheppLoganKernel}
  backproject edu.stanford.rsl.tutorial.cone.ConeBeamBackprojector.backprojectPixelDrivenCL(Grid3D)

Geometry is set up PROGRAMMATICALLY on the global Configuration's CircularTrajectory
(the pattern from IterativeReconstructionTestA + ConeBeamReconstructionExample);
both CL classes read the global Configuration.getGeometry() in their constructor,
so we rebuild the global geometry and then construct fresh projector/backprojector
per (voxel, detector-pixel) cell.

Outputs (PNGs) -> the scratchpad dir passed as OUTDIR below.
"""
from __future__ import annotations
import os
import sys
import numpy as np
import jpype

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import conrad_backend

OUTDIR = "/private/tmp/claude-501/-Users-maier-Documents-SPIONvsXRay/f77329d1-6870-4df8-b0b7-dd7a968fa1bb/scratchpad"

# ---------------- fixed physical geometry (held constant across the sweep) ----
SID_MM = 750.0          # source -> isocenter (axis) distance
SDD_MM = 1200.0         # source -> detector distance
N_VIEWS = 360           # projections over 2*pi
DET_LEN_MM = 400.0      # PHYSICAL detector length (u and v), held fixed [mm]
FOV_MM = 128.0          # PHYSICAL reconstructed FOV (cube edge), held fixed [mm]
SPHERE_R_MM = 40.0      # phantom sphere radius, held fixed [mm]

# ---------------- sweeps -------------------------------------------------------
VOXELS_MM = [4.0, 2.0, 1.0, 0.5]        # recon voxel; FOV fixed -> voxel count = FOV/voxel
DET_PIX_MM = [2.0, 1.0, 0.5, 0.25]      # detector pixel; det length fixed -> count = LEN/pix
# FOV 128 mm: 32^3 @4mm, 64^3 @2mm, 128^3 @1mm, 256^3 @0.5mm.
# Detector 400 mm: 200 @2mm, 400 @1mm, 800 @0.5mm, 1600 @0.25mm.


def _cls(pkg, name):
    conrad_backend.setup()
    return conrad_backend.class_getter(pkg).__getattr__(name)


def np3d_to_grid3d(vol):
    """numpy (Z,Y,X) float -> CONRAD Grid3D (built subgrid by subgrid, X fastest).

    CONRAD Grid3D is a stack of Grid2D(width=X, height=Y); depth = Z. We build
    each Grid2D from a row-major float[] (X fastest within a row)."""
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
    """CONRAD Grid3D -> numpy (Z,Y,X)."""
    sz = g3.getSize()
    X, Y, Z = int(sz[0]), int(sz[1]), int(sz[2])
    out = np.empty((Z, Y, X), dtype=np.float32)
    for z in range(Z):
        g2 = g3.getSubGrid(z)
        out[z] = np.array(g2.getBuffer()[:], dtype=np.float32).reshape(Y, X)
    return out


def setup_geometry(voxel_mm, det_pix_mm):
    """Program the global Configuration's CircularTrajectory for this cell.

    Physical detector length DET_LEN_MM and physical FOV_MM are held fixed:
        detector count = DET_LEN_MM / det_pix_mm  (u and v)
        recon count    = FOV_MM     / voxel_mm    (x, y, z)
    Returns dict with the numbers the FDK cosine/ramp filters need.
    """
    Configuration = _cls("edu.stanford.rsl.conrad.utils", "Configuration")
    CircularTrajectory = _cls("edu.stanford.rsl.conrad.geometry.trajectories",
                              "CircularTrajectory")
    Projection = _cls("edu.stanford.rsl.conrad.geometry", "Projection")
    CameraAxisDirection = Projection.CameraAxisDirection
    SimpleVector = _cls("edu.stanford.rsl.conrad.numerics", "SimpleVector")

    Configuration.loadConfiguration()
    config = Configuration.getGlobalConfiguration()

    n_det = int(round(DET_LEN_MM / det_pix_mm))     # detector width = height
    n_vox = int(round(FOV_MM / voxel_mm))           # recon cube edge

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
    # full 360 deg circular scan
    traj.setTrajectory(int(N_VIEWS), float(SID_MM), float(ang_inc), 0.0, 0.0,
                       CameraAxisDirection.DETECTORMOTION_PLUS,
                       CameraAxisDirection.ROTATIONAXIS_PLUS, rotationAxis)

    config.setGeometry(traj)
    Configuration.setGlobalConfiguration(config)

    return dict(n_det=n_det, n_vox=n_vox, det_pix_mm=det_pix_mm, voxel_mm=voxel_mm,
                focal=SDD_MM, maxU=n_det * det_pix_mm, maxV=n_det * det_pix_mm,
                deltaU=det_pix_mm, deltaV=det_pix_mm)


def make_sphere(n_vox, voxel_mm, radius_mm=SPHERE_R_MM):
    """Uniform sphere value 1.0, centered, radius fixed in mm."""
    c = (n_vox - 1) / 2.0
    ax = (np.arange(n_vox) - c) * voxel_mm
    Z, Y, X = np.meshgrid(ax, ax, ax, indexing="ij")
    r2 = X * X + Y * Y + Z * Z
    return (r2 <= radius_mm * radius_mm).astype(np.float32)


def fdk(vol_np, geo, kernel="ramlak"):
    """CONRAD cone-beam FDK: project -> cosine -> ramp -> backproject (all CL/CONRAD)."""
    CBP = _cls("edu.stanford.rsl.tutorial.cone", "ConeBeamProjector")
    CBBP = _cls("edu.stanford.rsl.tutorial.cone", "ConeBeamBackprojector")
    CosF = _cls("edu.stanford.rsl.tutorial.cone", "ConeBeamCosineFilter")
    if kernel == "ramlak":
        Ramp = _cls("edu.stanford.rsl.tutorial.filters", "RamLakKernel")
    else:
        Ramp = _cls("edu.stanford.rsl.tutorial.filters", "SheppLoganKernel")

    OG3 = _cls("edu.stanford.rsl.conrad.data.numeric.opencl", "OpenCLGrid3D")
    G3 = _cls("edu.stanford.rsl.conrad.data.numeric", "Grid3D")

    grid = np3d_to_grid3d(vol_np)
    grid.setSpacing(geo["voxel_mm"], geo["voxel_mm"], geo["voxel_mm"])

    cbp = CBP()
    # NOTE: the convenience projectRayDrivenCL(Grid3D) leaves the volume texture
    # un-uploaded on this Apple GPU and returns an ALL-ZERO sinogram (this was the
    # prior harness's bug: stage-1 forward projection was empty). The fast
    # OpenCLGrid3D path calls prepareForDeviceOperation() on the volume, forcing
    # the host->device upload before the projection kernel runs.
    nd = int(geo["n_det"])
    n_proj = int(N_VIEWS)
    gridCL = OG3(grid); gridCL.getDelegate().prepareForDeviceOperation()
    sino = OG3(G3(nd, nd, n_proj)); sino.getDelegate().prepareForDeviceOperation()
    cbp.fastProjectRayDrivenCL(sino, gridCL)
    sino.getDelegate().prepareForHostOperation()

    # FDK cosine weight (per projection) then ramp along u (per detector row v)
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

    cbbp = CBBP()
    rec = cbbp.backprojectPixelDrivenCL(sino)     # Grid3D volume
    return grid3d_to_np(rec)


def measure(rec, geo):
    """DC (mean interior r<0.6R central slice), center/edge ratio, axial shading."""
    n = geo["n_vox"]
    vox = geo["voxel_mm"]
    zc = n // 2
    sl = rec[zc]                                   # central axial slice (Y,X)
    c = (n - 1) / 2.0
    yy, xx = np.mgrid[0:n, 0:n]
    r_mm = np.sqrt((xx - c) ** 2 + (yy - c) ** 2) * vox
    interior = r_mm < 0.6 * SPHERE_R_MM
    dc = float(sl[interior].mean())
    # center value (small central disk) vs edge value (annulus just inside R)
    center = float(sl[r_mm < 0.15 * SPHERE_R_MM].mean())
    edge_band = (r_mm > 0.75 * SPHERE_R_MM) & (r_mm < 0.9 * SPHERE_R_MM)
    edge = float(sl[edge_band].mean())
    ce = center / edge if edge != 0 else float("nan")
    # axial (cone-angle) shading: DC of interior at several z slices
    zprof = []
    for z in range(n):
        s = rec[z]
        m = interior & (r_mm < 0.6 * SPHERE_R_MM)
        zprof.append(float(s[m].mean()))
    return dict(dc=dc, center=center, edge=edge, ce=ce,
                slice=sl, zprof=np.array(zprof))


def main():
    conrad_backend.setup(max_ram="12G")
    if not conrad_backend.opencl_available():
        print("OpenCL NOT available -- aborting (this study requires the CL pipeline).")
        sys.exit(1)
    print("OpenCL available. Running cone-beam FDK resolution sweep.")
    print(f"Fixed: SID={SID_MM} SDD={SDD_MM} views={N_VIEWS} det_len={DET_LEN_MM}mm "
          f"FOV={FOV_MM}mm sphere_R={SPHERE_R_MM}mm")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    kernels = ["ramlak", "shepplogan"]
    results = {k: {} for k in kernels}     # results[k][(voxel, detpix)] = measure dict
    failures = []

    for kernel in kernels:
        for voxel_mm in VOXELS_MM:
            for det_pix_mm in DET_PIX_MM:
                key = (voxel_mm, det_pix_mm)
                try:
                    geo = setup_geometry(voxel_mm, det_pix_mm)
                    vol = make_sphere(geo["n_vox"], voxel_mm)
                    rec = fdk(vol, geo, kernel=kernel)
                    m = measure(rec, geo)
                    m["geo"] = geo
                    results[kernel][key] = m
                    print(f"[{kernel}] voxel={voxel_mm} det={det_pix_mm} "
                          f"n_vox={geo['n_vox']} n_det={geo['n_det']} "
                          f"DC={m['dc']:.4f} c/e={m['ce']:.4f} "
                          f"center={m['center']:.4f} edge={m['edge']:.4f}")
                except Exception as e:
                    msg = repr(e)
                    low = msg.lower()
                    is_oom = any(t in low for t in (
                        "outofmemory", "out of memory", "cl_mem_object_allocation",
                        "cl_out_of_resources", "cl_invalid_buffer_size",
                        "mem_object_allocation", "heap space", "cannot allocate",
                        "clenqueue", "allocation failure", "cl_out_of_host_memory"))
                    tag = "OOM" if is_oom else "FAIL"
                    print(f"[{kernel}] voxel={voxel_mm} det={det_pix_mm} {tag}: {e}")
                    failures.append((kernel, key, tag, msg))
                    results[kernel][("status", key)] = tag

    # ---- image grids (rows=voxel, cols=detector pixel), shared window per kernel
    for kernel in kernels:
        cells = results[kernel]
        if not cells:
            continue
        real = [m for k, m in cells.items()
                if not (isinstance(k, tuple) and k and k[0] == "status")]
        if not real:
            continue
        vmax = max(m["slice"].max() for m in real) or 1.0
        nrow, ncol = len(VOXELS_MM), len(DET_PIX_MM)
        fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 4 * nrow),
                                 squeeze=False)
        for r, voxel_mm in enumerate(VOXELS_MM):
            for cix, det_pix_mm in enumerate(DET_PIX_MM):
                ax = axes[r][cix]
                key = (voxel_mm, det_pix_mm)
                if key in cells:
                    m = cells[key]
                    ax.imshow(m["slice"], cmap="gray", vmin=0, vmax=vmax)
                    ax.set_title(f"vox={voxel_mm}mm det={det_pix_mm}mm\n"
                                 f"DC={m['dc']:.3f} c/e={m['ce']:.3f}", fontsize=10)
                else:
                    tag = cells.get(("status", key), "FAIL")
                    ax.text(0.5, 0.5, tag, ha="center", va="center",
                            fontsize=14, color="red")
                    ax.set_title(f"vox={voxel_mm}mm det={det_pix_mm}mm", fontsize=10)
                ax.axis("off")
        fig.suptitle(f"Cone-beam FDK central axial slice -- {kernel} "
                     f"(shared window 0..{vmax:.3f})", fontsize=13)
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        p = os.path.join(OUTDIR, f"cone_res_grid_{kernel}.png")
        fig.savefig(p, dpi=110)
        plt.close(fig)
        print("wrote", p)

    # ---- line plots: DC vs det pix (one line per voxel), c/e vs det pix; both kernels
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    for col, kernel in enumerate(kernels):
        cells = results[kernel]
        axdc, axce = axes[0][col], axes[1][col]
        for voxel_mm in VOXELS_MM:
            xs, dcs, ces = [], [], []
            for det_pix_mm in DET_PIX_MM:
                key = (voxel_mm, det_pix_mm)
                if key in cells:
                    xs.append(det_pix_mm)
                    dcs.append(cells[key]["dc"])
                    ces.append(cells[key]["ce"])
            if xs:
                axdc.plot(xs, dcs, "o-", label=f"voxel={voxel_mm}mm")
                axce.plot(xs, ces, "s-", label=f"voxel={voxel_mm}mm")
        axdc.set_title(f"{kernel}: DC (interior mean) vs detector pixel")
        axdc.set_xlabel("detector pixel [mm]"); axdc.set_ylabel("DC")
        axdc.axhline(1.0, color="k", ls=":", lw=0.8); axdc.legend(); axdc.grid(alpha=0.3)
        axce.set_title(f"{kernel}: center/edge ratio vs detector pixel")
        axce.set_xlabel("detector pixel [mm]"); axce.set_ylabel("center/edge")
        axce.axhline(1.0, color="k", ls=":", lw=0.8); axce.legend(); axce.grid(alpha=0.3)
    fig.tight_layout()
    p = os.path.join(OUTDIR, "cone_res_lines.png")
    fig.savefig(p, dpi=120); plt.close(fig)
    print("wrote", p)

    # ---- horizontal center profiles (central slice, through center row)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for col, kernel in enumerate(kernels):
        cells = results[kernel]
        ax = axes[col]
        for key, m in cells.items():
            if isinstance(key, tuple) and key and key[0] == "status":
                continue
            voxel_mm, det_pix_mm = key
            sl = m["slice"]
            n = sl.shape[0]
            row = sl[n // 2]
            xmm = (np.arange(n) - (n - 1) / 2.0) * voxel_mm
            ax.plot(xmm, row, label=f"vox={voxel_mm} det={det_pix_mm}")
        ax.axvline(-SPHERE_R_MM, color="k", ls=":", lw=0.8)
        ax.axvline(SPHERE_R_MM, color="k", ls=":", lw=0.8)
        ax.axhline(1.0, color="gray", ls=":", lw=0.8)
        ax.set_title(f"{kernel}: horizontal center profile (central slice)")
        ax.set_xlabel("x [mm]"); ax.set_ylabel("recon value")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout()
    p = os.path.join(OUTDIR, "cone_res_profiles.png")
    fig.savefig(p, dpi=120); plt.close(fig)
    print("wrote", p)

    # ---- summary tables
    def table(metric):
        lines = []
        header = "voxel\\det | " + " | ".join(f"{d:>8.2f}" for d in DET_PIX_MM)
        lines.append(header)
        for voxel_mm in VOXELS_MM:
            row = [f"{voxel_mm:>8.2f} "]
            for det_pix_mm in DET_PIX_MM:
                key = (voxel_mm, det_pix_mm)
                for kernel in kernels:
                    pass
            lines.append(row)
        return lines

    print("\n================ SUMMARY ================")
    for kernel in kernels:
        cells = results[kernel]
        print(f"\n--- {kernel} : DC (interior mean) ---")
        print("           det=" + "  ".join(f"{d:.2f}mm" for d in DET_PIX_MM))
        for voxel_mm in VOXELS_MM:
            vals = []
            for det_pix_mm in DET_PIX_MM:
                m = cells.get((voxel_mm, det_pix_mm))
                if m:
                    vals.append(f"{m['dc']:>7.4f}")
                else:
                    vals.append(f"{cells.get(('status', (voxel_mm, det_pix_mm)), 'FAIL'):>7}")
            print(f"  vox={voxel_mm:.2f}mm  " + "  ".join(vals))
        print(f"--- {kernel} : center/edge ---")
        for voxel_mm in VOXELS_MM:
            vals = []
            for det_pix_mm in DET_PIX_MM:
                m = cells.get((voxel_mm, det_pix_mm))
                if m:
                    vals.append(f"{m['ce']:>7.4f}")
                else:
                    vals.append(f"{cells.get(('status', (voxel_mm, det_pix_mm)), 'FAIL'):>7}")
            print(f"  vox={voxel_mm:.2f}mm  " + "  ".join(vals))
        # axial shading note
        for voxel_mm in VOXELS_MM:
            for det_pix_mm in DET_PIX_MM:
                m = cells.get((voxel_mm, det_pix_mm))
                if m:
                    zp = m["zprof"]
                    nz = len(zp)
                    q = nz // 4
                    print(f"    [{kernel} vox={voxel_mm} det={det_pix_mm}] "
                          f"axial DC z=25%:{zp[q]:.4f} mid:{zp[nz//2]:.4f} "
                          f"z=75%:{zp[3*q]:.4f} (edge/mid={zp[q]/zp[nz//2]:.3f})")

    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(" ", f)


if __name__ == "__main__":
    main()
