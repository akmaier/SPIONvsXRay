"""E5 -- 3D cone-beam FDK verification, POLY-ENERGETIC and MATERIAL-SEPARATED.

Mirrors the 2D spectral detectability pipeline (src/run_detectability.py) but in a
full 3D cone-beam geometry, in real rabbit anatomy, at the E2 optimum spectrum
(results/spectral/optimum.json: 70 kVp / Al 5.0, PCD bins 10/37.5/47.5/70 keV).

Pipeline
--------
1. Anatomy: the RabbitCT reference reconstruction (256^3) re-materialized by
   src/rabbit_materials.segment into base-material occupancy volumes
   (adipose / muscle / bone, each scaled by its per-voxel density -> effective
   "grams-path"), plus an iron-loaded spherical tumor inserted in soft tissue.
2. Base-material 3D cone forward projection (CONRAD ConeBeamProjector,
   fastProjectRayDrivenCL(OpenCLGrid3D,OpenCLGrid3D) after prepareForDeviceOperation;
   the plain projectRayDrivenCL(Grid3D) returns all-zero on this Apple/JOCL setup) ->
   one path-integral sinogram per base material [g/cm^2] and one for the iron
   contrast [g Fe/cm^2].
3. Poly-energetic combine at the optimum spectrum, reusing the 2D spectral model:
     EID = energy-weighted single-log (with EID water precorrection)
     PCD = post-combination matched filter: bin counts summed in the I0 domain with
           CNR-optimal weights, single log, ONE water precorrection on the combined
           effective spectrum (bh_poly_for). Identical math to run_detectability.
4. FDK reconstruct: ConeBeamCosineFilter -> RamLak ramp -> ConeBeamBackprojector.
5. 4-5 runs varying the tumor iron level (SPION I fresh, cellularity sweep).

Volume 256^3, detector 1 mm (0.25 mm overflows the JVM 32-bit buffer). Recon slices
(EID vs PCD, a tumor slice) -> results/3d/; per-run tumor DeltaHU + CNR reported.

USAGE:  python scripts/run_e5_3d.py [--sanity] [--reps N] [--views V]
"""
from __future__ import annotations
import os
import sys
import argparse
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import conrad_backend
import conrad_ct
import rabbitct
import rabbit_materials as rm
import materials
import spectrum as spec

import jpype

ROOT = conrad_backend.REPO_ROOT
RCTD = str(ROOT / "data" / "rabbitct" / "rabbitct_512-v2.rctd")
VOLF = str(ROOT / "data" / "rabbitct" / "reference_256.vol")
OUTDIR = str(ROOT / "results" / "3d")

# ---- fixed cone geometry (mirrors scripts/cone_resolution_study.py) ----------
SID_MM = 750.0
SDD_MM = 1200.0
VOXEL_MM = 1.2          # 256^3 volume -> FOV = 307 mm (covers the rabbit body)
DET_PIX_MM = 1.0        # 0.25 mm overflows the JVM 32-bit buffer
N_DET = 512             # 512 mm detector; magnified object 307*1.6 ~ 491 mm < 512
N_VOX = 256
PCD_BINS = None         # filled from optimum.json

# ---- tumor / iron sweep ------------------------------------------------------
# SPION I fresh (0 h) = 8.23 pg Fe/cell (config.CELLULAR_LOADING "I_113").
# c_Fe [mg Fe/ml] = pg_fe * density[cells/cm^3] * 1e-9.
PG_FE_SPION_I_FRESH = 8.23
TUMOR_RADIUS_MM = 12.0
TUMOR_CENTER_VOX = (150, 150, 128)   # (z, y, x) neck soft tissue (rabbit3d.py)


def _cls(pkg, name):
    return conrad_backend.class_getter(pkg).__getattr__(name)


# ---------------------------------------------------------------------------
# Grid <-> numpy (from cone_resolution_study.py)
# ---------------------------------------------------------------------------
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


def setup_geometry(n_views):
    """Program the global CircularTrajectory for the fixed cone geometry."""
    Configuration = _cls("edu.stanford.rsl.conrad.utils", "Configuration")
    CircularTrajectory = _cls("edu.stanford.rsl.conrad.geometry.trajectories",
                              "CircularTrajectory")
    Projection = _cls("edu.stanford.rsl.conrad.geometry", "Projection")
    CameraAxisDirection = Projection.CameraAxisDirection
    SimpleVector = _cls("edu.stanford.rsl.conrad.numerics", "SimpleVector")

    Configuration.loadConfiguration()
    config = Configuration.getGlobalConfiguration()

    traj = CircularTrajectory(config.getGeometry())
    rotationAxis = SimpleVector(0.0, 0.0, 1.0)
    traj.setDetectorWidth(int(N_DET))
    traj.setDetectorHeight(int(N_DET))
    traj.setSourceToAxisDistance(float(SID_MM))
    traj.setSourceToDetectorDistance(float(SDD_MM))
    traj.setReconDimensions(int(N_VOX), int(N_VOX), int(N_VOX))
    traj.setReconDimensionX(int(N_VOX))
    traj.setReconDimensionY(int(N_VOX))
    traj.setReconDimensionZ(int(N_VOX))
    traj.setOriginInPixelsX((N_VOX - 1) / 2.0)
    traj.setOriginInPixelsY((N_VOX - 1) / 2.0)
    traj.setOriginInPixelsZ((N_VOX - 1) / 2.0)
    traj.setDetectorOffsetU(0.0)
    traj.setDetectorOffsetV(0.0)
    traj.setPixelDimensionX(float(DET_PIX_MM))
    traj.setPixelDimensionY(float(DET_PIX_MM))
    traj.setVoxelSpacingX(float(VOXEL_MM))
    traj.setVoxelSpacingY(float(VOXEL_MM))
    traj.setVoxelSpacingZ(float(VOXEL_MM))
    ang_inc = 360.0 / n_views
    traj.setAverageAngularIncrement(float(ang_inc))
    traj.setProjectionStackSize(int(n_views))
    traj.setDetectorUDirection(CameraAxisDirection.DETECTORMOTION_PLUS)
    traj.setDetectorVDirection(CameraAxisDirection.ROTATIONAXIS_PLUS)
    traj.setTrajectory(int(n_views), float(SID_MM), float(ang_inc), 0.0, 0.0,
                       CameraAxisDirection.DETECTORMOTION_PLUS,
                       CameraAxisDirection.ROTATIONAXIS_PLUS, rotationAxis)
    config.setGeometry(traj)
    Configuration.setGlobalConfiguration(config)
    return dict(n_det=N_DET, n_vox=N_VOX, focal=SDD_MM,
                maxU=N_DET * DET_PIX_MM, maxV=N_DET * DET_PIX_MM,
                deltaU=DET_PIX_MM, deltaV=DET_PIX_MM, voxel_mm=VOXEL_MM,
                n_views=n_views)


def cone_forward(vol_np, geo):
    """3D cone-beam forward projection -> sinogram (Z=view, Y=v, X=u) numpy.

    Uses fastProjectRayDrivenCL(OpenCLGrid3D, OpenCLGrid3D) after
    prepareForDeviceOperation (the plain projectRayDrivenCL returns all-zero on
    this Apple/JOCL setup -- see cone_resolution_study.py). Returns path integral
    in the volume's spacing units: value * voxel_mm gives the mm line integral of
    the input occupancy volume.
    """
    CBP = _cls("edu.stanford.rsl.tutorial.cone", "ConeBeamProjector")
    OG3 = _cls("edu.stanford.rsl.conrad.data.numeric.opencl", "OpenCLGrid3D")
    G3 = _cls("edu.stanford.rsl.conrad.data.numeric", "Grid3D")

    grid = np3d_to_grid3d(vol_np)
    grid.setSpacing(geo["voxel_mm"], geo["voxel_mm"], geo["voxel_mm"])
    cbp = CBP()
    nd = int(geo["n_det"]); nv = int(geo["n_views"])
    gridCL = OG3(grid); gridCL.getDelegate().prepareForDeviceOperation()
    sino = OG3(G3(nd, nd, nv)); sino.getDelegate().prepareForDeviceOperation()
    cbp.fastProjectRayDrivenCL(sino, gridCL)
    sino.getDelegate().prepareForHostOperation()
    out = grid3d_to_np(sino)
    del gridCL, sino, grid
    return out


def fdk_from_sino(sino_np, geo):
    """FDK: cosine weight -> RamLak ramp -> cone backprojection (all CONRAD)."""
    CBBP = _cls("edu.stanford.rsl.tutorial.cone", "ConeBeamBackprojector")
    CosF = _cls("edu.stanford.rsl.tutorial.cone", "ConeBeamCosineFilter")
    Ramp = _cls("edu.stanford.rsl.tutorial.filters", "RamLakKernel")
    OG3 = _cls("edu.stanford.rsl.conrad.data.numeric.opencl", "OpenCLGrid3D")

    sino = np3d_to_grid3d(sino_np)
    cbFilter = CosF(float(geo["focal"]), float(geo["maxU"]), float(geo["maxV"]),
                    float(geo["deltaU"]), float(geo["deltaV"]))
    ramK = Ramp(int(geo["n_det"]), float(geo["deltaU"]))
    n_proj = int(sino.getSize()[2]); maxV = int(geo["n_det"])
    for i in range(n_proj):
        sub = sino.getSubGrid(i)
        cbFilter.applyToGrid(sub)
        for j in range(maxV):
            ramK.applyToGrid(sub.getSubGrid(j))
    sinoCL = OG3(sino); sinoCL.getDelegate().prepareForDeviceOperation()
    cbbp = CBBP()
    rec = cbbp.backprojectPixelDrivenCL(sinoCL)
    return grid3d_to_np(rec)


# ---------------------------------------------------------------------------
# Anatomy: base-material occupancy volumes + iron tumor
# ---------------------------------------------------------------------------
def build_base_volumes():
    """Re-materialize the RabbitCT reference volume into per-material 'grams-path'
    occupancy volumes (density-weighted) and a tumor mask.

    For a poly-energetic line integral we need, per material m, the projection of
    rho_m(x)  [g/cm^3] so that tau(E) = sum_m (mu/rho)_m(E) * proj(rho_m).
    We build one volume per material holding that material's per-voxel density where
    the voxel belongs to that class (0 elsewhere). The cone projector returns the
    path integral of that volume; scaling by voxel_mm * 0.1 converts to [g/cm^2].
    """
    vol = rabbitct.load_reference_volume(VOLF, 256)
    seg = rm.segment(vol)
    label, dens, bfrac = seg["label"], seg["density"], seg["bone_frac"]

    # adipose (1), muscle (2), bone(3, blended muscle<->cortical by bone_frac).
    # For the bone class we split the density into a soft (muscle) part and a bone
    # part by bone_frac so both mass-attenuation curves apply where appropriate.
    v_adip = np.where(label == 1, dens, 0.0).astype(np.float32)
    v_musc = np.where(label == 2, dens, 0.0).astype(np.float32)
    # bone-class voxel: fraction bfrac is cortical bone, (1-bfrac) muscle-like
    bmask = label == 3
    v_bone = np.zeros_like(dens, np.float32)
    v_bone[bmask] = (dens[bmask] * bfrac[bmask]).astype(np.float32)
    v_musc[bmask] = (dens[bmask] * (1.0 - bfrac[bmask])).astype(np.float32)

    # soft-tissue mask for tumor placement (muscle class, away from bone/air)
    soft_mask = (label == 2)
    return dict(adipose=v_adip, muscle=v_musc, bone=v_bone), soft_mask, vol


def tumor_mask(soft_mask):
    """Spherical tumor inside soft tissue at TUMOR_CENTER_VOX."""
    L = N_VOX
    zz, yy, xx = np.mgrid[0:L, 0:L, 0:L]
    r = TUMOR_RADIUS_MM / VOXEL_MM
    cz, cy, cx = TUMOR_CENTER_VOX
    m = ((zz - cz) ** 2 + (yy - cy) ** 2 + (xx - cx) ** 2) < r ** 2
    m &= soft_mask
    return m


# ---------------------------------------------------------------------------
# Spectral model (mirrors run_detectability.polychromatic_accumulators etc.)
# ---------------------------------------------------------------------------
def load_optimum():
    import json
    with open(ROOT / "results" / "spectral" / "optimum.json") as f:
        opt = json.load(f)
    kvp = float(opt["kvp"])
    filters = [tuple(x) for x in opt["filters"]]
    edges = np.array(opt["pcd_bin_edges_kev"], float)
    return kvp, filters, edges, opt


def spectrum_at_optimum(kvp, filters):
    E, flux = spec.conrad_spectrum(kvp)
    if filters:
        flux = spec.apply_filters(E, flux, list(filters))
    s = flux / flux.sum()
    return E, s


def accumulators(base_paths_gcm2, iron_path_gcm2, E, s, edges, n0):
    """Noise-free detector accumulators over the 3D sinogram stack.

    base_paths_gcm2: dict material -> path integral [g/cm^2] (density-weighted).
    iron_path_gcm2:  iron-mass path integral [g Fe/cm^2] through the tumor.
    tau(E) = sum_m (mu/rho)_m(E) * path_m + oxide_contrast(E) * iron_path.
    Returns S_det, S_det_E2, S_air, C_det[b], C_air[b] (mirrors 2D model).
    """
    names = list(base_paths_gcm2.keys())
    # mass attenuation (mu/rho) [cm^2/g] per material
    marho = {m: materials.mass_attenuation(m, E) for m in names}
    ox = materials.oxide_contrast_massatten(E)          # cm^2 per g Fe
    shape = next(iter(base_paths_gcm2.values())).shape
    nb = len(edges) - 1

    S_det = np.zeros(shape); S_det_E2 = np.zeros(shape)
    S_air = float(np.sum(n0 * s * E))
    C_det = [np.zeros(shape) for _ in range(nb)]
    C_air = [0.0] * nb
    for j, Ej in enumerate(E):
        tau = np.zeros(shape)
        for m in names:
            tau += marho[m][j] * base_paths_gcm2[m]
        tau += ox[j] * iron_path_gcm2
        nph = n0 * s[j] * np.exp(-tau)
        S_det += nph * Ej
        S_det_E2 += nph * Ej * Ej
        b = int(np.searchsorted(edges, Ej, side="right") - 1)
        if 0 <= b < nb:
            C_det[b] += nph
            C_air[b] += n0 * s[j]
    return dict(S_det=S_det, S_det_E2=S_det_E2, S_air=S_air,
                C_det=C_det, C_air=C_air, E=E, s=s, edges=edges)


def pcd_weights(E, s, edges, body_diam_cm):
    """Count-domain matched-filter weights (run_detectability._pcd_weights)."""
    mu_tissue = materials.linear_attenuation_soft(E)
    ox = materials.oxide_contrast_massatten(E)
    Nt = s * np.exp(-mu_tissue * body_diam_cm)          # relative bkg counts/energy
    w = []
    for b in range(len(edges) - 1):
        m = (E >= edges[b]) & (E < edges[b + 1])
        S_b = float((Nt[m] * ox[m]).sum())
        V_b = float(Nt[m].sum())
        w.append(S_b / (V_b + 1e-12))
    w = np.clip(np.array(w), 0.0, None)
    return w / (w.sum() + 1e-12)


def bh_poly_for(acc, detector, w):
    """Water precorrection polynomial (run_detectability.bh_poly_for)."""
    E, s, edges = acc["E"], acc["s"], acc["edges"]
    mu_w = materials.linear_attenuation("water", E)
    if detector == "EID":
        return conrad_ct.water_precorrection_poly(E, s * E, mu_w)
    w_eff = np.zeros_like(s)
    for b in range(len(edges) - 1):
        m = (E >= edges[b]) & (E < edges[b + 1])
        w_eff[m] = w[b] * s[m]
    return conrad_ct.water_precorrection_poly(E, w_eff, mu_w)


def line_integral(acc, detector, seed, w, bh_polys, noiseless=False):
    """Poly-energetic combine -> line-integral sinogram (run_detectability.line_integral)."""
    rng = np.random.default_rng(seed)
    eps = 1e-6
    if detector == "EID":
        S = acc["S_det"]
        if not noiseless:
            S = S + rng.normal(0.0, np.sqrt(np.maximum(acc["S_det_E2"], 1e-30)))
        p = -np.log(np.clip(S, eps, None) / acc["S_air"])
        return np.polyval(bh_polys, p) if bh_polys is not None else p
    M = np.zeros_like(acc["S_det"]); M_air = 0.0
    for b, (cd, ca) in enumerate(zip(acc["C_det"], acc["C_air"])):
        cm = cd if noiseless else rng.poisson(np.maximum(cd, 0.0))
        M += w[b] * cm
        M_air += w[b] * ca
    p = -np.log(np.clip(M, eps, None) / max(M_air, eps))
    return np.polyval(bh_polys, p) if bh_polys is not None else p


# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------
def to_hu(rec, water_val):
    return 1000.0 * (rec - water_val) / (water_val + 1e-12)


def measure_iron_signal(rec, rec0, tmask, bg_shell):
    """Iron-only tumor signal, mirroring the 2D pipeline's iron_delta_hu.

    In real anatomy the tumor's own tissue differs structurally from any local
    background ROI, so a raw tumor-minus-shell difference is dominated by anatomy,
    not iron. Exactly as run_detectability measures the c0-CORRECTED iron signal
    (insert-with-iron minus the same insert at c0), we isolate iron as the tumor
    ROI mean of (iron recon - identical no-iron recon). HU calibration uses the
    soft-tissue background mean as the water/0-HU reference.
    """
    water = float(rec0[bg_shell].mean())                 # soft-tissue -> ~0 HU ref
    d_raw = float((rec - rec0)[tmask].mean())            # iron-only tumor signal
    d_hu = 1000.0 * d_raw / (water + 1e-12)
    return dict(delta_hu=d_hu, water=water, d_raw=d_raw)


def measure_noise(rec, bg_shell, water):
    """Quantum noise on the insert: std of the background soft-tissue ROI, in HU.

    The background shell is a homogeneous soft-tissue annulus, so its voxel-to-voxel
    std after the (noiseless) anatomy is subtracted is the reconstructed quantum
    noise. Measured on a NOISE-DIFFERENCED recon (noisy - noiseless) so residual
    anatomical texture cancels and only the noise remains, matching the 2D noise std.
    """
    noise_raw = float(rec[bg_shell].std())
    return 1000.0 * noise_raw / (water + 1e-12)


def bg_shell_mask(tmask, soft_mask):
    """Soft-tissue background shell around the tumor (for local noise/mean)."""
    from scipy.ndimage import binary_dilation
    d1 = binary_dilation(tmask, iterations=3)
    d2 = binary_dilation(tmask, iterations=8)
    shell = d2 & (~d1) & soft_mask
    return shell


# ---------------------------------------------------------------------------
def run(sanity=False, reps=1, n_views=180):
    os.makedirs(OUTDIR, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    conrad_backend.setup(max_ram="12G")
    if not conrad_backend.opencl_available():
        print("OpenCL NOT available -- aborting (E5 requires the CL cone pipeline).")
        sys.exit(1)

    kvp, filters, edges, opt = load_optimum()
    E, s = spectrum_at_optimum(kvp, filters)
    global PCD_BINS
    PCD_BINS = edges
    print(f"[E5] spectrum: {int(kvp)} kVp {opt['filter']}  PCD bins {list(edges)}")
    print(f"[E5] geometry: SID={SID_MM} SDD={SDD_MM} voxel={VOXEL_MM}mm "
          f"det_pix={DET_PIX_MM}mm N_det={N_DET} N_vox={N_VOX} views={n_views}")

    n0 = float(opt.get("score_n0", 1.0e5))

    print("[E5] re-materializing RabbitCT reference volume (256^3)...")
    base_vols, soft_mask, refvol = build_base_volumes()
    tmask = tumor_mask(soft_mask)
    n_tvox = int(tmask.sum())
    tumor_vol_cm3 = n_tvox * (VOXEL_MM * 0.1) ** 3
    print(f"[E5] tumor: {n_tvox} voxels ~ {tumor_vol_cm3:.2f} cm^3 at {TUMOR_CENTER_VOX}")
    if n_tvox == 0:
        print("[E5] ERROR: tumor mask empty (no soft tissue at center) -- aborting.")
        sys.exit(1)

    geo = setup_geometry(n_views)

    # ---- base-material cone projections (computed ONCE): path integral [g/cm^2]
    # projector returns path integral in voxel-spacing units; * voxel_mm*0.1 -> g/cm^2
    scale = VOXEL_MM * 0.1
    print("[E5] cone forward projection of base materials (adipose/muscle/bone)...")
    base_paths = {}
    for m, v in base_vols.items():
        base_paths[m] = cone_forward(v, geo) * scale
        print(f"       {m}: path max {base_paths[m].max():.3f} g/cm^2")

    # iron: unit-density iron occupancy in tumor -> path length; scale by c_fe later
    iron_occ = tmask.astype(np.float32)
    iron_path_unit = cone_forward(iron_occ, geo) * scale   # cm of tumor path
    print(f"[E5] iron unit-path (tumor chord) max {iron_path_unit.max():.3f} cm")

    body_diam_cm = 15.0     # rabbit torso ~150 mm for the PCD weight background
    w = pcd_weights(E, s, edges, body_diam_cm)
    print(f"[E5] PCD matched-filter weights: {np.round(w, 4)}")

    shell = bg_shell_mask(tmask, soft_mask)
    print(f"[E5] background shell voxels: {int(shell.sum())}")

    # ---- no-iron reference recon per detector (computed ONCE) -----------------
    # DeltaHU = tumor mean of (iron recon - identical no-iron recon), isolating iron
    # from the anatomy; noise = std of the background ROI on (noisy - noiseless),
    # removing residual anatomy so only quantum noise remains. Both mirror the 2D
    # pipeline's c0-corrected iron_delta_hu and its noise std.
    acc0 = accumulators(base_paths, iron_path_unit * 0.0, E, s, edges, n0)
    rec0 = {}
    for detector in ("EID", "PCD"):
        bh0 = bh_poly_for(acc0, detector, w)
        sino0 = line_integral(acc0, detector, 0, w, bh0, noiseless=True)
        rec0[detector] = fdk_from_sino(sino0, geo)
        wref = float(rec0[detector][shell].mean())
        print(f"[E5] no-iron reference recon [{detector}]: "
              f"bg(water) mean = {wref:.5f}")

    # ---- iron sweep -----------------------------------------------------------
    if sanity:
        cfe_levels = [("sanity_3e8", 3.0e8)]
    else:
        cfe_levels = [("low_3e7", 3.0e7), ("1e8", 1.0e8),
                      ("3e8", 3.0e8), ("1e9", 1.0e9)]
    runs = []
    for dens_label, density in cfe_levels:
        c_fe = PG_FE_SPION_I_FRESH * density * 1e-9      # mg Fe/ml
        c_fe_gcm3 = 1e-3 * c_fe                           # g Fe/cm^3
        iron_path = iron_path_unit * c_fe_gcm3            # g Fe/cm^2
        acc = accumulators(base_paths, iron_path, E, s, edges, n0)

        row = dict(label=dens_label, density=density, c_fe=c_fe)
        recs = {}       # noiseless signal recon per detector (for display + DeltaHU)
        for detector in ("EID", "PCD"):
            bh = bh_poly_for(acc, detector, w)
            # (a) noiseless signal recon -> iron DeltaHU (vs no-iron reference)
            sino_nl = line_integral(acc, detector, 0, w, bh, noiseless=True)
            rec_nl = fdk_from_sino(sino_nl, geo)
            recs[detector] = rec_nl
            msig = measure_iron_signal(rec_nl, rec0[detector], tmask, shell)
            water = msig["water"]
            # (b) noisy recons -> quantum noise (std of bg on noisy - noiseless)
            noises = []
            for seed in range(max(1, reps)):
                sino_n = line_integral(acc, detector, seed, w, bh, noiseless=False)
                rec_n = fdk_from_sino(sino_n, geo)
                nz = rec_n - rec_nl                        # anatomy cancels
                noises.append(1000.0 * float(nz[shell].std()) / (water + 1e-12))
            noise_hu = float(np.mean(noises))
            cnr = abs(msig["delta_hu"]) / (noise_hu + 1e-9)
            row[detector] = dict(delta_hu=msig["delta_hu"], noise_hu=noise_hu,
                                 cnr=cnr, water=water)
            print(f"[E5][{dens_label} c_Fe={c_fe:.2f} mg/ml][{detector}] "
                  f"DeltaHU={row[detector]['delta_hu']:+.1f}  "
                  f"noiseHU={noise_hu:.2f}  CNR={cnr:.2f}")

        # ---- save tumor-slice comparison PNG (EID vs PCD) --------------------
        # Left two: the anatomy recon (EID / PCD) at the tumor slice, in a HU-like
        # window, with the tumor outline. Right: the iron-only signal map
        # (recon - no-iron reference) in HU, where the +DeltaHU tumor is visible.
        zc = TUMOR_CENTER_VOX[0]
        from scipy.ndimage import binary_erosion
        tsl = tmask[zc]
        outline = tsl & ~binary_erosion(tsl)
        fig, ax = plt.subplots(1, 3, figsize=(14, 4.6))
        for a, (det, title) in zip(ax[:2],
                                   [("EID", "EID (energy-weighted)"),
                                    ("PCD", "PCD (post-comb matched filter)")]):
            sl = recs[det][zc]
            wr = row[det]["water"]
            a.imshow(sl, cmap="gray", vmin=wr * 0.55, vmax=wr * 1.7)
            a.imshow(np.ma.masked_where(~outline, outline), cmap="autumn",
                     vmin=0, vmax=1)
            a.set_title(f"{title}\nDeltaHU={row[det]['delta_hu']:+.0f} "
                        f"CNR={row[det]['cnr']:.1f}", fontsize=10)
            a.axis("off")
        # iron-only signal map for PCD (HU): (recon - no-iron ref)/water*1000
        iron_hu = 1000.0 * (recs["PCD"][zc] - rec0["PCD"][zc]) / (row["PCD"]["water"] + 1e-9)
        mx = max(5.0, float(np.percentile(np.abs(iron_hu), 99.9)))
        im = ax[2].imshow(iron_hu, cmap="magma", vmin=0, vmax=mx)
        ax[2].imshow(np.ma.masked_where(~outline, outline), cmap="winter",
                     vmin=0, vmax=1)
        fig.colorbar(im, ax=ax[2], fraction=0.046, label="iron signal [HU]")
        ax[2].set_title("PCD iron-only signal (recon - no-iron)", fontsize=10)
        ax[2].axis("off")
        fig.suptitle(f"E5 3D cone-beam FDK -- {int(kvp)}kVp {opt['filter']} -- "
                     f"SPION I fresh {dens_label} (c_Fe={c_fe:.2f} mg Fe/ml, "
                     f"tumor {tumor_vol_cm3:.1f} cm^3)", fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        p = os.path.join(OUTDIR, f"e5_recon_{dens_label}.png")
        fig.savefig(p, dpi=120); plt.close(fig)
        print(f"[E5] wrote {p}")
        row["png"] = p
        runs.append(row)

    # ---- summary + JSON -------------------------------------------------------
    import json
    print("\n================ E5 SUMMARY (3D cone-beam FDK) ================")
    print(f"{'run':>10} {'c_Fe':>7} | {'EID dHU':>8} {'EID CNR':>8} | "
          f"{'PCD dHU':>8} {'PCD CNR':>8} | {'CNR PCD/EID':>11}")
    for r in runs:
        gain = r["PCD"]["cnr"] / (r["EID"]["cnr"] + 1e-9)
        print(f"{r['label']:>10} {r['c_fe']:>7.2f} | "
              f"{r['EID']['delta_hu']:>8.1f} {r['EID']['cnr']:>8.2f} | "
              f"{r['PCD']['delta_hu']:>8.1f} {r['PCD']['cnr']:>8.2f} | {gain:>11.2f}")

    meta = dict(
        study="E5_3d_cone_fdk", spectrum=dict(kvp=kvp, filter=opt["filter"],
        pcd_bin_edges_kev=list(edges)),
        geometry=dict(sid_mm=SID_MM, sdd_mm=SDD_MM, voxel_mm=VOXEL_MM,
                      det_pix_mm=DET_PIX_MM, n_det=N_DET, n_vox=N_VOX, n_views=n_views),
        tumor=dict(center_vox=TUMOR_CENTER_VOX, radius_mm=TUMOR_RADIUS_MM,
                   n_voxels=n_tvox, volume_cm3=tumor_vol_cm3,
                   spion="SPION_I_fresh", pg_fe_per_cell=PG_FE_SPION_I_FRESH),
        n0=n0, reps=reps, pcd_weights=[float(x) for x in w],
        runs=[{k: (v if not isinstance(v, dict) else v)
               for k, v in r.items() if k != "png"} for r in runs])
    jpath = os.path.join(OUTDIR, "e5_results.json")
    with open(jpath, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"\n[E5] wrote {jpath}")
    return runs


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--sanity", action="store_true",
                    help="single noiseless sanity run at 3e8 cells")
    ap.add_argument("--reps", type=int, default=1, help="noise realizations per cell")
    ap.add_argument("--views", type=int, default=180, help="projection views over 360")
    args = ap.parse_args()
    run(sanity=args.sanity, reps=args.reps, n_views=args.views)
