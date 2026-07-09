"""CONRAD fan-beam projection + FBP reconstruction (GPU projector + CPU FBP).

OpenCL (Path A) is enabled on this Apple-Silicon machine via jogamp 2.6.0 +
an OpenCL.framework reexport shim (see scripts/install_opencl.sh and
conrad_backend). When available, the FORWARD projector uses CONRAD's
`projectRayDrivenCL` (validated: 0.03% vs CPU, ~4000x faster on the M1 GPU).
Backprojection stays on the CPU `backprojectPixelDriven` for now (the CL
backprojector has an unresolved convention mismatch — see fbp()). Everything
falls back to CPU if OpenCL is absent.

Uses edu.stanford.rsl.tutorial.fan.FanBeamProjector2D and
FanBeamBackprojector2D with RamLakRampFilter; numpy<->Grid2D via pyconrad.
"""
from __future__ import annotations
import numpy as np
import jpype

import conrad_backend

# --- geometry (fan-beam 2D analogue of the C-arm cone beam) ---
SDD_MM = 1200.0        # source-detector distance
SID_MM = 750.0         # source-isocenter distance = CONRAD focalLength
N_VIEWS = 500
FOV_MM = 200.0
# Virtual-detector bin size at the isocenter [mm]. CONRAD's tutorial.fan chain
# (FanBeamProjector/CosineFilter/RamLak-SheppLogan/Backprojector) is a self-consistent
# unit system calibrated for 1 grid-unit = 1 mm; the ramp kernels' deltaS scaling is
# only exact at deltaT=1 (SheppLoganKernel carries a 1/deltaS^2, RamLakKernel a
# 1/deltaS). Keeping deltaT=1 mm -- decoupled from the recon voxel -- makes the whole
# chain exact (a uniform water cylinder reconstructs flat to 0.1%). See DEVLOG.
DETECTOR_DT_MM = 1.0


def _cls(pkg, name):
    conrad_backend.setup()   # idempotent: ensures the JVM is up
    return conrad_backend.class_getter(pkg).__getattr__(name)


def np_to_grid2d(a: np.ndarray):
    """numpy (H,W) float -> CONRAD Grid2D (fast float[] constructor, row-major)."""
    a = np.ascontiguousarray(a, dtype=np.float32)
    G = _cls("edu.stanford.rsl.conrad.data.numeric", "Grid2D")
    h, w = a.shape
    jbuf = jpype.JArray(jpype.JFloat)(a.ravel(order="C"))
    return G(jbuf, int(w), int(h))


def grid2d_to_np(g) -> np.ndarray:
    w, h = int(g.getWidth()), int(g.getHeight())
    buf = np.array(g.getBuffer()[:], dtype=np.float32)
    return buf.reshape(h, w)


def fan_geometry(fov_mm=FOV_MM, n_pix=None, voxel_mm=None, short_scan=False):
    """Fan-beam geometry in CONRAD's convention (matches FanBeamReconstructionExample).

    CONRAD's tutorial.fan uses a *virtual detector at the isocenter*: the source
    sits at radius `focal` = source->isocenter distance (SID), and the detector
    line passes through the isocenter spanning +-maxT/2 [mm], sampled at deltaT
    [mm/bin] (n_det = maxT/deltaT). There is NO magnification or padding factor --
    the fan half-angle follows from the geometry alone,
        gamma = atan((maxT/2 - deltaT/2) / focal).

    Recon grid = n_pix x n_pix at `voxel_mm` (config.RECON_VOXEL_MM, 0.5 mm). The
    virtual detector is sized to the recon FOV (maxT = n_pix*voxel_mm) and sampled at
    deltaT = DETECTOR_DT_MM = 1 mm -- kept at 1 mm (NOT voxel) so CONRAD's ramp kernels
    stay exact (their deltaS scaling is only self-consistent at deltaT=1). The recon
    still uses `voxel_mm` spacing via the backprojector; a centered object up to the
    FOV reconstructs without truncation and a uniform water cylinder stays flat to 0.1%.
    short_scan=True uses the MINIMAL short scan (180 deg + 2*gamma) with Parker
    redundancy weighting in fbp (>= minimal is required for Parker to normalize
    correctly -- see the CONRAD weighting-comparison study); default is 360 deg.
    """
    from config import RECON_VOXEL_MM
    if n_pix is None:
        n_pix = 512
    if voxel_mm is None:
        voxel_mm = RECON_VOXEL_MM
    focal = SID_MM                                   # CONRAD focalLength = source->isocenter
    deltaT = DETECTOR_DT_MM                           # virtual-detector bin size at isocenter (1 mm)
    maxT = n_pix * voxel_mm                          # detector length = recon FOV [mm]
    gamma_max = float(np.arctan((maxT / 2.0 - 0.5 * deltaT) / focal))   # half fan angle
    maxBeta = (np.pi + 2.0 * gamma_max) if short_scan else 2.0 * np.pi
    deltaBeta = maxBeta / N_VIEWS
    return dict(focal=focal, maxBeta=maxBeta, deltaBeta=deltaBeta,
               maxT=maxT, deltaT=deltaT, imgN=n_pix, spacing=voxel_mm,
               voxel_mm=voxel_mm, gamma_max=gamma_max, short_scan=short_scan)


def use_cl():
    """Whether to use CONRAD's OpenCL projectors (Path A GPU) on this machine."""
    return conrad_backend.opencl_available()


def project(image_np, geo, gpu=None):
    g = np_to_grid2d(image_np)
    g.setSpacing(geo["spacing"], geo["spacing"])
    g.setOrigin(-(geo["imgN"] - 1) / 2.0 * geo["spacing"],
                -(geo["imgN"] - 1) / 2.0 * geo["spacing"])
    FP = _cls("edu.stanford.rsl.tutorial.fan", "FanBeamProjector2D")
    fp = FP(geo["focal"], geo["maxBeta"], geo["deltaBeta"], geo["maxT"], geo["deltaT"])
    gpu = use_cl() if gpu is None else gpu
    if gpu:
        try:
            return grid2d_to_np(fp.projectRayDrivenCL(g))     # GPU
        except Exception:
            pass                                              # fall back to CPU
    return grid2d_to_np(fp.projectRayDriven(g))               # CPU


def water_precorrection_poly(E, w_eff, mu_water_percm, Lmax_cm=30.0, deg=4):
    """Calibrated water beam-hardening precorrection polynomial.

    Given a detector's effective per-energy weighting w_eff(E) (EID: s*E;
    PCD-combined: bin-weight*s) and water linear attenuation mu(E) [1/cm], tabulate
    the polychromatic water line integral p_poly(L) = -log(sum w exp(-mu L)) and fit
    p_mono = poly(p_poly), where p_mono = mu_ref*L is the linear (monochromatic)
    response at the spectrum-weighted mean mu_ref. Applying poly() to a measured
    line integral linearizes the water response, so a water cylinder reconstructs
    flat (removing cupping) instead of relying on a guessed quadratic coefficient.
    """
    w = np.asarray(w_eff, float); w = w / w.sum()
    L = np.linspace(0.0, Lmax_cm, 400)
    p_poly = np.array([-np.log(np.sum(w * np.exp(-mu_water_percm * l))) for l in L])
    mu_ref = float(np.sum(w * mu_water_percm))        # low-dose slope dp/dL at L=0
    return np.polyfit(p_poly, mu_ref * L, deg)


def fbp(sino_np, geo, bh_correction=False, bh_poly=None, bh_c2=0.10):
    """Fan-beam FBP using the CONRAD API end-to-end (project-side; see README policy).

    Pipeline (identical order to CONRAD's FanBeamReconstructionExample):
      (optional water beam-hardening precorrection) -> ParkerWeights (short scan)
      -> CosineFilter -> SheppLoganKernel (ramp with roll-off) -> distance-weighted
      FanBeamBackprojector2D. Redundancy weighting, cosine, ramp and backprojection
      are all CONRAD classes; only the spectral water precorrection is numpy.

    Water precorrection: pass a calibrated `bh_poly` (from water_precorrection_poly);
    the legacy fixed `bh_c2` quadratic is only a fallback and is NOT spectrum-matched.
    """
    if bh_correction:
        if bh_poly is not None:
            sino_np = np.polyval(bh_poly, sino_np)    # calibrated water precorrection
        else:
            sino_np = sino_np + bh_c2 * sino_np ** 2   # legacy fallback (uncalibrated)
    focal, maxT, deltaT = float(geo["focal"]), float(geo["maxT"]), float(geo["deltaT"])
    maxBeta, deltaBeta = float(geo["maxBeta"]), float(geo["deltaBeta"])

    sino = np_to_grid2d(sino_np)                      # Grid2D(width=n_det, height=n_views)
    sino.setSpacing(deltaT, deltaBeta)
    NPO = _cls("edu.stanford.rsl.conrad.data.numeric", "NumericPointwiseOperators")

    # short-scan redundancy weighting (CONRAD ParkerWeights) BEFORE filtering. Do
    # NOT re-derive Parker from the 1982 paper -- its equations contain a known typo.
    if geo.get("short_scan"):
        PW = _cls("edu.stanford.rsl.tutorial.fan.redundancy", "ParkerWeights")
        NPO.multiplyBy(sino, PW(focal, maxT, deltaT, maxBeta, deltaBeta))

    # fan cosine weight + ramp with roll-off (Shepp-Logan), applied per view
    Cos = _cls("edu.stanford.rsl.tutorial.fan", "CosineFilter")
    SL = _cls("edu.stanford.rsl.tutorial.filters", "SheppLoganKernel")
    cos = Cos(focal, maxT, deltaT)
    ram = SL(int(round(maxT / deltaT)), deltaT)
    n_views = int(sino.getSize()[1])
    for th in range(n_views):
        cos.applyToGrid(sino.getSubGrid(th))
    for th in range(n_views):
        ram.applyToGrid(sino.getSubGrid(th))

    # Distance-weighted backprojection. Patched FanBeamBackprojector2D (conrad_ext,
    # shadows the jar): the 1/U^2 fan distance weighting is reinstated (CPU + CL),
    # and setSpacing(vox) feeds the CL kernel the recon pixel spacing that upstream
    # omitted -> configurable voxel size and backprojectPixelDrivenCL ~850x faster.
    vox = geo.get("voxel_mm", 1.0)
    BP = _cls("edu.stanford.rsl.tutorial.fan", "FanBeamBackprojector2D")
    bp = BP(focal, deltaT, deltaBeta, int(geo["imgN"]), int(geo["imgN"]))
    if use_cl():
        try:
            bp.setSpacing(float(vox))
            return grid2d_to_np(bp.backprojectPixelDrivenCL(sino))   # GPU, spacing-aware
        except Exception:
            pass
    return grid2d_to_np(bp.backprojectPixelDriven(sino))             # CPU fallback


if __name__ == "__main__":
    import os, matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    conrad_backend.setup()
    N = 256
    geo = fan_geometry(n_pix=N)
    # disk phantom
    yy, xx = np.mgrid[0:N, 0:N]
    disk = (((xx - N/2)**2 + (yy - N/2)**2) < (N*0.28)**2).astype(np.float32)
    disk += (((xx - N*0.62)**2 + (yy - N/2)**2) < (N*0.06)**2).astype(np.float32)
    sino = project(disk, geo)
    recon = fbp(sino, geo)
    print("sino shape", sino.shape, "range", float(sino.min()), float(sino.max()))
    print("recon shape", recon.shape, "range", float(recon.min()), float(recon.max()))
    fig, ax = plt.subplots(1, 3, figsize=(12, 4))
    ax[0].imshow(disk, cmap="gray"); ax[0].set_title("phantom (disk)")
    ax[1].imshow(sino, aspect="auto", cmap="gray"); ax[1].set_title("CONRAD fan sinogram")
    ax[2].imshow(recon, cmap="gray"); ax[2].set_title("CONRAD fan-beam FBP recon")
    for a in ax: a.axis("off")
    os.makedirs("results", exist_ok=True)
    fig.tight_layout(); fig.savefig("results/conrad_roundtrip.png", dpi=120)
    print("wrote results/conrad_roundtrip.png")
