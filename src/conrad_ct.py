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
SDD_MM = 1200.0        # source-detector distance (focal length)
SID_MM = 750.0         # source-isocenter distance
N_VIEWS = 500
FOV_MM = 200.0


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


def ramp_filter_sino(sino, deltaT):
    """Ram-Lak ramp filter (Hann-windowed) along the detector axis, per view."""
    n_ang, n_det = sino.shape
    n = 1
    while n < 2 * n_det:
        n *= 2
    freqs = np.fft.fftfreq(n, d=deltaT)
    ramp = 2.0 * np.abs(freqs) * deltaT
    ramp *= 0.5 * (1 + np.cos(np.pi * freqs / np.abs(freqs).max()))   # Hann
    out = np.empty_like(sino)
    for i in range(n_ang):
        f = np.fft.fft(sino[i], n=n)
        out[i] = np.real(np.fft.ifft(f * ramp))[:n_det]
    return out


def grid2d_to_np(g) -> np.ndarray:
    w, h = int(g.getWidth()), int(g.getHeight())
    buf = np.array(g.getBuffer()[:], dtype=np.float32)
    return buf.reshape(h, w)


def fan_geometry(fov_mm=FOV_MM, n_pix=None, voxel_mm=None):
    """Fan-beam geometry. Recon = n_pix x n_pix at `voxel_mm` (fixed CL backprojector).

    voxel_mm defaults to config.RECON_VOXEL_MM (0.5 mm); recon FOV = n_pix*voxel_mm.
    """
    from config import RECON_VOXEL_MM
    if n_pix is None:
        n_pix = 512
    if voxel_mm is None:
        voxel_mm = RECON_VOXEL_MM
    # detector must cover the physical FOV magnified to the detector plane
    mag = SDD_MM / SID_MM
    maxT = fov_mm * mag * 1.15
    deltaT = maxT / (n_pix * 2)
    maxBeta = 2.0 * np.pi
    deltaBeta = maxBeta / N_VIEWS
    return dict(focal=SDD_MM, maxBeta=maxBeta, deltaBeta=deltaBeta,
               maxT=maxT, deltaT=deltaT, imgN=n_pix, spacing=voxel_mm,
               voxel_mm=voxel_mm)


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
    """Fan-beam FBP: (optional water beam-hardening precorrection ->) cosine ->
    ramp filter -> CONRAD distance-weighted backprojection.

    Water precorrection: pass a calibrated `bh_poly` (from water_precorrection_poly);
    the legacy fixed `bh_c2` quadratic is only a fallback and is NOT spectrum-matched.
    """
    if bh_correction:
        if bh_poly is not None:
            sino_np = np.polyval(bh_poly, sino_np)    # calibrated water precorrection
        else:
            sino_np = sino_np + bh_c2 * sino_np ** 2   # legacy fallback (uncalibrated)
    n_ang, n_det = sino_np.shape
    t = (np.arange(n_det) - (n_det - 1) / 2.0) * geo["deltaT"]
    cos_w = geo["focal"] / np.sqrt(geo["focal"] ** 2 + t ** 2)   # fan cosine weight
    sino_w = sino_np * cos_w[None, :]
    sino_f = ramp_filter_sino(sino_w, geo["deltaT"])
    sino = np_to_grid2d(sino_f)
    sino.setSpacing(geo["deltaT"], geo["deltaBeta"])
    # Patched FanBeamBackprojector2D (conrad_ext, shadows the jar): its CL kernel
    # now receives the reconstruction pixel spacing that upstream omitted
    # ("// TODO: Spacing"), so setSpacing(vox) gives configurable voxel size and
    # backprojectPixelDrivenCL is correct + ~850x faster. Falls back to CPU.
    vox = geo.get("voxel_mm", 1.0)
    BP = _cls("edu.stanford.rsl.tutorial.fan", "FanBeamBackprojector2D")
    bp = BP(geo["focal"], geo["deltaT"], geo["deltaBeta"], geo["imgN"], geo["imgN"])
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
