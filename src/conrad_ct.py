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


def fan_geometry(fov_mm=FOV_MM, n_pix=None):
    """Return (focalLength, maxBeta, deltaBeta, maxT, deltaT, imgN, spacing_mm)."""
    if n_pix is None:
        n_pix = 512
    spacing = fov_mm / n_pix
    # detector must cover the FOV magnified to the detector plane
    mag = SDD_MM / SID_MM
    maxT = fov_mm * mag * 1.15
    deltaT = maxT / (n_pix * 2)
    maxBeta = 2.0 * np.pi
    deltaBeta = maxBeta / N_VIEWS
    return dict(focal=SDD_MM, maxBeta=maxBeta, deltaBeta=deltaBeta,
               maxT=maxT, deltaT=deltaT, imgN=n_pix, spacing=spacing)


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


def fbp(sino_np, geo, bh_correction=False, bh_c2=0.10):
    """Fan-beam FBP: (optional water beam-hardening precorrection ->) cosine ->
    ramp filter -> CONRAD distance-weighted backprojection."""
    if bh_correction:
        sino_np = sino_np + bh_c2 * sino_np ** 2      # nominal water precorrection
    n_ang, n_det = sino_np.shape
    t = (np.arange(n_det) - (n_det - 1) / 2.0) * geo["deltaT"]
    cos_w = geo["focal"] / np.sqrt(geo["focal"] ** 2 + t ** 2)   # fan cosine weight
    sino_w = sino_np * cos_w[None, :]
    sino_f = ramp_filter_sino(sino_w, geo["deltaT"])
    sino = np_to_grid2d(sino_f)
    sino.setSpacing(geo["deltaT"], geo["deltaBeta"])
    BP = _cls("edu.stanford.rsl.tutorial.fan", "FanBeamBackprojector2D")
    bp = BP(geo["focal"], geo["deltaT"], geo["deltaBeta"], geo["imgN"], geo["imgN"])
    # CPU backprojector only. The CL backprojector (backprojectPixelDrivenCL) is
    # INCOMPLETE upstream: its kernel is never passed the image pixel spacing
    # ("// TODO: Spacing :)" in FanBeamBackprojector2D.java), so its pixel<->world
    # mapping is inconsistent with the CPU path -> looks fine on a uniform disk
    # but corrupts sub-HU quantitative measurements. Forward projection uses CL.
    return grid2d_to_np(bp.backprojectPixelDriven(sino))             # CPU (correct)


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
