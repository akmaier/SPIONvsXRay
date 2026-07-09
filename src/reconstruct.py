"""M5 — parallel-beam filtered backprojection + ROI measurement.

Reconstructs the EID and PCD line-integral sinograms from simulate.py into
attenuation images, with optional water beam-hardening correction (M6 factor).
HU and CNR are computed from tumor vs. reference ROIs; absolute FBP scaling
cancels in HU (defined relative to the reconstructed background), so no
calibration is needed.
"""
from __future__ import annotations
import numpy as np
from scipy.ndimage import rotate

import conrad_backend
from config import PHANTOM


def ramp_filter(sino):
    """Apply a ramp (Ram-Lak) filter with a Hann window along the detector axis."""
    n_ang, n_det = sino.shape
    n = 1
    while n < 2 * n_det:
        n *= 2
    freqs = np.fft.fftfreq(n)
    ramp = 2.0 * np.abs(freqs)
    ramp *= 0.5 * (1 + np.cos(2 * np.pi * freqs))   # Hann window
    out = np.empty_like(sino)
    for i in range(n_ang):
        f = np.fft.fft(sino[i], n=n)
        out[i] = np.real(np.fft.ifft(f * ramp))[:n_det]
    return out


def fbp(sino, angles_deg, N=None):
    """Parallel-beam filtered backprojection -> attenuation image (a.u.)."""
    n_ang, n_det = sino.shape
    if N is None:
        N = n_det
    filt = ramp_filter(sino)
    recon = np.zeros((N, N))
    for i, a in enumerate(angles_deg):
        # smear the filtered projection into an image, then rotate to the view
        strip = np.tile(filt[i], (N, 1))
        recon += rotate(strip, -a, reshape=False, order=1, mode="constant", cval=0.0)
    recon *= np.pi / (2.0 * n_ang)
    return recon


def water_precorrection(p, coeffs=(1.0, 0.10)):
    """Simple polynomial beam-hardening (water) precorrection: p -> p + c2*p^2."""
    c1, c2 = coeffs
    return c1 * p + c2 * p * p


def reconstruct(sino, angles, bh_correction=False, N=None):
    s = water_precorrection(sino) if bh_correction else sino
    return fbp(s, angles, N=N)


def _mask(N, voxel_cm, center_mm, radius_mm):
    # build_components: X -> axis 0 (rows), Y -> axis 1 (cols); mgrid yy=rows, xx=cols
    yy, xx = np.mgrid[0:N, 0:N]
    row_c = center_mm[0] / (voxel_cm * 10.0) + N / 2.0   # X -> row
    col_c = center_mm[1] / (voxel_cm * 10.0) + N / 2.0   # Y -> col
    r = radius_mm / (voxel_cm * 10.0)
    return (xx - col_c) ** 2 + (yy - row_c) ** 2 <= r ** 2


def measure(img, voxel_cm, tumor_radius_mm=8.0, ref_radius_mm=8.0):
    """ΔHU and CNR of the tumor vs a soft-tissue reference ROI.

    Reference is placed at the SAME radius from iso-center as the tumor so
    radial beam-hardening cupping cancels in ΔHU.
    """
    N = img.shape[0]
    tumor_c = (PHANTOM.tumor_center_mm[0], PHANTOM.tumor_center_mm[1])   # (25, 0)
    ref_c = (0.0, -25.0)   # same radius (25 mm), soft tissue, clear of bone/tumor
    tm = _mask(N, voxel_cm, tumor_c, tumor_radius_mm)
    rm = _mask(N, voxel_cm, ref_c, ref_radius_mm)
    mu_t, mu_r = img[tm].mean(), img[rm].mean()
    sd_r = img[rm].std()
    d_hu = 1000.0 * (mu_t - mu_r) / mu_r
    cnr = (mu_t - mu_r) / (sd_r + 1e-12)
    hu_noise = 1000.0 * sd_r / mu_r
    return dict(delta_hu=float(d_hu), cnr=float(cnr), hu_noise=float(hu_noise),
                mu_tumor=float(mu_t), mu_ref=float(mu_r))


def figures(outdir=None):
    import os, matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import simulate
    if outdir is None:
        outdir = str(conrad_backend.REPO_ROOT / "results" / "recon")
    os.makedirs(outdir, exist_ok=True)

    N = 512
    voxel_cm = 20.0 / N
    proj = simulate.project_materials(20.0, "homogeneous", n_views=500, N=N)
    det = simulate.detector_sinograms(proj, add_noise=True, seed=1)
    img_eid = reconstruct(det["eid"], det["angles"], bh_correction=False, N=N)
    meas = measure(img_eid, voxel_cm)
    print(f"[check] EID recon (c_form=20): dHU={meas['delta_hu']:.2f}, "
          f"CNR={meas['cnr']:.2f}, noise={meas['hu_noise']:.1f} HU")

    # zoom on tumor
    cx = int(PHANTOM.tumor_center_mm[0] / (voxel_cm * 10) + N / 2)
    cy = int(PHANTOM.tumor_center_mm[1] / (voxel_cm * 10) + N / 2)
    w = 90
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.4))
    axes[0].imshow(img_eid, cmap="gray"); axes[0].set_title("EID FBP reconstruction")
    axes[0].axis("off")
    axes[1].imshow(img_eid[cy - w:cy + w, cx - w:cx + w], cmap="gray")
    axes[1].set_title(f"tumor zoom  dHU={meas['delta_hu']:.1f}, CNR={meas['cnr']:.2f}")
    axes[1].axis("off")
    fig.tight_layout(); fig.savefig(f"{outdir}/recon_tumor_zoom.png", dpi=130); plt.close(fig)
    print("[ok] wrote", outdir, "-> recon_tumor_zoom.png")
    return img_eid, meas


if __name__ == "__main__":
    figures()
