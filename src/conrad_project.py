"""M4 (CONRAD-native) — fan-beam base-material sinograms via projectRayDrivenCL.

Each material is rasterized into an indicator grid on the detector-matched grid
(1 px = deltaT mm, so FanBeamProjector2D geometry is in mm) and forward-projected
with CONRAD's OpenCL grid projector FanBeamProjector2D.projectRayDrivenCL, giving
base-material path lengths [mm]. The projector geometry matches conrad_ct.fbp
(FanBeamBackprojector2D), so the sinograms reconstruct directly. Path lengths are
combined polychromatically (per-material CONRAD attenuation over the real 90 kVp
spectrum) into EID + multi-bin PCD line-integral sinograms with Poisson noise.
Sinogram layout: [views, det], matching conrad_ct.
"""
from __future__ import annotations
import numpy as np

import conrad_backend
import conrad_ct
import conrad_phantom
import spectrum as spec
from config import SPECTRUM, DETECTORS

CG = conrad_backend.class_getter
N0 = SPECTRUM.photons_per_pixel


def project_base_materials(inserts, geo=None, n_pix=512):
    """CONRAD-native fan-beam base-material sinograms (OpenCL projector).

    Rasterizes each material into an indicator grid on the detector-matched grid
    (1 px = deltaT mm, so FanBeamProjector2D geometry is in mm) and projects each
    with projectRayDrivenCL. Returns ({material_name: sino[views, det]} in mm, geo),
    matching the recon geometry in `geo`. Inserts override the water body, matching
    the PriorityRayTracer scene ordering (inserts added after the body).
    """
    conrad_backend.setup()
    if geo is None:
        geo = conrad_ct.fan_geometry(n_pix=n_pix)
    n = int(round(geo["maxT"] / geo["deltaT"]))          # detector-matched, 1 px = deltaT mm
    vox = geo["deltaT"]
    yy, xx = np.mgrid[0:n, 0:n]
    x = (xx - n / 2.0) * vox
    y = (yy - n / 2.0) * vox
    assigned = np.zeros((n, n), bool)
    masks = {}
    for ins in inserts:                                  # inserts override body
        cx, cy = ins["center_mm"]; r = ins["radius_mm"]; nm = ins["name"]
        m = ((x - cx) ** 2 + (y - cy) ** 2 <= r * r) & ~assigned
        masks.setdefault(nm, np.zeros((n, n), np.float32))[m] = 1.0
        assigned |= m
    body = ((x * x + y * y) <= conrad_phantom.BODY_RADIUS_MM ** 2) & ~assigned
    masks.setdefault("water", np.zeros((n, n), np.float32))[body] = 1.0
    FP = CG("edu.stanford.rsl.tutorial.fan").FanBeamProjector2D
    fp = FP(geo["focal"], geo["maxBeta"], geo["deltaBeta"], geo["maxT"], geo["deltaT"])
    base = {}
    for nm, mask in masks.items():
        g = conrad_ct.np_to_grid2d(mask)
        g.setSpacing(vox, vox)
        g.setOrigin(-(n * vox) / 2.0, -(n * vox) / 2.0)
        base[nm] = conrad_ct.grid2d_to_np(fp.projectRayDrivenCL(g)) * vox   # path length [mm]
    return base, geo


def _material_mu(names, energies):
    """Per-material linear attenuation mu(E) [1/cm] from CONRAD, for given names."""
    DB = CG("edu.stanford.rsl.conrad.physics.materials.database").MaterialsDB
    AT = CG("edu.stanford.rsl.conrad.physics.materials.utils").AttenuationType
    at = AT.TOTAL_WITH_COHERENT_ATTENUATION
    mu = {}
    for name in names:
        mat = DB.getMaterialWithName(name)
        mu[name] = np.array([float(mat.getAttenuation(float(e), at)) for e in energies])
    return mu


def detector_sinograms(base_sinos, kvp=None, filters=(), add_noise=True, seed=0):
    """Combine base-material sinograms polychromatically -> EID + PCD sinograms."""
    rng = np.random.default_rng(seed)
    if kvp is None:
        E, flux, _ = spec.standard_spectrum()
    else:
        E, flux = spec.conrad_spectrum(kvp)
    if filters:
        flux = spec.apply_filters(E, flux, list(filters))
    s = flux / flux.sum()

    names = list(base_sinos.keys())
    mu = _material_mu(names, E)                         # [1/cm]
    shape = next(iter(base_sinos.values())).shape

    S_det = np.zeros(shape); S_det_E2 = np.zeros(shape)
    S_air = float(np.sum(N0 * s * E))
    edges = np.array(DETECTORS.pcd_bin_edges_kev)
    nb = len(edges) - 1
    C_det = [np.zeros(shape) for _ in range(nb)]
    C_air = [0.0] * nb

    # path lengths mm -> cm
    L = {name: base_sinos[name] * 0.1 for name in names}
    for j, Ej in enumerate(E):
        tau = np.zeros(shape)
        for name in names:
            tau += L[name] * mu[name][j]
        n = N0 * s[j] * np.exp(-tau)
        S_det += n * Ej
        S_det_E2 += n * Ej * Ej
        b = int(np.searchsorted(edges, Ej, side="right") - 1)
        if 0 <= b < nb:
            C_det[b] += n
            C_air[b] += N0 * s[j]

    if add_noise:
        S_meas = S_det + rng.normal(0.0, np.sqrt(np.maximum(S_det_E2, 1e-30)))
        C_meas = [rng.poisson(np.maximum(c, 0.0)) for c in C_det]
    else:
        S_meas, C_meas = S_det, C_det
    eps = 1e-6
    p_eid = -np.log(np.clip(S_meas, eps, None) / S_air)
    p_pcd = [-np.log(np.clip(cm, eps, None) / max(ca, eps)) for cm, ca in zip(C_meas, C_air)]
    return dict(eid=p_eid, pcd=p_pcd, edges=edges)


def measure_inserts(recon, geo, inserts, roi_mm=8.0, bg_inner_mm=15.0, bg_outer_mm=22.0):
    """Per-insert ΔHU/noise vs a LOCAL annular background around each insert.

    A ring just outside each insert (all within the water body) shares that
    insert's local cupping/streak level, so subtracting it cancels the
    bone-streak bias that a single global reference cannot. ΔHU = 1000·(μ_insert −
    μ_localBG)/μ_localBG; noise = σ(localBG) in HU (for CNR).
    """
    N = recon.shape[0]
    sp = geo.get("voxel_mm", 1.0)   # recon grid spacing (mm/px), set by fbp
    yy, xx = np.mgrid[0:N, 0:N]
    out = []
    for ins in inserts:
        cx, cy = ins["center_mm"]
        col, row = cx / sp + N / 2.0, cy / sp + N / 2.0
        r2 = (xx - col) ** 2 + (yy - row) ** 2
        roi = r2 <= (roi_mm / sp) ** 2
        bg = (r2 >= (bg_inner_mm / sp) ** 2) & (r2 <= (bg_outer_mm / sp) ** 2)
        mu_i = float(recon[roi].mean())
        mu_b = float(recon[bg].mean())
        sd_b = float(recon[bg].std())
        out.append({**ins, "mu": mu_i, "mu_bg": mu_b,
                    "delta_hu": 1000.0 * (mu_i - mu_b) / mu_b,
                    "hu_noise": 1000.0 * sd_b / mu_b})
    # c0 correction: subtract the zero-iron (SPION_c0) insert's ΔHU, which is a
    # geometric insert-vs-annulus bias common to all inserts -> pure iron ΔHU.
    c0 = next((m["delta_hu"] for m in out if m["name"] == "SPION_c0"), 0.0)
    for m in out:
        m["iron_delta_hu"] = m["delta_hu"] - c0
    return out


if __name__ == "__main__":
    import time, os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    scene, inserts = conrad_phantom.build_phantom()
    geo = conrad_ct.fan_geometry(n_pix=512)
    t = time.time()
    base, geo = project_base_materials(inserts, geo)
    print(f"[native] base-material projection {time.time()-t:.1f}s, "
          f"sino {next(iter(base.values())).shape}, materials {len(base)}")

    det = detector_sinograms(base, add_noise=False)
    recon = conrad_ct.fbp(det["eid"], geo)
    meas = measure_inserts(recon, geo, inserts)
    print("[native] per-insert IRON ΔHU (local-annulus, c0-corrected; noise-free EID):")
    for m in meas:
        if m["c_form"] is not None:
            print(f"  {m['name']:10s} c_Fe={m['c_fe']:.3f}  iron_dHU={m['iron_delta_hu']:+.2f}")

    outdir = str(conrad_backend.REPO_ROOT / "results" / "native")
    os.makedirs(outdir, exist_ok=True)
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.6))
    ax[0].imshow(base["water"], aspect="auto", cmap="viridis")
    ax[0].set_title("base-material sinogram: water [mm]"); ax[0].set_xlabel("detector"); ax[0].set_ylabel("view")
    ax[1].imshow(det["eid"], aspect="auto", cmap="gray")
    ax[1].set_title("EID line-integral sinogram")
    lo, hi = np.percentile(recon, [30, 99.5])
    ax[2].imshow(recon, cmap="gray", vmin=lo, vmax=hi)
    ax[2].set_title("CONRAD fan-beam FBP reconstruction")
    for a in (ax[2],):
        a.axis("off")
    fig.tight_layout(); fig.savefig(f"{outdir}/native_pipeline.png", dpi=130); plt.close(fig)
    print("[ok] wrote", outdir, "-> native_pipeline.png")
