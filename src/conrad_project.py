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


def _disk_chords(inserts, geo):
    """Exact fan-beam projection of the all-cylinder phantom (closed-form chords).

    Every object is a centered/offset cylinder, so the line integral through a disk
    is the analytic chord 2*sqrt(r^2 - d^2), d = perpendicular ray-to-center
    distance. Inserts sit inside the water body and override it, so
    water_path = body_chord - sum(insert_chords). Exact (no rasterisation bias);
    the ground-truth forward used to validate the anti-aliased grid projection.
    Path lengths [mm], sinogram [views, det] matching FanBeamProjector2D.
    """
    focal, maxT, deltaT = geo["focal"], geo["maxT"], geo["deltaT"]
    maxBeta, deltaBeta = geo["maxBeta"], geo["deltaBeta"]
    n_views = int(maxBeta / deltaBeta)             # match CONRAD's (int) truncation
    n_det = int(maxT / deltaT)
    beta = (np.arange(n_views) * deltaBeta)[:, None]                # [views,1]
    t = ((np.arange(n_det) + 0.5) * deltaT - maxT / 2.0)[None, :]   # [1,det] bin centers [mm]
    cb, sb = np.cos(beta), np.sin(beta)
    ax, ay = focal * cb, focal * sb                                 # source     [views,1]
    px, py = t * sb, -t * cb                                        # det. point [views,det]
    dx, dy = px - ax, py - ay                                       # ray direction A->P
    L = np.hypot(dx, dy)

    def chord(cx, cy, r):
        d = np.abs(dx * (ay - cy) + dy * (cx - ax)) / L            # |(C-A) x dir| / |dir|
        return np.where(d < r, 2.0 * np.sqrt(np.maximum(r * r - d * d, 0.0)), 0.0)

    base = {}
    insert_sum = np.zeros((n_views, n_det))
    for ins in inserts:
        c = chord(ins["center_mm"][0], ins["center_mm"][1], ins["radius_mm"])
        base[ins["name"]] = base.get(ins["name"], np.zeros_like(c)) + c
        insert_sum += c
    body = chord(0.0, 0.0, conrad_phantom.BODY_RADIUS_MM)
    base["water"] = np.maximum(body - insert_sum, 0.0)
    return {k: v.astype(np.float32) for k, v in base.items()}


def _rasterize_aa(inserts, n, vox, ss=8):
    """Anti-aliased rasterisation of the all-cylinder phantom on an n x n grid at
    `vox` mm/px, inserts overriding the water body. Each pixel gets the FRACTION of
    its area inside a material (ss x ss subsamples), so edges are band-limited and
    path lengths are sub-pixel accurate (fixes the hard-0/1 quantisation bias)."""
    fine = n * ss
    coord = ((np.arange(fine) + 0.5) / ss - 0.5 - n / 2.0) * vox    # subpixel centers [mm]
    X, Y = np.meshgrid(coord, coord)                               # [fine,fine], X=col, Y=row
    assigned = np.zeros((fine, fine), bool)
    fine_masks = {}
    for ins in inserts:
        cx, cy = ins["center_mm"]; r = ins["radius_mm"]; nm = ins["name"]
        m = ((X - cx) ** 2 + (Y - cy) ** 2 <= r * r) & ~assigned
        fine_masks.setdefault(nm, np.zeros((fine, fine), bool))[m] = True
        assigned |= m
    body = ((X * X + Y * Y) <= conrad_phantom.BODY_RADIUS_MM ** 2) & ~assigned
    fine_masks.setdefault("water", np.zeros((fine, fine), bool))[body] = True
    return {nm: fm.reshape(n, ss, n, ss).mean(axis=(1, 3)).astype(np.float32)
            for nm, fm in fine_masks.items()}


def _project_conrad_analytic(scene, geo):
    """Exact per-material path-length sinograms via CONRAD's FanBeamAnalyticProjector2D
    (PriorityRayTracer.castRay + accumulatePathLenghtForEachMaterial -> a per-material
    MultiChannelGrid2D). Priority/overlap is handled by the ray tracer. Path lengths
    in mm, sinogram [views, det], geometry identical to conrad_ct.fbp."""
    FAP = CG("edu.stanford.rsl.tutorial.fan").FanBeamAnalyticProjector2D
    proj = FAP(geo["focal"], geo["maxBeta"], geo["deltaBeta"], geo["maxT"], geo["deltaT"])
    mc = proj.projectRayDrivenMaterials(scene)
    mats = proj.getMaterials()                       # List<Material>, index == channel
    return {str(mats[c].getName()): conrad_ct.grid2d_to_np(mc.getChannel(c))
            for c in range(int(mats.size()))}


def project_base_materials(inserts, geo=None, n_pix=512, method="aa", scene=None):
    """CONRAD-native fan-beam base-material sinograms (path lengths [mm]).

    method="aa" (default): anti-aliased fractional-coverage rasterisation projected
      with CONRAD's OpenCL grid projector FanBeamProjector2D.projectRayDrivenCL.
      Sub-pixel-accurate path lengths (fixes the hard-0/1 quantisation bias that
      corrupted low-c_Fe iron dHU) with band-limited edges (low ramp overshoot).
    method="analytic": exact per-material path lengths via CONRAD's
      FanBeamAnalyticProjector2D (ray tracer) when a PrioritizableScene is given;
      falls back to the numpy _disk_chords oracle (disks only) if scene is None.
    method="disk_chords": the numpy closed-form disk-chord oracle directly (used to
      unit-test the CONRAD analytic projector).

    Returns ({material_name: sino[views, det]} in mm, geo), matching the recon
    geometry in `geo`. Inserts override the water body.
    """
    conrad_backend.setup()
    if geo is None:
        geo = conrad_ct.fan_geometry(n_pix=n_pix)
    if method == "disk_chords":
        return _disk_chords(inserts, geo), geo
    if method == "analytic":
        return (_project_conrad_analytic(scene, geo) if scene is not None
                else _disk_chords(inserts, geo)), geo
    if method != "aa":
        raise ValueError(f"unknown projection method {method!r} (use 'aa'/'analytic'/'disk_chords')")
    n = int(round(geo["maxT"] / geo["deltaT"]))          # detector-matched, 1 px = deltaT mm
    vox = geo["deltaT"]
    masks = _rasterize_aa(inserts, n, vox)
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
