"""M4 (rev.) — polychromatic fan-beam projection via CONRAD, EID + multi-bin PCD.

Projection uses CONRAD's CPU fan-beam projector (src/conrad_ct.py). Material
decomposition: project the ED-phantom soft/bone/iron component maps ONCE through
CONRAD, then combine polychromatically for any detector/noise:

    tau(E) = L_soft*mu_soft(E) + L_bone*mu_bone(E) + A_iron*oxide(E)

Sinograms are saved to results/sinograms/ (npz).
"""
from __future__ import annotations
import os
import numpy as np

import conrad_backend
import conrad_ct
import phantom
import materials
import spectrum as spec
from config import SPECTRUM, DETECTORS

N0 = SPECTRUM.photons_per_pixel


def project_materials(model="homogeneous", N=512):
    """CONRAD fan-beam projection of the ED-phantom component maps.

    Returns material sinograms in physical units:
      L_soft, L_bone [cm];  A_iron [g/cm^2];  plus the phantom + geometry.
    """
    ph = phantom.build_ed_phantom(model, N=N)
    geo = conrad_ct.fan_geometry(n_pix=N)
    # CONRAD ray integral = value * path[mm]
    L_soft = conrad_ct.project(ph.soft, geo) / 10.0                 # mm -> cm
    L_bone = conrad_ct.project(ph.bone, geo) / 10.0                 # mm -> cm
    A_iron = conrad_ct.project(ph.iron, geo) * 1e-4                 # (mg/cm^3)*mm -> g/cm^2
    return dict(L_soft=L_soft, L_bone=L_bone, A_iron=A_iron, geo=geo, phantom=ph)


def _spectrum(kvp=None, filters=()):
    if kvp is None:
        E, flux, _ = spec.standard_spectrum()
    else:
        E, flux = spec.conrad_spectrum(kvp)
    if filters:
        flux = spec.apply_filters(E, flux, list(filters))
    return E, flux / flux.sum()


def detector_sinograms(proj, kvp=None, filters=(), add_noise=True, seed=0):
    """Material sinograms -> EID + PCD-bin line-integral sinograms (with noise)."""
    rng = np.random.default_rng(seed)
    E, s = _spectrum(kvp, filters)
    mu_soft = materials.linear_attenuation("water", E)
    mu_bone = materials.linear_attenuation("bone", E)
    oxide = materials.oxide_contrast_massatten(E)

    shape = proj["L_soft"].shape
    S_det = np.zeros(shape); S_det_E2 = np.zeros(shape)
    S_air = float(np.sum(N0 * s * E))
    edges = np.array(DETECTORS.pcd_bin_edges_kev)
    nb = len(edges) - 1
    C_det = [np.zeros(shape) for _ in range(nb)]
    C_air = [0.0] * nb

    for i, Ei in enumerate(E):
        tau = mu_soft[i] * proj["L_soft"] + mu_bone[i] * proj["L_bone"] + oxide[i] * proj["A_iron"]
        n = N0 * s[i] * np.exp(-tau)
        S_det += n * Ei
        S_det_E2 += n * Ei * Ei
        b = int(np.searchsorted(edges, Ei, side="right") - 1)
        if 0 <= b < nb:
            C_det[b] += n
            C_air[b] += N0 * s[i]

    if add_noise:
        S_meas = S_det + rng.normal(0.0, np.sqrt(np.maximum(S_det_E2, 1e-30)))
        C_meas = [rng.poisson(np.maximum(c, 0.0)) for c in C_det]
    else:
        S_meas, C_meas = S_det, C_det

    eps = 1e-6
    p_eid = -np.log(np.clip(S_meas, eps, None) / S_air)
    p_pcd = [-np.log(np.clip(cm, eps, None) / max(ca, eps)) for cm, ca in zip(C_meas, C_air)]
    return dict(eid=p_eid, pcd=p_pcd, edges=edges, geo=proj["geo"])


def save_sinograms(model="homogeneous", N=512, outdir=None, seed=1):
    if outdir is None:
        outdir = str(conrad_backend.REPO_ROOT / "results" / "sinograms")
    os.makedirs(outdir, exist_ok=True)
    proj = project_materials(model, N=N)
    det = detector_sinograms(proj, add_noise=True, seed=seed)
    path = f"{outdir}/sino_{model}.npz"
    np.savez_compressed(path, L_soft=proj["L_soft"], L_bone=proj["L_bone"],
                        A_iron=proj["A_iron"], eid=det["eid"],
                        pcd=np.array(det["pcd"]), edges=det["edges"])
    print(f"[ok] saved {path}: material + EID + {len(det['pcd'])}-bin PCD sinograms "
          f"shape {det['eid'].shape}")
    return proj, det


if __name__ == "__main__":
    save_sinograms("homogeneous")
