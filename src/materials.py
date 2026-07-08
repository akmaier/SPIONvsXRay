"""M1 — X-ray materials for the SPION CT study.

Uses CONRAD's material database via pyconrad for mass attenuation coefficients
(cm^2/g), and models the iron-loaded tumor by the mixture rule:

    mu_tumor(E) = rho_soft * (mu/rho)_soft(E) + c_Fe[g/cm^3] * (mu/rho)_Fe(E)

i.e. iron is *added* to the soft-tissue matrix (not diluted in water), so a
zero-iron tumor equals the background exactly (SPEC §5.2). Soft tissue is
approximated by water (CONRAD ships no ICRU-44 soft-tissue XML); this is a
standard CT proxy and keeps the c=0 baseline clean.
"""
from __future__ import annotations
import numpy as np

import conrad_backend
from config import C_FORM_LEVELS, tumor_iron_conc, SPECTRUM

# Densities [g/cm^3]
RHO_SOFT = 1.0      # water proxy for soft tissue
RHO_WATER = 1.0

# The particles are iron OXIDE (magnetite Fe3O4), not pure iron. Per gram of
# iron, magnetite carries this mass of bound oxygen (4*O / 3*Fe):
_M_FE, _M_O = 55.845, 15.999
MAGNETITE_O_PER_FE = (4.0 * _M_O) / (3.0 * _M_FE)   # = 0.382 g O per g Fe


def oxide_contrast_massatten(energies_kev: np.ndarray) -> np.ndarray:
    """Effective mass attenuation [cm^2/g] of magnetite PER GRAM OF IRON.

    = (mu/rho)_Fe + 0.382*(mu/rho)_O, so the tumor contrast for iron mass
    concentration c_Fe is  d_mu = c_Fe * oxide_contrast_massatten(E).
    Oxygen is nearly tissue-equivalent, so it adds only a few % over pure Fe.
    """
    return (mass_attenuation("iron", energies_kev)
            + MAGNETITE_O_PER_FE * mass_attenuation("oxygen", energies_kev))

_AT = None
_DB = None
_DENS: dict = {}


def _api():
    global _AT, _DB
    conrad_backend.setup()
    if _DB is None:
        _DB = conrad_backend.class_getter(
            "edu.stanford.rsl.conrad.physics.materials.database").MaterialsDB
        _AT = conrad_backend.class_getter(
            "edu.stanford.rsl.conrad.physics.materials.utils").AttenuationType
    return _DB, _AT


def density(material_name: str) -> float:
    if material_name not in _DENS:
        DB, _ = _api()
        _DENS[material_name] = float(DB.getMaterial(material_name).getDensity())
    return _DENS[material_name]


def linear_attenuation(material_name: str, energies_kev: np.ndarray) -> np.ndarray:
    """Linear attenuation coefficient mu [1/cm] at the material's own density.

    NOTE: CONRAD's Material.getAttenuation(E, TYPE) returns LINEAR attenuation
    [1/cm] (verified: water@60keV=0.206 with rho=1; iron@60keV=9.48=0.206_massx7.87).
    """
    DB, AT = _api()
    mat = DB.getMaterial(material_name)
    at = AT.TOTAL_WITH_COHERENT_ATTENUATION
    return np.array([float(mat.getAttenuation(float(e), at)) for e in energies_kev])


def mass_attenuation(material_name: str, energies_kev: np.ndarray) -> np.ndarray:
    """Total mass attenuation coefficient (mu/rho) [cm^2/g] = linear / density."""
    return linear_attenuation(material_name, energies_kev) / density(material_name)


def energy_grid() -> np.ndarray:
    return np.arange(SPECTRUM.e_min_kev, SPECTRUM.e_max_kev + 1e-6, SPECTRUM.e_delta_kev)


def linear_attenuation_soft(energies_kev: np.ndarray) -> np.ndarray:
    """Linear attenuation of the soft-tissue matrix [1/cm]."""
    return RHO_SOFT * mass_attenuation("water", energies_kev)


def linear_attenuation_bone(energies_kev: np.ndarray) -> np.ndarray:
    DB, _ = _api()
    rho = float(DB.getMaterial("bone").getDensity())
    return rho * mass_attenuation("bone", energies_kev)


def linear_attenuation_tumor(c_form: float, energies_kev: np.ndarray) -> np.ndarray:
    """Linear attenuation [1/cm] of the iron-loaded tumor for a formulation conc.

    Mixture rule: soft-tissue matrix + iron at c_Fe (mg Fe/ml == 1e-3 g/cm^3).
    """
    c_fe_gcm3 = 1e-3 * tumor_iron_conc(c_form)          # mg/ml -> g/cm^3 iron
    mu_soft = linear_attenuation_soft(energies_kev)
    # magnetite (Fe3O4): iron + bound oxygen, per gram of iron
    return mu_soft + c_fe_gcm3 * oxide_contrast_massatten(energies_kev)


def hu(mu: np.ndarray, mu_water: np.ndarray) -> np.ndarray:
    """Hounsfield units relative to water."""
    return 1000.0 * (mu - mu_water) / mu_water


def sanity_check(outdir: str = None):
    """Compute mu(E) and HU-vs-concentration; save figures + a table."""
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if outdir is None:
        outdir = str(conrad_backend.REPO_ROOT / "results" / "materials")
    os.makedirs(outdir, exist_ok=True)

    E = energy_grid()
    mu_water = linear_attenuation_soft(E)
    mu_bone = linear_attenuation_bone(E)

    # sanity vs NIST mass attenuation (mu/rho) [cm^2/g] at 60 keV
    fe60 = float(mass_attenuation("iron", np.array([60.0]))[0])
    w60 = float(mass_attenuation("water", np.array([60.0]))[0])
    print(f"[check] (mu/rho) iron@60keV = {fe60:.4f} cm^2/g (NIST ~1.205); "
          f"water@60keV = {w60:.4f} (NIST ~0.206)")

    # mu(E) for water, bone, and tumor at a few concentrations
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(E, mu_water, "k-", lw=2, label="soft tissue (water)")
    ax.plot(E, mu_bone, "-", color="#888", lw=1.5, label="cortical bone")
    for c in [1.0, 10.0, 20.0]:
        ax.plot(E, linear_attenuation_tumor(c, E), lw=1.5,
                label=f"tumor c_form={c:g} ({tumor_iron_conc(c):.3f} mg Fe/ml)")
    ax.set_xlabel("energy [keV]"); ax.set_ylabel(r"linear attenuation $\mu$ [1/cm]")
    ax.set_yscale("log"); ax.set_title("Material attenuation vs energy")
    ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(f"{outdir}/mu_vs_energy.png", dpi=130); plt.close(fig)

    # HU vs concentration at the spectrum mean energy (approx effective energy)
    E_eff = 0.5 * (SPECTRUM.e_min_kev + SPECTRUM.e_max_kev)  # placeholder until M3 spectrum
    Ee = np.array([E_eff])
    muw = linear_attenuation_soft(Ee)
    rows, hus = [], []
    for c in C_FORM_LEVELS:
        mut = linear_attenuation_tumor(c, Ee)
        h = float(hu(mut, muw)[0])
        hus.append(h)
        rows.append((c, tumor_iron_conc(c), h))
    fig, ax = plt.subplots(figsize=(6, 4))
    cfe = [tumor_iron_conc(c) for c in C_FORM_LEVELS]
    ax.plot(cfe, hus, "o-", color="#c0392b")
    ax.set_xlabel("tumor iron concentration [mg Fe/ml]")
    ax.set_ylabel(f"tumor HU @ {E_eff:.0f} keV")
    ax.set_title("Predicted tumor contrast vs iron load")
    ax.grid(alpha=0.3); fig.tight_layout()
    fig.savefig(f"{outdir}/hu_vs_conc.png", dpi=130); plt.close(fig)

    with open(f"{outdir}/hu_table.csv", "w") as f:
        f.write("c_form_mg_per_ml,c_fe_mg_per_ml,tumor_HU_at_Eeff\n")
        for c, cfe_, h in rows:
            f.write(f"{c},{cfe_:.4f},{h:.3f}\n")
    print("[ok] wrote", outdir, "-> mu_vs_energy.png, hu_vs_conc.png, hu_table.csv")
    print("[HU vs c_Fe]", [(round(c, 1), round(h, 2)) for c, _, h in rows])
    return rows


if __name__ == "__main__":
    sanity_check()
