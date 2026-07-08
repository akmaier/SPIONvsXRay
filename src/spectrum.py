"""M3 — real CONRAD polychromatic X-ray spectra.

Wraps edu.stanford.rsl.conrad.physics.PolychromaticXRaySpectrum. The no-arg
constructor is the CONRAD *standard* spectrum (10-150 keV, 0.5 keV steps, mean
~55.4 keV, tungsten anode with characteristic lines). Parameterized spectra
(varying kVp) drive the filter/kVp optimization sweep. Added filtration is
applied via CONRAD material attenuation.
"""
from __future__ import annotations
import numpy as np

import conrad_backend
import materials


def _spectrum_class():
    conrad_backend.setup()
    return conrad_backend.class_getter("edu.stanford.rsl.conrad.physics").PolychromaticXRaySpectrum


def _arrays(sp):
    E = np.array(list(sp.getPhotonEnergies()), dtype=float)
    flux = np.array(list(sp.getPhotonFlux()), dtype=float)
    return E, flux


def standard_spectrum():
    """CONRAD standard spectrum. Returns (E[keV], flux, info dict)."""
    P = _spectrum_class()
    sp = P()
    E, flux = _arrays(sp)
    info = dict(e_min=float(sp.getMin()), e_max=float(sp.getMax()),
                delta=float(sp.getDelta()), kvp=float(sp.getPeakVoltage()),
                e_avg=float(sp.getAveragePhotonEnergy()))
    return E, flux, info


def conrad_spectrum(kvp: float, e_min=10.0, e_max=150.0, delta=0.5, mas=1.0):
    """Parameterized CONRAD spectrum at a given peak voltage. Returns (E, flux)."""
    P = _spectrum_class()
    sp = P(float(e_min), float(e_max), float(delta), float(kvp), float(mas))
    return _arrays(sp)


def apply_filters(E, flux, filters):
    """Attenuate a spectrum by (material_name, thickness_mm) pairs."""
    out = flux.copy()
    for name, t_mm in filters:
        mu = materials.linear_attenuation(name, E)   # 1/cm
        out = out * np.exp(-mu * (t_mm / 10.0))
    return out


def normalized(E, flux):
    tot = np.trapezoid(flux, E)
    return flux / tot if tot > 0 else flux


if __name__ == "__main__":
    E, flux, info = standard_spectrum()
    print("CONRAD standard spectrum:", info)
    print("energies:", E[0], "..", E[-1], "n=", E.size)
    nz = E[flux > flux.max() * 1e-3]
    print(f"effective range (flux>0.1% peak): {nz.min():.0f}-{nz.max():.0f} keV")
    print("peak flux at E =", float(E[np.argmax(flux)]), "keV")
