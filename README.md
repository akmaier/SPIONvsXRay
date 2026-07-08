# SPIONvsXRay

**Can superparamagnetic iron-oxide nanoparticles (SPIONs) be seen in CT?**

This project simulates the X-ray absorption and cone-beam CT reconstruction of
PAA-coated iron-oxide nanoparticle suspensions at varying iron concentrations,
to determine whether — and at what concentration — they become visible in a
standard C-arm CT. The simulation uses the [CONRAD](https://www5.cs.fau.de/conrad/)
framework via [`pyconrad`](https://pypi.org/project/pyconrad/).

The nanoparticles are those described in Heinen et al., *Ultrasonics
Sonochemistry* 130 (2026) 107876 (PAA-coated SPION clusters, iron-oxide cores
~8–12 nm, used at ~1–10 mg Fe/ml). *The article PDF is git-ignored (copyright).*

## Idea in one picture

```
SPION suspensions          multi-vial phantom        C-arm cone-beam scan
@ 0,1,2,5,10 mg Fe/ml  ──►  (water + inserts)   ──►  500 projections, std spectrum
                                                             │
                                          HU vs. concentration, CNR  ◄── FDK reco
```

Because X-ray attenuation is dominated by **iron content**, the tumor is modeled
as **iron-loaded soft tissue** — the same soft-tissue matrix as the rabbit
background with Fe added at `c` mg/ml (not diluted in water), so a zero-iron
tumor is indistinguishable from background and all contrast is attributable to
iron. The phantom is imaged polychromatically, reconstructed, and the contrast
(Hounsfield Units) and contrast-to-noise ratio are measured in the tumor ROI.

## Method summary

1. **Materials** — register a CONRAD material per iron concentration.
2. **Phantom** — a water cylinder with rod inserts, one concentration each.
3. **Spectrum** — CONRAD standard polychromatic X-ray spectrum.
4. **Acquisition** — 500 projections in the default C-arm cone-beam geometry.
5. **Reconstruction** — cone-beam FDK.
6. **Analysis** — HU-vs-concentration curve, CNR, detectability threshold.

See [`SPEC.md`](SPEC.md) for the full plan and milestones, and
[`DEVLOG.md`](DEVLOG.md) for progress.

## Status

🚧 Planning complete; environment setup next. This repo currently contains the
plan (`SPEC.md`), this README, and the development log (`DEVLOG.md`). Code lands
under `src/` per the milestones in the SPEC.

## Getting started (planned)

CONRAD requires **Java 8** (already present) and a pyconrad-compatible Python
(3.9/3.10 recommended — the system Python 3.14 is too new for JPype):

```bash
python3.10 -m venv .venv && source .venv/bin/activate
pip install pyconrad
python -c "import pyconrad; pyconrad.setup_pyconrad()"   # downloads CONRAD.jar, starts JVM
```

## Experimental design (finalized)

A factorial **effects study** on an in-vivo-like rabbit phantom (soft tissue +
bone) with a single **8 cm³** SPION tumor, **20 cm FOV**, **500 projections**,
**70 000 photons/pixel**:

| Factor | Levels |
|--------|--------|
| Formulation conc. [mg SPION/ml] | 0, 0.5, 1, 2, 5, 10, 20 (7) |
| Detector | energy-integrating (EID), photon-counting multi-bin (PCD) |
| Beam-hardening correction | off, on |
| Noise realizations @ 70 000 ph/px | 10 |

**Dose model:** a *delivered mass*, not a fixed vial concentration — anchored at
**6 mg SPIONs for the 10 mg/ml formulation**, spread over the 8 cm³ tumor and
scaled with concentration. With a magnetite core (72.4% Fe), tumor iron works out
to `c_Fe = 0.0543 × c_form`, so the top dose is only **~0.5 mg Fe/ml** — ~10×
below iodine CT enhancement. Whether that is detectable is a key study outcome.
See [`SPEC.md`](SPEC.md) §5.2 for the full conversion table.

→ **14 forward simulations → ≈ 308 reconstructed/analyzed volumes.**
Detectability reported as **both** ΔHU and **CNR** (Rose CNR ≥ 3–5), per detector
and per BH-correction state. See [`SPEC.md`](SPEC.md) §5 for full parameters.

Remaining defaults (non-blocking, read from the installed CONRAD and logged at
implementation): standard-spectrum kVp/filtration, C-arm SID/SDD/detector/angular
range, PCD bin thresholds, reconstruction voxel size (Fe₃O₄ assumed for
reporting).

## Repository layout

```
SPEC.md          Full experimental plan and milestones
README.md        This file
DEVLOG.md        Running development log
paper_plan.md    SPIE Medical Imaging paper plan + Maier writing-style analysis
paper/template/  Official SPIE Proceedings LaTeX template (spie.cls v3.4, verified)
src/             Simulation & analysis code (to be added)
results/         Generated outputs (git-ignored)
```

## License / attribution

Simulation built on CONRAD (FAU Pattern Recognition Lab). The reference article
PDF is not redistributed (git-ignored).
