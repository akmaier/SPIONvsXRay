# SPEC — Simulating the CT Visibility of Magnetic Nanoparticles (SPIONs)

## 1. Motivation & Research Question

The reference article (Heinen et al., *Ultrasonics Sonochemistry* 130 (2026) 107876)
describes **polyacrylic-acid (PAA)-coated superparamagnetic iron-oxide
nanoparticles (SPIONs)** engineered for biomedical use (magnetic drug targeting,
cell seeding/tracking, hyperthermia). The nanoparticles have iron-oxide cores
(~8 nm and ~12 nm) and are used at **iron concentrations of ~1–10 mg Fe/ml**;
cellular uptake is on the order of ~8 pg Fe/cell.

**Central question:** *At which iron concentrations do these SPIONs produce
enough X-ray attenuation contrast to be visible in a standard C-arm cone-beam
CT?*

We answer this by **simulating** X-ray absorption and CT reconstruction of SPION
suspensions at several concentrations, using the CONRAD framework via
[`pyconrad`](https://pypi.org/project/pyconrad/), and quantifying the resulting
contrast (Hounsfield Units / attenuation) and detectability.

## 2. Scientific Background

- X-ray attenuation of a dilute nanoparticle suspension is driven almost entirely
  by the **iron content** (Z=26, K-edge 7.1 keV); the water matrix and the thin
  PAA/carbon coating contribute negligibly at diagnostic energies.
- The nanoparticle core is an iron oxide — **magnetite Fe₃O₄** (ρ ≈ 5.17 g/cm³,
  Fe mass fraction 0.724) or **maghemite γ-Fe₂O₃** (ρ ≈ 4.9 g/cm³, Fe fraction
  0.700). For the X-ray model, only the delivered **mg Fe/ml** matters; core
  polymorph mainly affects how we convert particle mass ↔ Fe mass.
- A suspension at concentration `c` [mg Fe/ml] is modeled as
  **water + dissolved iron** at the corresponding mass fraction, following the
  CONRAD custom-material "mixture" recipe (analogous to the Ultravist/iopromide
  example in the tutorial).

## 3. Approach Overview

```
Materials  →  Phantom  →  Polychromatic forward projection  →  FDK reco  →  ROI analysis
(SPION @ c)   (vials)     (500 proj, C-arm, std spectrum)      (3D vol)    (HU vs c, CNR)
```

1. **Custom materials** — one CONRAD material per iron concentration
   (`water + Fe @ c mg/ml`) plus a pure-water reference.
2. **Phantom** — a water cylinder ("sample holder") containing several rod/vial
   inserts, each assigned a different SPION concentration (a multi-vial
   contrast-detectability phantom, inspired by the CONRAD numerical phantoms).
3. **Spectrum** — the CONRAD standard polychromatic X-ray spectrum
   (`PolychromaticXRaySpectrum`), per the spectral-absorption tutorial.
4. **Forward projection** — polychromatic line-integral absorption
   (`PolychromaticAbsorptionModel`) rendered into **500 projections** in a
   **standard C-arm cone-beam geometry** (CONRAD default trajectory).
5. **Reconstruction** — cone-beam FDK (filtered back-projection).
6. **Analysis** — mean attenuation / HU in each insert ROI, contrast vs.
   concentration curve, contrast-to-noise ratio (CNR), and a detectability
   threshold (Rose criterion, CNR ≳ 3–5).

## 4. Tools & Environment

| Component | Choice | Status |
|-----------|--------|--------|
| Simulation/reco | CONRAD via **pyconrad 0.8.0** (pip) | to install |
| JVM | **Java 8** (required by CONRAD) | ✓ present (1.8.0_291) |
| Python | pyconrad tested on 3.6–3.10 | ⚠ local is **3.14** — needs a compatible venv |
| Platform | macOS (darwin, Apple Silicon) | ⚠ pyconrad notes Mac caveats |
| Reference impls | `edu.stanford.rsl.tutorial.physics.*` (Java, called from Python) | verify shipped in CONRAD.jar |

CONRAD is Java; pyconrad starts a JVM (`pyconrad.setup_pyconrad()`) and exposes
Java classes to Python (`ClassGetter` / `pyconrad.ClassGetter`), with NumPy
interop for grids. **The exact Python-side API must be verified against the
installed pyconrad version** before writing the pipeline.

## 5. Experimental Configuration (proposed defaults — to confirm, see §7)

- **Concentrations swept:** `0, 1, 2, 5, 10 mg Fe/ml` (matches article range),
  optionally add a low-end `0.5` and a high `20` to bracket the detection limit.
- **Core model:** magnetite Fe₃O₄ for particle-mass conversion; X-ray model uses
  Fe mass fraction directly.
- **Background/matrix:** water (proxy for aqueous suspension / soft tissue).
- **Spectrum:** CONRAD standard spectrum (nominal ~120 kVp C-arm); confirm exact
  min/max/delta energy, peak voltage and filtration from the tutorial defaults.
- **Geometry:** CONRAD default C-arm cone-beam `Configuration`/`Trajectory`,
  **500 projections**; angular range and SID/SDD/detector taken from the default
  (to be recorded once read from config).
- **Reconstruction:** FDK, isotropic voxels; volume sized to the phantom.
- **Noise:** first pass **noise-free** (pure contrast ceiling); second pass add
  Poisson photon noise at a specified dose to assess realistic detectability.

## 6. Deliverables

- `src/materials.py` — build/register SPION materials for a list of concentrations.
- `src/phantom.py` — construct the multi-vial phantom with material labels.
- `src/simulate.py` — spectrum + 500-projection C-arm forward projection.
- `src/reconstruct.py` — FDK reconstruction.
- `src/analyze.py` — ROI HU/attenuation, contrast curve, CNR, detectability.
- `run_experiment.py` — end-to-end driver; writes figures + a results table.
- `results/` — projection stacks, reconstructions, plots (git-ignored).
- Findings summarized back into `README.md` / `DEVLOG.md`.

## 7. Milestones

- **M0 — Environment**: create a Python 3.9/3.10 venv, `pip install pyconrad`,
  `setup_pyconrad()` succeeds, JVM up, CONRAD.jar present. *(de-risks Mac/3.14.)*
- **M1 — Materials**: register water + SPION@c materials; sanity-check linear
  attenuation vs. c at a fixed energy against NIST/known Fe values.
- **M2 — Phantom**: build multi-vial phantom, verify material map visually.
- **M3 — Forward model**: standard spectrum + single polychromatic projection.
- **M4 — Full scan**: 500-projection C-arm acquisition.
- **M5 — Reconstruction**: FDK volume, correct geometry/scaling.
- **M6 — Analysis**: HU-vs-concentration curve, CNR, detectability threshold.
- **M7 — Noise study**: repeat with photon noise; report detection limit.

## 8. Open Questions / Missing Information (blocking or affecting design)

See the "Missing information" section in `README.md`. Key unknowns: iron-oxide
polymorph & density, exact concentration list and units, standard-spectrum kVp
and dose/noise model, phantom & C-arm geometry parameters (SID/SDD, detector,
angular range, object scale), and the quantitative detectability criterion.
