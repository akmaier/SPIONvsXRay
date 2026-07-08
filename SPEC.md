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

## 5. Finalized Experimental Configuration

### 5.1 Phantom (rabbit-scale, in-vivo-like)

- **Body:** homogeneous **ICRU soft-tissue** cylinder, ~10–12 cm diameter
  (rabbit trunk), comfortably inside the **20 cm FOV**.
- **Bone insert:** a cortical-bone rod (spine/rib surrogate) to create a genuine
  **beam-hardening** source — this makes the BH-correction on/off comparison
  meaningful.
- **Tumor:** a single **8 cm³** sphere (radius ≈ 1.24 cm, ⌀ ≈ 2.48 cm) of
  iron-loaded soft tissue holding the scan's delivered SPION dose spread
  uniformly (see §5.2); sits in soft tissue offset from the bone (melanoma-like
  accumulation site).
- **Reference ROI:** adjacent tumor-free soft tissue for ΔHU / CNR contrast.

### 5.2 Materials & dose model

- **Core model:** magnetite **Fe₃O₄** (ρ ≈ 5.17 g/cm³, **Fe mass fraction
  0.724**). "6 mg of SPIONs" = whole-particle mass; iron content = 0.724 × mass.
- **Attenuation uses the OXIDE, not pure Fe:** contrast is modeled as magnetite —
  per gram of iron, `Δμ = c_Fe·[(μ/ρ)_Fe + 0.382·(μ/ρ)_O]` (0.382 g bound O per
  g Fe). Oxygen is near tissue-equivalent, so it adds ~+7% over pure iron
  (realistic-dose tumor HU 2.85 → 3.04) — small but included for correctness.
- **Iron-loaded tumor:** the tumor host medium is the **same ICRU soft tissue**
  as the phantom background, with Fe **added** (density raised by 0.001·c_Fe
  g/cm³, mass negligible; high-Z iron drives the contrast). The iron is **not
  diluted in water** — this keeps **c_Fe = 0 ≡ background (ΔHU = 0)**, so all
  contrast is attributable to iron and the tumor-free soft tissue is a true
  iron-free reference.

**Dose model — delivered mass, not a fixed suspension concentration.** The
independent variable is the article's **formulation concentration** `c_form`
[mg SPION/ml]. The delivered SPION mass scales with it, anchored so
`c_form = 10 mg/ml → 6 mg` delivered into the 8 cm³ tumor:

```
m_SPION   = 6 mg × (c_form / 10)          # delivered particle mass
c_tumor   = m_SPION / 8 cm³               # tumor particle concentration
c_Fe      = 0.724 × c_tumor               # tumor IRON concentration (drives X-ray)
          = 0.0543 × c_form  [mg Fe/ml]
```

| `c_form` (mg/ml) | delivered SPION (mg) | tumor SPION (mg/ml) | **tumor Fe (mg/ml)** |
|---:|---:|---:|---:|
| 0    | 0    | 0      | **0**     |
| 0.5  | 0.3  | 0.0375 | **0.027** |
| 1    | 0.6  | 0.075  | **0.054** |
| 2    | 1.2  | 0.15   | **0.109** |
| 5    | 3.0  | 0.375  | **0.271** |
| 10   | 6.0  | 0.75   | **0.543** |
| 20   | 12.0 | 1.5    | **1.086** |

**Note:** even the top dose yields only ~0.5 mg Fe/ml in the tumor — roughly an
order of magnitude below iodine CT enhancement (~2–15 mg/ml). Quantifying whether
this is detectable is a central outcome of the study.

### 5.3 Spectrum & dose

- **Spectrum:** CONRAD **standard** polychromatic X-ray spectrum (nominal
  ~120 kVp C-arm); exact min/max/delta energy + filtration to be read from the
  installed default and recorded at M3.
- **Dose:** **70 000 photons per detector pixel** (unattenuated I₀); Poisson
  noise applied at the projection level.

### 5.4 Detectors (both simulated)

- **EID** — energy-integrating: photon energies weight the signal.
- **PCD** — photon-counting, **energy-resolved multi-bin**. Default **3 bins**
  spanning the spectrum (thresholds set at ≈⅓/⅔ of the energy range once the
  standard spectrum is known), plus the summed count image. Detectability is
  evaluated per-bin and via optimal energy weighting.

### 5.5 Geometry & reconstruction

- **Geometry:** standard C-arm cone-beam (CONRAD default
  `Configuration`/`Trajectory`), **500 projections**, **20 cm FOV**; SID/SDD,
  detector pixel size/count and angular range taken from the default and recorded
  at M4.
- **Reconstruction:** cone-beam **FDK**, isotropic voxels sized to the FOV
  (e.g. 512³ over 20 cm ≈ 0.39 mm), with **water beam-hardening correction
  toggled off/on**.

### 5.6 Factorial design (the "effects study")

| Factor | Levels | n |
|--------|--------|---|
| Formulation conc. `c_form` [mg SPION/ml] (→ tumor Fe, see §5.2) | 0, 0.5, 1, 2, 5, 10, 20 | 7 |
| Detector | EID, PCD (multi-bin) | 2 |
| Beam-hardening correction | off, on | 2 |
| Noise realization @ 70 000 ph/px | repeats | R = 10 |

**Run counts** (concentration lives in the phantom, so it is a per-scan factor;
BH-correction and noise are post-projection steps):

- **Polychromatic forward projections (compute-heavy, 500 views each):**
  `7 conc × 2 detectors = 14`.
- **Noisy sinogram realizations:** `14 × 10 = 140`.
- **Reconstructions (× BH off/on):** `140 × 2 = 280`.
- **Noise-free reference volumes (contrast ceiling):** `14 × 2 = 28`.
- **➡ ≈ 308 analyzed volumes from 14 forward simulations.**

Total scales linearly with the concentration list and `R`.

### 5.7 Detectability metrics

Per tumor ROI vs. soft-tissue reference, for every factor cell:
- **ΔHU** = mean(tumor) − mean(reference).
- **CNR** = ΔHU / σ(reference), from the R noise realizations.
- **Detection threshold:** lowest `c` with **CNR ≥ 3–5** (Rose) and a reported
  minimum meaningful ΔHU — reported separately for EID vs each PCD bin, and for
  BH-correction off vs on.

### 5.8 Spectral optimization (filters, kVp, PCD thresholds)

Iron has **no usable K-edge** (7.1 keV, far below the diagnostic window), so all
contrast is **photoelectric** and lives at low energy. Computed with
`src/spectral.py` (matched-filter / optimal energy-weighting framework, fixed
70 000 air photons/pixel, realistic 6 mg dose through ~10 cm tissue):

- **Optimal monochromatic energy ≈ 30 keV** (lower = more contrast, but the body
  extinguishes flux below ~30 keV).
- **kVp:** lower is better — 60 kVp gives **1.26×**, 80 kVp **1.15×** the ideal
  CNR of the 120 kVp baseline (fixed air flux).
- **Filters:** hardening filters HURT iron contrast — Cu 0.3 mm → 0.70×,
  Sn 0.5 mm → **0.42×**. (They are for spectral *separation*, not low-Z contrast.)
  Softening / thin quasi-mono filters give little at fixed air flux; the win is
  keeping low energy, not adding a K-edge filter.
- **Detector weighting is the biggest lever:** EID captures only **51%** of the
  ideal CNR; an ideal photon-counter with **optimal low-energy weighting ≈ 2×**
  the EID CNR (≈ 4× dose-equivalent).
- **Optimal PCD thresholds (120 kVp, through rabbit):** 2-bin split at **40 keV**
  (1.78× vs EID, 90% of ideal); 3-bin at **35 / 47.5 keV** (1.87× vs EID, 95%).
  Thresholds isolate the photoelectric-rich low band from the Compton high band;
  **not** placed at a K-edge. Refine against the real spectrum at M3.

Added study dimension: sweep `KVP_LEVELS` and `FILTER_CONFIGS` (config.py) and
report EID vs PCD-unweighted vs PCD-optimal-weighted CNR.

### 5.9 Study B — vascular (vessel) tumor phantom

Realistic alternative to the homogeneous tumor: the contrast agent resides in
**tumor vasculature — 150 µm vessels occupying 10 % of the tumor volume** — not
uniformly dissolved. Same delivered iron mass, now confined to 10 % of the
volume → **10× local vessel concentration**.

Key physics: the CT voxel (~390 µm) **cannot resolve** 150 µm vessels, so each
tumor voxel partial-volume-averages vessel + tissue. With mass conserved the
**mean tumor iron per voxel is unchanged**, so first-order the mean ΔHU/CNR
matches Study A. The study tests the **second-order** effects that the
homogeneous model misses:
1. **Beam-hardening nonlinearity** — rays through 10× iron harden more than
   proportionally, so the projection mean ≠ projection of the mean (interacts
   with the BH-correction factor).
2. **Structural noise** — per-voxel vessel fraction varies (binomial), adding
   texture that can degrade detectability.
3. **Resolution/partial-volume** sensitivity.

Modeling: build the tumor as a **sub-resolution voxelized two-compartment
texture** (10 % vessel voxels at 10× conc, 90 % iron-free tissue), forward-
project at fine resolution, reconstruct at the CT grid, then apply the **same
detectability analysis**. Study B runs the **full factorial** (7 conc × 2 det ×
2 BH × 10 noise) **× all filters/kVp**, mirroring Study A, so the two phantoms
are directly comparable. (Requires the M4–M6 projector/FDK pipeline.)

## 6. Deliverables

- `src/materials.py` — build/register SPION materials for a list of concentrations.
- `src/phantom.py` — construct the rabbit-scale soft-tissue + bone phantom with
  the 8 cm³ SPION tumor insert.
- `src/simulate.py` — spectrum + 500-projection C-arm forward projection.
- `src/reconstruct.py` — FDK reconstruction.
- `src/analyze.py` — ROI HU/attenuation, contrast curve, CNR, detectability.
- `run_experiment.py` — end-to-end driver; writes figures + a results table.
- `src/build_dashboard.py` — regenerate the GitHub Pages audit dashboard
  (`docs/`) from `results/`: copies figures into `docs/assets/`, fills the
  detectability table, and stamps the commit hash.
- `results/` — projection stacks, reconstructions, plots (git-ignored).
- `docs/` — public **audit dashboard** (GitHub Pages) showing spectra, materials,
  phantom, projections, reconstructions, detectability curves, and annotated
  code snippets for auditing. Committed (images live in `docs/assets/`).
- `LICENSE` — MIT.
- Findings summarized back into `README.md` / `DEVLOG.md`.

## 7. Milestones

- **M0 — Environment**: create a Python 3.9/3.10 venv, `pip install pyconrad`,
  `setup_pyconrad()` succeeds, JVM up, CONRAD.jar present. *(de-risks Mac/3.14.)*
- **M1 — Materials**: register water + SPION@c materials; sanity-check linear
  attenuation vs. c at a fixed energy against NIST/known Fe values.
- **M2 — Phantom**: build rabbit-scale soft-tissue + bone phantom with the 8 cm³
  iron-loaded tumor; verify material map visually.
- **M3 — Forward model**: standard spectrum + single polychromatic projection.
- **M4 — Full scan**: 500-projection C-arm acquisition.
- **M5 — Reconstruction**: FDK volume, correct geometry/scaling.
- **M6 — Analysis**: HU-vs-concentration curve, CNR, detectability threshold.
- **M7 — Noise study**: repeat with photon noise; report detection limit.
- **M8 — Dashboard**: `src/build_dashboard.py` populates the `docs/` audit
  dashboard from `results/` (figures + detectability table + commit stamp);
  GitHub Pages serves it publicly. *(Scaffold + Pages already live.)*

## 8. Status of Open Questions

**Resolved** (user, 2026-07-08): dose = 70 000 ph/px; tumor = 8 cm³; phantom =
rabbit-scale soft-tissue **+ bone** insert; FOV = 20 cm; detectors = **both**
EID and multi-bin PCD; beam-hardening correction = **both** off/on;
detectability = **both** CNR and ΔHU; concentrations = 0/0.5/1/2/5/10/20
mg Fe/ml. → **Sufficient to implement.**

**Defaults recorded at implementation time** (non-blocking, will be read from the
installed CONRAD and logged): standard-spectrum kVp/filtration & energy sampling
(M3), C-arm SID/SDD/detector/angular range (M4), PCD bin thresholds (from the
spectrum range), reconstruction voxel size, iron-oxide polymorph fixed to Fe₃O₄
for reporting.
