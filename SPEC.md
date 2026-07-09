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
  PAA/carbon coating contribute negligibly at diagnostic energies. We nonetheless
  **simulate the full particle** (core + coating) for completeness — the coating is
  included as mass but, being low-Z, barely moves μ (see §5.2).
- The nanoparticle is a **PAA-coated magnetite (Fe₃O₄) cluster** (Heinen et al.):
  iron-oxide cores ~8/12 nm assembling into ~70–250 nm hydrodynamic clusters, with a
  **polyacrylic-acid (PAA, monomer (C₃H₄O₂)ₙ) coating**. The core is magnetite
  (ρ ≈ 5.17 g/cm³, Fe mass fraction 0.724) — maghemite γ-Fe₂O₃ (ρ ≈ 4.9, Fe 0.700) is
  the alternative polymorph. For the X-ray model, only the delivered **mg Fe/ml**
  drives contrast; the polymorph mainly affects particle-mass ↔ Fe-mass conversion.
- **Coating fraction is per formulation, from the article's TGA** (supplement
  Table A.1 / Fig A.1). The inorganic (iron-oxide) residual is **87.7 % (SPION I,
  12 nm cores)** and **66.4 % (SPION II, 8 nm cores)**; on heating the residual
  oxidizes magnetite → hematite (+3.4 % mass, article Eq A.1 uses 0.966·w_residual),
  so the original magnetite fraction is `R_m/1.034` = **84.8 % (SPION I)** and
  **64.2 % (SPION II)**. The complement is PAA: **φ ≈ 0.15 (SPION I)**, **φ ≈ 0.36
  (SPION II)** — SPION II carries roughly twice the coating. See §5.2.
- The nanoparticle-loaded tumor at iron concentration `c_Fe` [mg Fe/ml] is modeled as
  **ICRU soft tissue + magnetite (Fe₃O₄) core + PAA coating** added to the tissue
  matrix (NOT diluted in water — see §5.2), so a zero-iron tumor equals the
  soft-tissue background (ΔHU = 0) and all contrast is attributable to iron. Built
  with the CONRAD custom-material "mixture" recipe (PAA registered as a C/H/O
  material).

## 3. Approach Overview

```
Materials  →  Phantom  →  Polychromatic forward projection  →  Fan FBP  →  ROI analysis
(SPION @ c)   (soft-     (500 views, C-arm short scan, 90 kVp) (2D slice)  (HU vs c, CNR)
              tissue)
```

1. **Custom materials** — one CONRAD material per iron concentration (**magnetite
   Fe₃O₄ in ICRU soft tissue**) plus the iron-free soft-tissue reference.
2. **Phantom** — a soft-tissue cylinder with rod/vial inserts, each a different SPION
   concentration, plus a cortical-bone rod (a multi-vial contrast-detectability
   phantom).
3. **Spectrum** — the CONRAD standard **90 kVp** polychromatic X-ray spectrum
   (`PolychromaticXRaySpectrum`).
4. **Forward projection** — per-material path lengths (anti-aliased grid, or exact
   analytic `FanBeamAnalyticProjector2D`), combined polychromatically into EID + PCD
   line integrals over **500 views** of a **C-arm fan short scan**.
5. **Reconstruction** — CONRAD-native **fan FBP** (Parker → cosine → ramp →
   distance-weighted backprojection), always with water precorrection (§5.5).
6. **Analysis** — per-insert iron ΔHU (c0-corrected, local annulus), CNR over noise
   realizations, and a Rose detection threshold (CNR ≳ 3–5).

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

### 5.1 Phantom (round geometric, rabbit-scale)

**Decision 2026-07-08:** use a **single round geometric phantom**; the
realistic-anatomy arm is dropped. Rationale: there is **no digital rabbit
phantom** (ROBY = rat, MOBY = mouse; the Segars/RADAR family has no rabbit — JNM
2010;51:471), and a real rabbit CT was not pursued. The geometric phantom is
clean, reproducible, and isolates the detector/spectral variables. **Paper
limitation:** idealized cylindrical anatomy → idealized beam hardening/scatter;
no organ heterogeneity.

- **Body:** homogeneous **ICRU soft-tissue** (water proxy) cylinder, ~10–12 cm
  diameter (rabbit trunk), inside the **20 cm FOV**.
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
  0.724**). "6 mg of SPIONs" = whole-particle mass; iron content = 0.724 × magnetite
  mass.
- **Attenuation uses the OXIDE, not pure Fe:** contrast is modeled as magnetite —
  per gram of iron, `Δμ = c_Fe·[(μ/ρ)_Fe + 0.382·(μ/ρ)_O]` (0.382 g bound O per
  g Fe). Oxygen is near tissue-equivalent, so it adds ~+7% over pure iron
  (realistic-dose tumor HU 2.85 → 3.04) — small but included for correctness.
- **Full-nanoparticle model (core + PAA coating):** the whole particle is
  magnetite core + a **polyacrylic-acid (PAA, monomer (C₃H₄O₂)ₙ) coating**, so the
  tumor material per ml = soft-tissue matrix + magnetite (mass = c_Fe/0.724) + PAA
  (mass = φ·c_NP). From the tumor iron concentration we compute the
  **whole-nanoparticle** concentration
  `c_NP = c_Fe / (0.724·(1 − φ))`, where **φ = PAA coating mass fraction**.
  - **φ is per formulation, from the article's TGA** (supplement Table A.1 / Fig A.1),
    replacing the earlier single estimate. The inorganic (iron-oxide) residual R_m is
    **87.7 % (SPION I, 12 nm cores)** and **66.4 % (SPION II, 8 nm cores)**; the residual
    oxidizes magnetite → hematite (+3.4 % mass; article Eq A.1 uses 0.966·w_residual),
    so the original magnetite fraction = R_m/1.034 = **84.8 % (SPION I)** ⇒ coating
    **φ ≈ 0.15**, and **64.2 % (SPION II)** ⇒ coating **φ ≈ 0.36**. SPION II carries
    much more PAA than SPION I, so each formulation's full-particle model uses **its own
    coating**: `c_NP ≈ 1.625·c_Fe` (SPION I), `c_NP ≈ 2.16·c_Fe` (SPION II). *(The
    earlier central estimate φ = 0.15 turned out ≈ exact for SPION I.)*
  - **φ sets ONLY the reported mg SPION/ml** (particle basis). The PAA is low-Z
    (C/H/O), tissue/water-equivalent, so it is **negligible for μ**: registering it and
    adding its mass moves the monochromatic iron ΔHU at the top dose by only ≈ 3.6 %
    at φ = 0.15 (3.62 → 3.75 HU @ 60 keV). *Iron* (mg Fe/ml) and *particle* (mg SPION/ml)
    concentrations are reported distinctly.
- **Nanoparticle-loaded tumor:** the tumor host medium is the **same ICRU soft
  tissue** as the phantom background, with the **full particle** (magnetite + PAA)
  **added** (density raised by ~0.001·c_Fe g/cm³, mass negligible; high-Z iron drives
  the contrast). The particle is **not diluted in water** — this keeps
  **c_Fe = 0 ≡ background (ΔHU = 0)**, so all contrast is attributable to iron and the
  tumor-free soft tissue is a true iron-free reference.

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

Reported concentrations: **tumor Fe** (mg Fe/ml, drives μ) and the **whole-particle**
`c_NP` (mg SPION/ml) with its PAA-coating mass, using the formulation's own coating
— `c_NP = 1.625·c_Fe` for SPION I (φ ≈ 0.15), `c_NP = 2.16·c_Fe` for SPION II
(φ ≈ 0.36). The delivered-mass table below is tabulated at φ = 0.15 (SPION I).

| `c_form` (mg/ml) | delivered SPION (mg) | **tumor Fe (mg/ml)** | tumor c_NP (mg SPION/ml) | tumor PAA (mg/ml) |
|---:|---:|---:|---:|---:|
| 0    | 0    | **0**     | 0      | 0      |
| 0.5  | 0.3  | **0.027** | 0.044  | 0.0066 |
| 1    | 0.6  | **0.054** | 0.088  | 0.0132 |
| 2    | 1.2  | **0.109** | 0.176  | 0.0265 |
| 5    | 3.0  | **0.271** | 0.441  | 0.0662 |
| 10   | 6.0  | **0.543** | 0.882  | 0.1324 |
| 20   | 12.0 | **1.086** | 1.765  | 0.2647 |

**Note:** even the top dose yields only ~0.5 mg Fe/ml in the tumor — roughly an
order of magnitude below iodine CT enhancement (~2–15 mg/ml). Quantifying whether
this is detectable is a central outcome of the study.

**Concentration basis — cellular-loading anchor.** The tumor iron band is grounded
in the article's **measured cellular loading** (Heinen et al., Fig 5, B16-F10
melanoma, all configurations incl. fresh; 0 h / 24 h): SPION I-113 nm
**8.23 / 3.78**, SPION II-115 **3.86 / 1.52**, II-98 **3.60 / 1.37**, II-76
**2.51 / 0.85** pg Fe/cell. At a tumor cell density ~10⁸–10⁹ cells/cm³ this gives
~**0.25–8 mg Fe/ml** (realistic ~1–2.5 mg Fe/ml for SPION I), which brackets the
swept tumor concentrations. This *justifies* the range; it does not change the
delivered-mass sweep above.

### 5.3 Spectrum & dose

- **Spectrum:** CONRAD **standard** polychromatic X-ray spectrum, **90 kVp**
  (`PolychromaticXRaySpectrum`, W anode, 10–150 keV @ 0.5 keV, mean ~55 keV), fixed
  for the detectability study. *(An earlier draft said "~120 kVp" — that was wrong;
  kVp/filtration are explored separately in §5.8.)*
- **Dose:** a study **factor** — **low = 20 000 / high = 100 000** unattenuated
  photons per detector pixel (I₀), bracketing a realistic stationary C-arm CBCT
  low/high-dose span (~0.15–1.2 µGy/frame at the detector; the earlier fixed 70 000
  sits mid/high). Poisson (PCD) / Gaussian-quantum (EID) noise applied at the
  projection level.

### 5.4 Detectors (both simulated)

- **EID** — energy-integrating: photon energies weight the signal.
- **PCD** — photon-counting, **energy-resolved multi-bin**. Default **3 bins**
  spanning the spectrum (thresholds set at ≈⅓/⅔ of the energy range once the
  standard spectrum is known), plus the summed count image. Detectability is
  evaluated per-bin and via optimal energy weighting.

### 5.5 Geometry & reconstruction

The detectability study reconstructs a **2D fan-beam** central slice (the fan is the
mid-plane of the C-arm cone), using CONRAD's `tutorial.fan` classes end-to-end.

- **Geometry:** SID = 750 mm (= CONRAD `focalLength`), SDD = 1200 mm, **500 views**
  over a **short scan** (180°+2γ, Parker-weighted). CONRAD's virtual-detector
  convention: detector at the isocenter, sampled at **deltaT = 1 mm** (kept at 1 mm
  so the ramp kernels' `deltaS` scaling is exact — SheppLogan carried a 1/deltaS²
  bug, fixed upstream this session), `maxT` = recon FOV. Recon **512×512 @ 0.5 mm**
  (256 mm FOV).
- **Reconstruction (CONRAD-native fan FBP):** forward projection — anti-aliased grid
  `FanBeamProjector2D.projectRayDrivenCL` (default `aa`) or exact analytic
  `FanBeamAnalyticProjector2D` (per-material, contributed upstream this session) →
  `ParkerWeights` (short scan) → `CosineFilter` → `SheppLoganKernel` (ramp + roll-off)
  → **distance-weighted** `FanBeamBackprojector2D` (the 1/U² weight was reinstated
  upstream this session). **Water beam-hardening precorrection (`WaterPrecorrectionTool`)
  is always applied.**
- **3D cone-beam FDK** (`ConeBeamProjector`/`ConeBeamBackprojector`, verified this
  session: DC scale-invariant across voxel/detector, only axial cone shading) is a
  **separate track** (RabbitCT geometry), not the detectability reconstructor.

### 5.6 Detectability study — factor structure

A multi-factor effects study (NOT a mathematical factorial). Factors:

| Factor | Levels | n |
|--------|--------|---|
| Formulation conc. `c_form` [mg SPION/ml] (→ tumor Fe, see §5.2) | 0, 0.5, 1, 2, 5, 10, 20 | 7 |
| Detector | EID, PCD (multi-bin) | 2 |
| Bone (cortical-bone rod) | absent, present | 2 |
| Dose (unattenuated I₀, ph/px) | low 20 000, high 100 000 | 2 |
| Tumor distribution | uniform, vessel (§5.9) | 2 |
| Noise realizations | repeats | R = 30 |

→ 2 det × 2 bone × 2 dose × 2 distribution × 6 iron levels (c0 = reference) =
**96 cells × R**.

**Held fixed — no longer factors:**
- **Water beam-hardening precorrection: ALWAYS on** (via CONRAD `WaterPrecorrectionTool`,
  §5.5). On the homogeneous phantom it is ~inert against the local-annulus metric; it
  matters once bone streaking / multi-material correction is in play (discussion point).
- **Short scan: ALWAYS on** (Parker-weighted 180°+2γ).
- **Spectrum: fixed** at the 90 kVp CONRAD standard; kVp/filtration is the separate
  spectral-optimization analysis (§5.8), not a factor here.

### 5.7 Detectability metrics

Per tumor ROI vs. soft-tissue reference, for every factor cell:
- **ΔHU** = mean(tumor) − mean(reference).
- **CNR** = ΔHU / σ(reference), from the R noise realizations.
- **Detection threshold:** lowest `c` with **CNR ≥ 3–5** (Rose) and a reported
  minimum meaningful ΔHU — reported per cell: detector (EID vs PCD), bone
  (absent/present), dose (low/high), and distribution (uniform/vessel). Water
  precorrection is always on, so it is not a reported split.

### 5.8 Spectral optimization (filters, kVp, PCD thresholds)

Iron has **no usable K-edge** (7.1 keV, far below the diagnostic window), so all
contrast is **photoelectric** and lives at low energy. Computed with
`src/spectral.py` on the **real CONRAD standard spectrum (90 kVp)** (matched-filter
/ optimal energy-weighting framework, fixed 70 000 air photons/pixel, realistic
6 mg dose through ~10 cm tissue):

- **Optimal monochromatic energy ≈ 30 keV** (lower = more contrast, but the body
  extinguishes flux below ~30 keV).
- **kVp:** lower is better — 60 kVp gives **1.34×**, 80 kVp **1.09×** the ideal
  CNR of the 90 kVp standard; 120 kVp is **0.85×** (harder = worse).
- **Filters:** hardening filters HURT iron contrast — Cu 0.3 mm → 0.82×,
  Sn 0.5 mm → **0.58×**. (They are for spectral *separation*, not low-Z contrast.)
  Extra Al softening gives little (0.96×); the win is keeping low energy, not
  adding a K-edge filter.
- **Detector weighting still helps:** on the softer 90 kVp standard, EID already
  captures **72%** of ideal CNR, so optimal photon-counting weighting gives
  **≈1.35× the EID CNR** (vs ~2× on a harder 120 kVp beam).
- **Optimal PCD thresholds (real 90 kVp spectrum, through rabbit):** 2-bin split
  at **47.5 keV** (1.31× vs EID, 94% of ideal); 3-bin at **37.5 / 50 keV**
  (1.35× vs EID, 97%). Thresholds isolate the photoelectric-rich low band from
  the Compton high band; **not** at a K-edge.

Added study dimension: sweep `KVP_LEVELS` and `FILTER_CONFIGS` (config.py) and
report EID vs PCD-unweighted vs PCD-optimal-weighted CNR.

### 5.9 Tumor distribution factor — two experiments (cellular vs vascular)

The **tumor-distribution factor** (§5.6) encodes two biological phases of the same
delivered iron, run as one directly-comparable factor:

- **Study A — homogeneous (cellular uptake).** Iron internalised by the tumor cells
  gives a ~uniform tumor iron distribution (the uniform level). Concentrations are
  sampled from the article's **measured cellular loading** (Fig 5, all configurations
  incl. fresh; 0 h / 24 h): SPION I-113 nm **8.23 / 3.78**, SPION II-115
  **3.86 / 1.52**, II-98 **3.60 / 1.37**, II-76 **2.51 / 0.85** pg Fe/cell. Each
  configuration converts to a tumor mg Fe/ml via the tumor **cell density** and uses
  **its formulation's coating** (SPION I φ ≈ 0.15, SPION II φ ≈ 0.36; §5.2). Cite
  Heinen et al.
- **Study B — vascular / "fresh" delivery.** Freshly-injected SPIONs still in the
  **blood carrier inside the vessels, before cellular uptake** — contrast confined to
  **150 µm vessels occupying 10 % of the tumor volume** (the vessel level),
  heterogeneous, at the **injection concentration**. Referenced to the Genç/Lyer
  flow-accumulation context (suspension ~0.84 mg Fe/mL).

For Study B the same delivered iron mass is confined to 10 % of the volume →
**10× local vessel concentration**.

> **Two open decisions (TODO — pending the user):**
> (a) **Study A** tumor **cell density** used to convert pg Fe/cell → mg Fe/ml;
> (b) **Study B** **injection concentration** (article suspension 1–10 mg Fe/ml, or a
> specific injected dose). Both remain placeholders until the user fixes them.

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
detectability analysis**. Uniform vs vessel is a **factor of the single study**
(§5.6) — not a separate run — so the two distributions are directly comparable
within one design at every detector / bone / dose / concentration cell.

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

**Resolved** (user, 2026-07-08, updated 2026-07-09): tumor = 8 cm³; phantom =
rabbit-scale **ICRU soft-tissue** body + cortical-**bone** insert; FOV = 20 cm;
detectors = **both** EID and multi-bin PCD; detectability = **both** CNR and ΔHU;
concentrations = 0/0.5/1/2/5/10/20 mg Fe/ml.
**Design updates (2026-07-09):** water beam-hardening precorrection = **always on**
(was off/on); **short scan = always on**; **dose** = low 20 000 / high 100 000 ph/px
**factor** (was fixed 70 000); tumor **distribution** = uniform / vessel **factor**;
**bone** = absent / present **factor**; spectrum fixed at **90 kVp** (kVp/filtration
is the separate §5.8 optimization). → **Sufficient to implement.**

**Defaults recorded at implementation time** (non-blocking, will be read from the
installed CONRAD and logged): standard-spectrum kVp/filtration & energy sampling
(M3), C-arm SID/SDD/detector/angular range (M4), PCD bin thresholds (from the
spectrum range), reconstruction voxel size, iron-oxide polymorph fixed to Fe₃O₄
for reporting.
