# SPIONvsXRay

**Can superparamagnetic iron-oxide nanoparticles (SPIONs) be seen in CT?**

A fully reproducible, **CONRAD-native** simulation study of the X-ray absorption
and cone-beam CT reconstruction of iron-oxide nanoparticles at realistic
biological iron loads — to determine whether, and at what concentration, they
become visible in a C-arm CT, and how detector type, beam-hardening correction,
and spectral shaping shift that limit. Built on [CONRAD](https://www5.cs.fau.de/conrad/)
via [`pyconrad`](https://pypi.org/project/pyconrad/).

The nanoparticles are those in Heinen et al., *Ultrasonics Sonochemistry* 130
(2026) 107876: **PAA-coated magnetite (Fe₃O₄) SPION clusters**, cores ~8/12 nm,
hydrodynamic clusters ~70–250 nm, used at ~1–10 mg Fe/ml.
*The article PDF is git-ignored (copyright).*

**Full-nanoparticle model.** We simulate the *whole particle* — magnetite core
**and** PAA coating (polyacrylic acid, monomer C₃H₄O₂) — not iron alone. The
independent variable is the tumor **iron** concentration `c_Fe` (mg Fe/ml); from
it we estimate the **whole-nanoparticle** concentration
`c_NP = c_Fe / (0.724·(1−φ))`, where 0.724 is the Fe mass fraction of magnetite
and **φ is the PAA coating mass fraction**. **φ is now taken per formulation from
the article's supplementary TGA** (Heinen et al., supplement Table A.1 / Fig A.1),
not a single estimate. The TGA inorganic (iron-oxide) residual is **87.7 % for
SPION I** (12 nm cores) and **66.4 % for SPION II** (8 nm cores); on heating the
residual oxidizes magnetite → Fe₂O₃ (+3.4 % mass — identical for the maghemite and hematite polymorphs, article Eq A.1), so the
original magnetite fraction is `R_m/1.034` = **84.8 % (SPION I)** and
**64.2 % (SPION II)**. The complement is PAA: **φ ≈ 0.15 for SPION I** and
**φ ≈ 0.36 for SPION II** — so SPION II carries roughly twice the coating. (The
earlier single central estimate φ = 0.15 turned out ≈ exact for SPION I.) Each
Study A configuration therefore carries **its formulation's coating** — SPION I
gives `c_NP ≈ 1.625·c_Fe`, SPION II `c_NP ≈ 2.16·c_Fe`. *(The current pipeline
still applies a single φ = 0.15 — exact for SPION I; per-formulation φ enters with
the Study A per-configuration sampling, in progress — see [`HANDOFF.md`](HANDOFF.md).)*
φ sets **only** the reported
mg SPION/ml; the PAA is low-Z (C/H/O), tissue/water-equivalent, so it is
**negligible for the X-ray μ** — adding it shifts the monochromatic iron ΔHU at
the top dose by only ≈ 3.6 % at φ = 0.15 (3.62 → 3.75 HU @ 60 keV). *Iron*
(mg Fe/ml) and *particle* (mg SPION/ml) concentrations are kept distinct
throughout; all contrast numbers are reported as **mg Fe/ml in the tumor**.

---

## Study design

Two phantoms, one comparable C-arm geometry:

- **Disc phantom study (2D fan-beam).** An **ED-phantom-style** calibration
  phantom: a soft-tissue cylinder with the 7 concentration inserts arranged on a
  ring at equal radius (so beam-hardening cupping is common-mode) plus a
  cortical-bone insert. Each insert is a 2.5 cm (≈8 cm³-equivalent) disk of
  **nanoparticle-loaded soft tissue** (magnetite Fe₃O₄ core **+ PAA coating**
  added to the tissue matrix, so a zero-iron insert equals background and all
  contrast is iron).
- **Rabbit case (3D, planned).** The real **RabbitCT** rabbit volume (see below)
  as anatomy, imaged with the real RabbitCT C-arm geometry — the disc study uses
  a 2D slice of the same acquisition for direct comparability.

**Dose model** — a *delivered mass*, not a fixed vial concentration: anchored at
**6 mg SPIONs for the 10 mg/ml formulation**, spread over an 8 cm³ tumor, scaled
with concentration. Magnetite (72.4 % Fe) ⇒ `c_Fe = 0.0543 × c_form`, so the top
dose is only **~0.5 mg Fe/ml** (~10× below iodine CT enhancement). This tumor band
is **grounded in the article's measured cellular loading** (Fig 5, B16-F10 melanoma,
0 h / 24 h): SPION I-113 nm **8.23 / 3.78**, SPION II-115 **3.86 / 1.52**, II-98
**3.60 / 1.37**, II-76 **2.51 / 0.85** pg Fe/cell; at a tumor cell density
~10⁸–10⁹ cells/cm³ this is ~0.25–8 mg Fe/ml (realistic ~1–2.5 mg Fe/ml for SPION I)
— bracketing the sweep (Heinen et al.).

**Two experiments (one comparable design).** The delivered iron is imaged under
two biological distributions, run as a single directly-comparable factor:
- **Study A — homogeneous (cellular uptake).** Iron internalised by tumor cells ⇒
  ~uniform tumor iron distribution. Concentrations are sampled from the article's
  **measured cellular loading** (Fig 5, all configurations incl. fresh): SPION I-113 nm
  **8.23 / 3.78** pg Fe/cell (0 h / 24 h), SPION II-115 **3.86 / 1.52**, II-98
  **3.60 / 1.37**, II-76 **2.51 / 0.85**. Each configuration converts to a tumor
  mg Fe/ml via the tumor cell density and uses **its formulation's coating**
  (SPION I φ ≈ 0.15, SPION II φ ≈ 0.36). *(Heinen et al.)*
- **Study B — vascular / "fresh" delivery.** Freshly-injected SPIONs still in the
  **blood carrier inside the vessels**, before cellular uptake ⇒ contrast confined
  to the 150 µm vessels (~10 % of tumor volume, 10× local conc), heterogeneous, at
  the injection concentration. Referenced to the Genç/Lyer flow-accumulation context
  (suspension ~0.84 mg Fe/mL).

> **Two open decisions (TODO — pending the user):**
> (a) **Study A** tumor **cell density** used to convert pg Fe/cell → mg Fe/ml;
> (b) **Study B** **injection concentration** (article suspension 1–10 mg Fe/ml, or
> a specific injected dose). Both are placeholders until fixed.

**Detectability study** (a multi-factor effects study, *not* a mathematical
factorial): formulation conc. {0, 0.5, 1, 2, 5, 10, 20 mg/ml} × detector {EID,
multi-bin PCD} × bone {absent, present} × dose {low 20 000, high 100 000
photons/px} × distribution {homogeneous, vessel}, with **30 noise realizations**
per cell. Water beam-hardening precorrection and the minimal short scan are
**always on** (not factors). Detectability reported as **ΔHU and CNR** (Rose ≥ 3–5).
See [`SPEC.md`](SPEC.md) for full parameters.

## Key findings so far

- **SPIONs sit right at the CT detection limit at the realistic dose.** Iron
  tumor contrast is ~3.4 HU (EID) / ~4.7 HU (PCD) at the realistic 6 mg dose and
  ~7–9 HU at 2× dose. Detectability crosses the Rose criterion (CNR ≈ 5) *exactly*
  at the realistic dose (c_Fe ≈ 0.54 mg Fe/ml: CNR 5.4 EID / 5.1 PCD @ 70 000
  photons/px, 0.5 mm voxels) — the 6 mg SPION load is a borderline, not a
  comfortable, detection. *(GPU 0.5 mm detectability run, 30 noise realizations/cell,
  70 000 photons reference dose; `results/factorial/`. These headline numbers predate
  the dose/distribution/per-formulation redesign — pending regeneration.)*
- **Neither photon-counting nor beam-hardening correction lowers the threshold.**
  All four cells (EID/PCD × BH off/on) share the same detection limit
  (~0.54 mg Fe/ml). PCD's matched-filter energy weighting delivers its ideal
  ~1.3× CNR gain only at 2× dose (c_Fe 1.09: PCD CNR 11.5 vs EID 8.8 = 1.31×);
  at the realistic limit the gain is swamped by noise. BH correction changes iron
  detectability negligibly — iron contrast is not a beam-hardening artifact.
- **Iron has no usable K-edge** (7.1 keV) → contrast is photoelectric and lives
  at low energy. Optimal mono energy ≈ 30 keV; **lower kVp helps** (60 kVp =
  1.34× the 90 kVp ideal CNR), **hardening filters hurt** (Sn → 0.58×); the
  ideal-observer PCD ceiling is ≈ 1.35× EID (optimal 3-bin thresholds ~37.5/50 keV).
- The CONRAD magnetite material and an independent NIST-based model agree to 5
  decimals — the physics is cross-validated.

## Pipeline (all CONRAD-native)

> **Policy: all reconstruction uses the CONRAD API — no hand-rolled recon math.**
> Projection (`FanBeamProjector2D`), redundancy/short-scan weighting
> (`redundancy.ParkerWeights` / `SilverWeights` / …), filtering
> (`CosineFilter`, `RamLakKernel` and its roll-off variants), and backprojection
> (`FanBeamBackprojector2D`) are all CONRAD classes, used in CONRAD's geometry
> convention (e.g. `focalLength` derived from `maxT`, `deltaT`, and the fan angle —
> **no ad-hoc magnification or padding factors**). We do not mix in custom ramp
> filters, custom Parker weights, or custom geometry. Only the polychromatic /
> spectral bookkeeping (per-energy attenuation, EID/PCD detector combination, noise)
> lives outside CONRAD's fan classes.


| Module | Role |
|--------|------|
| `src/config.py` | single source of truth for all parameters |
| `src/conrad_backend.py` | JVM/CONRAD bootstrap (JAVA_HOME, OpenCL, extension classpath) |
| `src/materials.py` | attenuation via CONRAD's material DB; magnetite oxide model |
| `src/spectrum.py` | real CONRAD 90 kVp polychromatic spectrum (+ kVp/filter variants) |
| `src/spectral.py` | spectral optimization (optimal energy weighting, PCD thresholds) |
| `src/conrad_phantom.py` | CONRAD `AnalyticPhantom` (ED phantom + registered full-particle SPION materials: magnetite core + PAA coating) |
| `src/conrad_project.py` | `PriorityRayTracer` base-material sinograms + polychromatic EID/PCD + noise |
| `src/conrad_ct.py` | CONRAD fan-beam projection + FBP (GPU when available) + calibrated water precorrection |
| `src/run_detectability.py` | the detectability study: per-insert ΔHU + CNR + detection thresholds |

## Spectral processing (EID vs. PCD)

Both detectors start from the same polychromatic forward model. For each energy
`E`, the detected photons through the object are `N(E) = N0·s(E)·exp(−τ)` with `s`
the normalized 90 kVp spectrum and `τ = Σ_material μ_material(E)·L_material` (analytic
base-material path lengths × CONRAD's energy-dependent `μ(E)`). The detectors differ
in how they collapse `N(E)` to a line integral and how beam hardening is corrected.

**Energy-integrating detector (EID).** Signal is the energy-weighted sum
`S = Σ_E N(E)·E`. Quantum noise is Poisson per energy, so `S` has variance
`Σ_E N(E)·E²` (the exact second moment of the energy-weighted Poisson sum — *not*
`Poisson(S)`; we add Gaussian noise with that variance). Line integral `p = −log(S/S_air)`.
Beam hardening: a **single** water precorrection polynomial calibrated for the `s·E`
spectrum (`conrad_ct.water_precorrection_poly`) linearizes the water response so a
water cylinder reconstructs flat.

**Photon-counting detector (PCD).** Photons are sorted into 3 energy bins
(thresholds 37.5 / 50 keV). Per bin `b`: Poisson counts `C_b = Σ_{E∈bin} N(E)`. Bins
are combined with **matched-filter weights** `w_b ∝ S_b/V_b`
(`S_b = Σ_bin N_t·c`, `V_b = Σ_bin N_t`; `N_t` = transmitted background photons,
`c` = per-energy iron contrast) — the CNR-optimal weighting that reaches the
`Σ S_b²/V_b` detectability ceiling (`src/spectral.py`). Combination is in the
**count domain** — `M = Σ_b w_b·C_b`, then a single `−log(M/M_air)` — which is robust
to the photon-starved low-energy bin (a per-bin `−log` diverges when its counts hit
zero). Beam hardening is **energy-dependent / per-bin**: each bin is corrected on
*its own* spectrum with a bin-specific water precorrection polynomial, applied as a
**corrected count** `C_b_corr = C_air_b·exp(−poly_b(p_b))` before the count-domain
combination. This linearizes each narrower bin's hardening on its own spectrum while
keeping the count-domain noise robustness (a pure per-bin *log-domain* combination is
equivalent for water cupping but noisier, so it is not used in the noisy study).

The uncorrected water precorrection (a fixed `p + 0.10·p²`) is deprecated — it was
~7× too strong and over-corrected cupping into capping; the calibrated polynomial
flattens a monochromatic-equivalent water disk to ~0.1%. The Poisson/energy noise
model is verified against CONRAD's `StatisticsUtil` RNG and `PolychromaticAbsorptionModel`,
and GPU versions (Philox RNG + OpenCL EID/PCD spectral detectors) are upstreamed to CONRAD.

## Getting started

CONRAD needs **Java 8**. Verified environment (Apple Silicon M1):
Python 3.12 + pyconrad 0.8.0 + **JPype1 1.5.0** (newer JPype dropped Java 8) +
a native **arm64 Zulu JDK 8** (bundled locally).

```bash
python3.12 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt          # pyconrad + JPype1==1.5.0 + numerics
bash scripts/install_jdk8.sh             # arm64 Zulu JDK 8 -> ./.jdk8 (no sudo)
bash scripts/install_opencl.sh           # optional: OpenCL (jogamp 2.6.0 + shim)
bash scripts/build_conrad_ext.sh         # optional: build the fixed GPU backprojector
python src/conrad_backend.py             # smoke test: JVM up, CONRAD Grid2D round-trips
```

`conrad_backend` auto-configures `JAVA_HOME`, the OpenCL shim, and the extension
classpath; the pipeline falls back to CPU if OpenCL is absent.

## CONRAD extension — fixed OpenCL fan backprojector (`conrad_ext/`)

CONRAD's OpenCL fan backprojector shipped incomplete: its kernel was never passed
the reconstruction pixel spacing (`// TODO: Spacing` in the source) plus some
`//FIXME` detector-geometry hacks, so it couldn't do configurable voxel sizes and
corrupted quantitative reconstructions. We fixed it **in the original class**
`edu.stanford.rsl.tutorial.fan.FanBeamBackprojector2D` (+ its `.cl` kernel): added
`setSpacing`, passed the spacing to the kernel, and corrected the 2D detector
geometry to match the CPU path. Result: correct **0.5 mm voxels** and **~850×
faster** GPU reconstruction (matches CPU to ~1e-5). Source is in `conrad_ext/`
(git-tracked, contributable upstream); `scripts/build_conrad_ext.sh` compiles it
against the bundled jar and it is classpath-shadowed at runtime.

## RabbitCT dataset — realistic rabbit anatomy & geometry (published)

The realistic rabbit is the **RabbitCT** benchmark (a post-mortem C-arm cone-beam
scan from the FAU Pattern Recognition Lab). The original challenge site is offline;
with the authors' agreement we **re-published the dataset on Zenodo** under
**CC BY 4.0**:

➡ **DOI: [10.5281/zenodo.21267885](https://doi.org/10.5281/zenodo.21267885)**
&nbsp;·&nbsp; record: <https://zenodo.org/records/21267885>
&nbsp;·&nbsp; concept DOI: 10.5281/zenodo.21267884

Please cite the original paper: Rohkohl C, Keck B, Hofmann HG, Hornegger J.
*Technical Note: RabbitCT — an open platform for benchmarking 3D cone-beam
reconstruction algorithms.* **Medical Physics** 36(9):3940–3944, 2009.
DOI: [10.1118/1.3180956](https://doi.org/10.1118/1.3180956). The scan is
**post-mortem** (no live-animal procedure).

Re-publication tooling is in [`publish/rabbitct/`](publish/rabbitct/) (data
descriptor, Zenodo metadata, and a python3 uploader run from the LME server). The
data itself is git-ignored (`data/rabbitct/`, fetched from LME via SSH).

## Audit dashboard

A public **GitHub Pages dashboard** shows all intermediate results — spectra,
material attenuation, phantom, projections, reconstructions, detectability curves
— with **annotated code snippets for auditing**, regenerated from `results/`.

➡ **https://akmaier.github.io/SPIONvsXRay/** (source: `docs/`)

## Repository layout

```
SPEC.md            Full experimental plan and milestones
HANDOFF.md         Continuation guide for the next contributor (start here)
DEVLOG.md          Running development log (reverse-chronological)
paper_plan.md      SPIE Medical Imaging paper plan + Maier writing-style analysis
paper/template/    Official SPIE Proceedings LaTeX template (spie.cls v3.4)
src/               CONRAD-native simulation & analysis pipeline
conrad_ext/        Patched CONRAD class (fixed OpenCL fan backprojector) + build
scripts/           Environment setup (JDK 8, OpenCL, conrad_ext build)
publish/rabbitct/  RabbitCT Zenodo re-publication package
docs/              GitHub Pages audit dashboard (index.html + assets/)
data/, results/    Git-ignored: fetched data, generated outputs
.jdk8/, .jogamp/   Git-ignored: bundled arm64 JDK 8 + OpenCL shim
```

## License / attribution

Project code: **MIT License** ([`LICENSE`](LICENSE)). Simulation built on CONRAD
(FAU Pattern Recognition Lab). RabbitCT data: **CC BY 4.0** (cite the Med Phys
paper above). The reference SPION article PDF is not redistributed (git-ignored).
