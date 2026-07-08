# Development Log — SPIONvsXRay

Reverse-chronological log of progress. Newest entries on top.

---

## 2026-07-08 — SPIE Medical Imaging paper plan

- Downloaded the **official SPIE Proceedings LaTeX template** (spie.cls v3.4,
  2015/08/14) from the arXiv source of the SPIE style guide; kept spie.cls,
  spiebib.bst, the guide (`spie_template.tex`) and its two demo figures under
  `paper/template/`. **Verified it compiles** locally (`latexmk -pdf`, exit 0).
- Wrote `paper_plan.md`: target = SPIE Medical Imaging, *Physics of Medical
  Imaging* track (~8 pp.); working titles; full section-by-section structure;
  figure/table plan; result-dependent-claim placeholders; template usage notes.
  Flagged TODOs: confirm SPIE MI 2027 dates + track, confirm co-authors (SEON
  nanoparticle group).
- Added a **writing-style analysis of Andreas Maier** (from *A Gentle
  Introduction to Deep Learning in Medical Image Processing*, arXiv:1810.05401):
  11 observed traits with quoted evidence, distilled into concrete "Maier voice"
  rules for drafting this paper, plus anti-patterns to avoid.

**Next:** run the simulation (M0→M7) to produce the figures/tables, then draft
`paper/spie_manuscript.tex` from the template applying the style rules.

---

## 2026-07-08 — Dose model: delivered mass, anchored at 6 mg / 10 mg-ml

- **Decision (user):** the dose is a *delivered SPION mass*, not a fixed vial
  concentration. Anchor: **10 mg/ml formulation → 6 mg SPIONs** spread over the
  8 cm³ tumor; other formulation levels scale the delivered mass proportionally.
  "6 mg SPIONs" = whole-particle mass; magnetite core → 72.4% Fe.
- **Conversion:** `c_Fe = 0.724 × (6·c_form/10) / 8 = 0.0543 × c_form` mg Fe/ml.
  Tumor iron per level: c_form 0.5/1/2/5/10/20 → Fe 0.027/0.054/0.109/0.271/
  0.543/1.086 mg/ml.
- **Implication:** even the top dose is only ~0.5 mg Fe/ml — ~10× below iodine
  CT enhancement. Detectability at these low iron loads is the central question;
  expect SPIONs to be borderline/invisible except at the high end. Updated
  SPEC §5.1/§5.2/§5.6 and README.

---

## 2026-07-08 — Material model corrected: iron in tissue, not water

- Caught an inconsistency: modeling the tumor as *water + Fe* would give a
  spurious ~40 HU water-vs-tissue offset at c = 0, confounding ΔHU/CNR.
- **Decision (user):** model the tumor as **iron-loaded ICRU soft tissue** (same
  matrix as background + Fe @ c mg/ml). The contrast is **not diluted in water**.
  → c = 0 ≡ background (ΔHU = 0); all contrast attributable to iron; tumor-free
  soft tissue is a true iron-free reference. Updated SPEC §5.2 and README.

---

## 2026-07-08 — Experimental design finalized

**Decided (user input)**
- Dose: **70 000 photons/pixel** (I₀) → Poisson noise at projection level.
- Tumor: **8 cm³** sphere (⌀ ≈ 2.48 cm), single SPION insert.
- Phantom: **rabbit-scale** soft-tissue cylinder (~10–12 cm) **+ bone insert**
  (bone chosen so the beam-hardening on/off factor is actually meaningful).
- **20 cm FOV**, 500 projections, standard C-arm geometry.
- Detectors: **both** EID and **energy-resolved multi-bin PCD** (default 3 bins).
- Beam-hardening correction: **both** off and on.
- Detectability: **both** ΔHU and CNR (Rose CNR ≥ 3–5).
- Concentrations: **0, 0.5, 1, 2, 5, 10, 20 mg Fe/ml** (7 levels).

**Design → run count**
- Factorial: 7 conc × 2 detectors × 2 BH × 10 noise reps.
- **14** polychromatic forward projections (500 views) → 140 noisy sinograms →
  **280** reconstructions + **28** noise-free references ≈ **308** analyzed volumes.

**Status:** all information needed to implement is now in hand; remaining items
are CONRAD-default values to be read and logged during M3/M4 (spectrum kVp,
SID/SDD, detector geometry, PCD bin thresholds). Updated `SPEC.md` §5/§8 and
`README.md`.

**Next:** start **M0** — Python 3.9/3.10 venv + `pip install pyconrad` +
`setup_pyconrad()` smoke test on this Mac.

---

## 2026-07-08 — Project kickoff & planning

**Done**
- Read the reference article (Heinen et al., *Ultrason. Sonochem.* 130 (2026)
  107876): PAA-coated SPIONs, iron-oxide cores ~8/12 nm, used at ~1–10 mg Fe/ml,
  cellular loading ~8 pg Fe/cell. Established that CT contrast is driven by iron
  content → model suspensions as *water + Fe @ c mg/ml*.
- Reviewed CONRAD tutorials: custom materials (WeightedAtomicComposition /
  Mixture / MaterialsDB), spectral absorption (PolychromaticXRaySpectrum +
  PolychromaticAbsorptionModel). Noted the reference impls are Java classes in
  `edu.stanford.rsl.tutorial.physics.*`, called from Python via pyconrad.
- Checked local environment:
  - Java **1.8.0_291** ✓ (CONRAD needs Java 8).
  - Python **3.14.4** ⚠ (too new for pyconrad/JPype — will use a 3.9/3.10 venv).
  - Platform macOS/Apple Silicon ⚠ (pyconrad notes Mac caveats).
  - pyconrad **0.8.0** on PyPI, not yet installed.
- Wrote `SPEC.md` (full plan, milestones M0–M7), `README.md` (overview +
  missing-information list), this `DEVLOG.md`, and `.gitignore` (ignores the
  copyrighted PDF and regenerable outputs).
- Confirmed git remote `git@github.com:akmaier/SPIONvsXRay.git`, branch `main`.

**Decisions**
- Experiment = multi-vial contrast-detectability phantom (water cylinder + rod
  inserts, one iron concentration each), 500-projection C-arm cone-beam scan
  with CONRAD standard spectrum, FDK reconstruction, then HU/CNR analysis.
- First pass noise-free (contrast ceiling); second pass adds Poisson noise.

**Open questions (see README §"Missing information")**
- Iron-oxide polymorph/density; exact concentration list & units; standard
  spectrum kVp & dose/noise model; phantom & C-arm geometry parameters;
  detectability criterion.

**Next**
- M0: create Python 3.9/3.10 venv, `pip install pyconrad`, verify
  `setup_pyconrad()` starts the JVM and downloads CONRAD.jar on this Mac.
- Verify the Python-side pyconrad API (ClassGetter, grid/NumPy interop) and that
  the tutorial physics classes are present in the shipped CONRAD.jar.
