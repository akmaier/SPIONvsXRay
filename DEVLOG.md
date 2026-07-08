# Development Log — SPIONvsXRay

Reverse-chronological log of progress. Newest entries on top.

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
