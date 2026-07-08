# Development Log — SPIONvsXRay

Reverse-chronological log of progress. Newest entries on top.

---

## 2026-07-08 — Two phantom backgrounds: round geometric + ROBY

- Decision (user): use **both** a round geometric phantom (cylinder + bone) and
  the **ROBY** digital rabbit (Segars/Duke XCAT family) for anatomical realism.
- ROBY is **licensed** (Duke/XCAT) → cannot auto-download. Plan: implement the
  round phantom fully now; add a ROBY adapter (`src/phantom.py`) that ingests a
  user-supplied ROBY voxel volume (organ→material map) and embeds the tumor.
  Batch runs on `round` until ROBY files are provided under `data/roby/`.
- config.py: `PHANTOM_BACKGROUNDS=['round','roby']`, `TUMOR_MODELS=
  ['homogeneous','vessel']`, `ROBY` paths, vessel params. SPEC §5.1 updated.

**ACTION FOR USER:** to enable the ROBY arm, drop the ROBY-generated volume +
log into `data/roby/` (they're licensed; I can't fetch them).

---

## 2026-07-08 — Correctness: model iron OXIDE (magnetite), not pure Fe

- Caught: attenuation used pure elemental Fe. Fixed to magnetite Fe₃O₄ — added
  the bound oxygen (0.382 g O per g Fe) via `oxide_contrast_massatten(E) =
  (μ/ρ)_Fe + 0.382·(μ/ρ)_O` in `materials.py`, propagated to `spectral.py`.
- Effect: realistic-dose tumor HU **2.85 → 3.04** (+7%); optimal PCD thresholds
  unchanged (35/47.5 keV), gains 1.83× vs EID. Dashboard figures refreshed.
- Open: soft-tissue matrix still water proxy (CONRAD lacks ICRU-44 XML); could
  build ICRU-44 by composition later. Rabbit background still geometric (see below).

---

## 2026-07-08 — Spectral optimization + Study B (vessel phantom) designed

**Contrast optimization (`src/spectral.py`, matched-filter framework):**
- Iron has no usable K-edge (7.1 keV) → contrast is photoelectric, lives low.
  Optimal mono energy ≈ **30 keV**.
- kVp: 60 kVp → **1.26×**, 80 kVp → 1.15× the 120 kVp baseline ideal CNR.
- Hardening filters HURT: Cu 0.3 mm → 0.70×, **Sn 0.5 mm → 0.42×**.
- **Detector weighting dominates:** EID = 51% of ideal CNR; optimal PCD weighting
  ≈ **2× EID** (≈4× dose-equivalent).
- **Optimal PCD thresholds:** 2-bin @ **40 keV** (1.78× vs EID); 3-bin @
  **35 / 47.5 keV** (1.87× vs EID, 95% of ideal). Updated config.py defaults.
- Figures (`spectrum.png`, `pcd_bins.png`) added to `docs/assets/` → dashboard
  Spectrum panel now shows real data. Added kVp/filter sweep to the study.

**Study B — vascular tumor phantom (designed, SPEC §5.9):** contrast confined to
150 µm vessels at 10 % tumor volume (10× local conc, mass conserved). Vessels are
sub-resolution (390 µm voxels) → mean ΔHU ≈ Study A; the study isolates
second-order effects (beam-hardening nonlinearity, structural noise, partial
volume). Runs the full factorial × filters/detectors through the M4–M6 pipeline
once built.

---

## 2026-07-08 — M1 COMPLETE ✅ (materials + first real result)

- Introspected the CONRAD material API from the live JVM: `MaterialsDB.getMaterial`
  serves elements + `water`/`bone`/`iron`; `Material.getAttenuation(E_keV, TYPE)`.
- **Units gotcha caught & fixed:** `getAttenuation` returns **linear** μ [1/cm],
  not mass. Verified: water@60keV=0.206 (ρ=1), iron@60keV=9.48=1.205×7.87. So
  `mass_attenuation = getAttenuation / density`. After the fix iron μ/ρ@60keV =
  **1.2050** = NIST 1.205 exactly.
- `src/materials.py`: soft tissue = water proxy (no ICRU-44 XML in CONRAD);
  iron-loaded tumor via the mixture rule `μ_soft + c_Fe[g/cm³]·(μ/ρ)_Fe`.
- **Key finding:** at the realistic 6 mg dose (0.543 mg Fe/ml) the tumor is only
  **~2.9 HU** above background (≈5.7 HU at 2× dose) — quantitatively at the CT
  detection floor. Figures (`mu_vs_energy.png`, `hu_vs_conc.png`) copied into
  `docs/assets/` — the dashboard's Materials panels now show real data.

**Next (M2):** rabbit phantom (soft tissue + bone + 8 cm³ tumor); then M3 spectrum.

---

## 2026-07-08 — M0 COMPLETE ✅ (pyconrad bridge working)

Working environment on Apple Silicon nailed down after three blockers:
1. **Arch mismatch** — arm64 Python can't load the system x86_64 `libjvm.dylib`.
   → installed a native **arm64 Zulu JDK 8** into `.jdk8/` (`scripts/install_jdk8.sh`).
2. **JPype too new** — JPype1 1.7.1 requires Java 9+, but CONRAD needs Java 8.
   → pinned **JPype1==1.5.0** (last Java-8-compatible release).
3. **Env wiring** — `src/conrad_backend.py` auto-sets `JAVA_HOME` to the bundled
   JDK and starts the JVM idempotently.

Result: `pyconrad.setup_pyconrad()` succeeds and a CONRAD `Grid2D` round-trips
through the JVM (`CONRAD bridge OK: True`) with **no manual environment setup**.
Stack: Python 3.12 + pyconrad 0.8.0 + JPype1 1.5.0 + arm64 Zulu JDK 8.

Pinned `JPype1==1.5.0` in `requirements.txt`; updated README getting-started.

**Next (M1):** custom materials (soft tissue, cortical bone, iron-loaded tumor)
and a μ(E)/HU sanity check.

---

## 2026-07-08 — Implementation start (M0) + 30-min loop

- **Paper metadata confirmed:** SPIE MI track = *Physics of Medical Imaging*;
  co-authors = Maier, Stefan Lyer, Lukas Heinen, Rainer Tietze (SEON /
  Universitätsklinikum Erlangen). Updated `paper_plan.md`.
- **Environment probe (M0):** available Pythons = 3.12/3.13/3.14 (arm64) + old
  3.6 (x86_64); Java 8 present but **x86_64** (Oracle). ⚠ Apple-Silicon arch
  mismatch risk: JPype needs Python and JVM at the same arch. Testing pyconrad on
  arm64 py3.12 first; arm64 Temurin/Zulu 8 is the fallback if the x86_64 JVM
  can't be bridged.
- **Scaffolding (backend-agnostic):** `requirements.txt`, `src/` package,
  `src/config.py` = single source of truth for every SPEC parameter. Verified
  `python src/config.py` reproduces the dose table (c_form→mg Fe/ml) and the
  7×2×2×10 = 280 factorial.
- Started a **30-min work/report loop**; pyconrad install running in background.

**Loop cadence:** each iteration advances a milestone, commits + pushes, reports
to chat. Next: resolve the pyconrad/JVM arch question (M0), then M1 materials.

---

## 2026-07-08 — MIT license + GitHub Pages audit dashboard

- Added **MIT `LICENSE`** (© 2026 Andreas Maier); README + paper plan updated.
- Built a **GitHub Pages audit dashboard** under `docs/` (`index.html`,
  `assets/`, `.nojekyll`): sections for spectrum, materials/attenuation, phantom,
  projections, reconstructions, detectability, an annotated **code-audit**
  section, and reproducibility. Self-contained (inline CSS/JS); figure panels
  show a "pending" placeholder until `docs/assets/*.png` exist, then auto-render.
- Planned `src/build_dashboard.py` to regenerate `docs/` from `results/`
  (added to SPEC deliverables + new milestone **M8**).
- Enabled GitHub Pages (source = `main`/`docs/`) via `gh api`.
  Live at <https://akmaier.github.io/SPIONvsXRay/>.

**Next:** unchanged — run M0→M7, then M8 fills the dashboard with real figures.

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
