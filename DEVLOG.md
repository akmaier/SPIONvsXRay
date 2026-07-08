# Development Log — SPIONvsXRay

Reverse-chronological log of progress. Newest entries on top.

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
