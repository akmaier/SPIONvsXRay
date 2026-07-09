# HANDOFF — SPIONvsXRay

For the next agent picking up this project. Read this, then `README.md` (state) and
`DEVLOG.md` (chronological detail). `SPEC.md` = the experimental plan. Commit and
push after each milestone; the user works on `main` (`git@github.com:akmaier/SPIONvsXRay.git`).

## 1. What this is
Simulation study: **can SPIONs (magnetite iron-oxide nanoparticles) be seen in
CT** at realistic biological iron loads? Fully **CONRAD-native** via `pyconrad`.
Aiming at a **SPIE Medical Imaging** paper (Physics of Medical Imaging track;
authors Maier, Lyer, Heinen, Tietze — see `paper_plan.md`).

## 2. Current state — WORKING end-to-end
The CONRAD-native 2D disc-phantom pipeline runs at **0.5 mm on the GPU**:
`conrad_phantom` (ED phantom, `AnalyticPhantom` + magnetite materials) →
`conrad_project` (`PriorityRayTracer` base-material sinograms + polychromatic
EID/PCD + Poisson noise) → `conrad_ct.fbp` (CONRAD fan FBP) →
`conrad_project.measure_inserts`. Sanity run: `python src/conrad_project.py`
gives monotonic per-insert iron ΔHU (0 → **+7.14 HU** at c_form=20), noise-free EID.

**Headline finding (GPU 0.5 mm factorial, 30 reps/cell, `results/factorial/`):**
SPIONs sit **right at the detection limit** at the realistic dose — iron ΔHU ≈ 3.4
(EID) / 4.7 (PCD), CNR ≈ 5.4/5.1 (just crosses Rose 5) at c_Fe ≈ 0.54 mg Fe/ml.
Detection threshold ≈ 0.54 mg Fe/ml is the **same for all four** EID/PCD × BH
off/on cells: neither photon-counting nor BH correction lowers it. PCD delivers
its ideal ~1.3× CNR gain only at 2× dose (c_Fe 1.09: PCD 11.5 vs EID 8.8). Iron
has no usable K-edge → low-energy/photoelectric contrast; lower kVp helps (60 kVp
1.34×), hardening filters hurt (Sn 0.58×), ideal-observer PCD ceiling ≈ 1.35× EID
(optimal 3-bin thresholds ≈ 37.5/50 keV). See `src/spectral.py`.

## 3. Environment — the hard-won parts (don't re-derive)
- macOS **Apple Silicon (M1 Max)**. ⚠ **Disk ~99 % full (~14 GB free)** — do NOT
  fetch the full 17 GB RabbitCT set.
- `.venv` = Python **3.12**; `pip install -r requirements.txt`. pyconrad 0.8.0 +
  **JPype1 must stay 1.5.0** (1.6+ drops Java 8; CONRAD needs Java 8).
- CONRAD needs Java 8 and the JVM must be **arm64** (can't load the system x86_64
  JVM). `scripts/install_jdk8.sh` puts an **arm64 Zulu JDK 8** in `.jdk8/`.
- **OpenCL (Path A):** works on the M1 GPU via **jogamp 2.6.0** + a 1-line
  `OpenCL.framework` reexport shim (`.jogamp/libOpenCL.dylib`), installed by
  `scripts/install_opencl.sh`. (Root cause of the original failure: a gluegen
  aarch64 CPU-detection guard, not the hardware.)
- **`conrad_ext/`** = our patched `FanBeamBackprojector2D` (adds `setSpacing` +
  passes spacing to the CL kernel + fixes the 2DCL detector geometry). Build with
  `scripts/build_conrad_ext.sh`; `conrad_backend` puts `conrad_ext/out` on
  `CONRAD_DEV_DIRS` so it **shadows** the jar class. This gives correct 0.5 mm +
  ~850× GPU backprojection.
  - ⚠ **This fix (plus a GPU `exp`/`log` operator fix) is now upstreamed to
    akmaier/CONRAD `master` (`71bb20c`)** and baked into a patched
    `conrad_1.1.0.jar` (`publish/conrad/`, git-ignored). Once that jar is
    re-hosted on FAU and picked up by pyconrad (see `publish/conrad/RELEASE.md`),
    `conrad_ext` is redundant — retire it by dropping the `CONRAD_DEV_DIRS` shadow
    in `conrad_backend`. **Until then, keep `conrad_ext`** (our env still runs the
    currently-hosted, unpatched jar). Rebuild the jar with
    `scripts/rebuild_conrad_jar.sh`.
- `conrad_backend.setup()` wires ALL of the above (JAVA_HOME, DYLD, CONRAD_DEV_DIRS)
  and starts the JVM idempotently. Always import + `setup()` before CONRAD calls.
- **Reference CONRAD source** cloned read-only at `~/Documents/CONRAD` (grep it).

### Gotchas
- `Material.getAttenuation(E, TYPE)` returns **linear μ [1/cm]**, not mass; mass = linear/density.
- numpy 2.5: use `np.trapezoid` (not `np.trapz`).
- pyconrad **skips `loadConfiguration()`** (avoids a Swing/ImageJ hang) → there is
  no global CONRAD Configuration by default; the tutorial fan classes take geometry explicitly.
- Recon voxel size = `geo["voxel_mm"]` (0.5), applied via the patched backprojector's
  `setSpacing`; `measure_inserts` reads `geo["voxel_mm"]` for ROI px mapping.
- Long JVM runs: launch in the background; filter stdout noise with
  `grep -viE "vtk|UserWarning|warnings.warn|stack guard|could not install|Replacing File|globalWork|CLDevice|Context:|Device:"`.

## 4. Open work (priority order)
1. ✅ **DONE — PCD per-bin bug fixed** in `src/run_factorial.py`. Root cause was
   two-fold: (a) `_pcd_weights` used contrast-only weighting that UP-weighted the
   photon-starved low-energy bin; (b) bins were combined as `Σ w_b·(−log C_b)`, so
   a per-bin log diverged when starved-bin counts hit 0 → image-wide streaks
   (CNR ≈ 0/negative). Fix: combine COUNTS then take a single log
   (`p=−log(Σ w_b C_b / Σ w_b C_air,b)`) with matched-filter weights
   `w_b ∝ S_b/V_b` (count-weighted mean iron contrast). Verified against a
   detector-level Monte Carlo (reaches the 1.29–1.34× ideal ceiling of `spectral.py`).
2. ✅ **DONE — full factorial re-run** on the GPU 0.5 mm pipeline, 30 reps/cell,
   persisted to `results/factorial/` (`factorial.csv` + `factorial.json`).
   Threshold ≈ 0.54 mg Fe/ml, same for all cells (see headline above).
   **Still TODO:** add the **filter/kVp sweep** (currently 90 kVp only) —
   `polychromatic_accumulators(base, kvp=..., filters=...)` already accepts them.
### RabbitCT .rctd format (CRACKED — see `src/rabbitct.py`)
The `.rctd` binary was undocumented in the distributed sources; reverse-engineered:
`[24-byte header: uint32 version=2, S_x=1248, S_y=960, numProj=496; float32 R_L=1.2076,
O_L=24.0]` then per projection `[96-byte 3×4 double matrix A_n, COLUMN-major] +
[S_x·S_y float32 image]`. The matrix maps world[mm]→detector[px]; decomposing gives
the real Siemens C-arm short scan: **SID 745 mm, SDD ~1196 mm (0.308 mm pitch),
1248×960 detector, 496 views over 198° (0.40°/step), rotation about z**, principal
point (558.5, 486.1) (offset detector). `src/rabbitct.py` provides
`geometry()`, `read_matrices()`, `load_reference_volume()` (256³ anatomy loads +
renders correctly). A tumor-insertion PoC works. **Next: full 3D** — build a CONRAD
`ProjectionTableTrajectory` from the 496 matrices (or set SID/SDD/angles on a
config), insert the SPION tumor into `reference_256.vol`, cone-beam forward-project
+ FDK, run the same detectability analysis. The GPU spectral detectors + RNG now in
CONRAD make the on-device version feasible.

### Older open items
3. **RabbitCT C-arm geometry** for the Conrad config: extract the 496 3×4
   projection matrices from `data/rabbitct/rabbitct_512-v2.rctd` (format in
   `data/rabbitct/develop/rabbitct_develop/include/rabbitct.h`; struct
   `RabbitCtGlobalData`, matrix `A_n`, voxel `R_L`). Use it for the rabbit case,
   and set the 2D disc geometry to match (SID/SDD/detector/angular range).
4. **Rabbit case (3D):** use `data/rabbitct/reference_256.vol` (256³ float32 rabbit
   anatomy) as the phantom, insert SPION tumor(s), project with the RabbitCT
   cone-beam geometry, reconstruct. CONRAD has RabbitCT support as `LolaBunny`
   (`edu.stanford.rsl.conrad.reconstruction.LolaBunnyBackprojector`).
5. **Dashboard build** `src/build_dashboard.py` (M8): regenerate `docs/` from
   `results/` (figures + detectability table + commit stamp). Figures are copied
   manually now. Dashboard is live at https://akmaier.github.io/SPIONvsXRay/.
6. **Paper** (`paper_plan.md`): draft `paper/spie_manuscript.tex` from
   `paper/template/` once results are in. Confirm SPIE MI 2027 dates/track.

### Cleanup
- `src/simulate.py` and `src/reconstruct.py` are the **superseded** numpy/hybrid
  pipeline (untracked). The CONRAD-native path is `conrad_project.py` +
  `conrad_ct.py` + `conrad_phantom.py`. Remove the old two.

## 5. Known limitations (state honestly in the paper)
- 2D central-slice disc study (not full 3D) for the phantom; rabbit case is the 3D arm.
- Bone-insert **streak artifacts** (~2.5 HU) contaminate low-c measurement;
  mitigated by a **local-annulus background + c0-insert subtraction** in
  `measure_inserts`. Could improve with more views / ramp apodization.
- Soft tissue = water proxy (CONRAD ships no ICRU-44 soft-tissue XML).

## 6. External resources & access
- **RabbitCT data on LME:** two-hop SSH `ssh -J maier@cluster.i5.cs.fau.de
  maier@lme31.cs.fau.de` (key-based, works). Data at
  `/disks/data1/share/webdata/fileadmin/Forschung/Software/RabbitCT`. lme31 has
  python3 (no curl/jq); cluster.i5 can't see the share.
- **RabbitCT published on Zenodo:** DOI **10.5281/zenodo.21267885** (CC BY 4.0);
  cite Rohkohl/Keck/Hofmann/Hornegger, Med Phys 36(9):3940–3944, 2009
  (10.1118/1.3180956). Upload tooling in `publish/rabbitct/`; token in the
  git-ignored `publish/rabbitct/zenodo_secrets.env`.
- **Dashboard/Pages:** `docs/` on `main`, served at the URL above.

## 7. Decisions already locked (don't relitigate)
Magnetite Fe₃O₄ (not pure Fe); delivered-mass dose (`c_Fe = 0.0543·c_form`, 6 mg
@ 10 mg/ml); iron in soft-tissue matrix (c=0 ≡ background); ED-phantom circular
inserts + bone; detectability by ΔHU **and** CNR (Rose 3–5); 0.5 mm voxels; both
EID and multi-bin PCD; beam-hardening off/on; RabbitCT for the realistic rabbit;
CL backprojector fixed **in the original class** (no new "Fixed" class);
RabbitCT published **v2 only**.
