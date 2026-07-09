# Development Log — SPIONvsXRay

Reverse-chronological log of progress. Newest entries on top.

---

## 2026-07-09 — RabbitCT .rctd format cracked; real 3D geometry extracted

- Reverse-engineered the undocumented RabbitCT `.rctd` binary (`src/rabbitct.py`):
  24-byte header + per-projection [96-byte column-major 3×4 double matrix + 1248×960
  float32 image]. Decomposed the 496 matrices into the real Siemens C-arm short scan:
  **SID 745 mm, SDD ~1196 mm (0.308 mm pitch), 496 views over 198° (0.40°/step)**,
  rotation about z, offset detector. Loaded `reference_256.vol` (256³ rabbit head/neck
  anatomy — renders correctly) and inserted a SPION tumor (3D phantom PoC).
- Added the real rabbit geometry + anatomy figure to the paper (future-work bridge to
  the 3D study). Documented the format + next steps in HANDOFF §4.
- **Next: full 3D** — CONRAD cone-beam forward projection (ProjectionTable trajectory
  from the matrices) + FDK of the SPION-in-rabbit phantom, same detectability analysis.

## 2026-07-09 — CONRAD GPU spectral detectors (Layer 2) upstreamed

- Added on-device polychromatic detectors so the whole spectral projection stays
  on the GPU, pushed to **akmaier/CONRAD master (`fdf5365`)**:
  - `physics/detector/SpectralDetector.cl` — fused `energyIntegratingDetector` /
    `photonCountingDetector` kernels (self-contained Philox/Poisson).
  - `OpenCLSpectralDetector` (shared base: samples the spectrum, builds
    per-material linear-attenuation tables identical to
    `PolychromaticAbsorptionModel`, uploads once) + two subclasses
    `OpenCLEnergyIntegratingDetector` and `OpenCLPhotonCountingDetector` (kept
    separate — different detector concepts, per AM).
  - Noise is applied **per energy** on the photon counts (before energy weighting
    / binning) — physically correct: EID variance = ΣE²·n, which Poisson(ΣE·n)
    cannot reproduce (AM's correctness requirement).
- **Verified on the M1 GPU** vs a reference from CONRAD's own spectrum + material
  attenuation: EID matches to 1.4e-7, PCD bins exactly; noisy EID variance
  = 1.25e8 ≈ ΣE²·n (1.26e8), not ~mean (1.9e6).
- Regenerated the full re-release jar (`publish/conrad/`, sha `6d41971c3f5d`) with
  all four contributions (exp/log, backprojector, RNG, spectral detectors);
  `scripts/rebuild_conrad_jar.sh` updated. Hosting still pending (see RELEASE.md).

## 2026-07-09 — CONRAD GPU noise generators (Layer 1) upstreamed

Goal (AM): give CONRAD OpenCL-native noise generators, photon counter, and
detector models so the whole spectral projection can run on-device. Staged
Layer 1 first (noise ops), verified, pushed to **akmaier/CONRAD master (`631fbed`)**.

- **No GPU RNG existed in CONRAD** — noise was CPU-only (`StatisticsUtil`).
- Added a **Philox4x32-10 counter-based RNG** (Random123) in
  `PointwiseOperators.cl`: stateless, reproducible via a seed, no per-thread
  state buffers. On it: `poisson` kernel (element = mean λ; Knuth for λ<30,
  Hörmann **PTRS** transformed rejection for λ≥30 — exact across the range) and
  `standardNormal` (Box–Muller).
- CPU reference impls in `NumericGridOperator` (seeded `Random` + Lanczos
  lgamma), GPU overrides in `OpenCLGridOperators`, static wrappers in
  `NumericPointwiseOperators` — so `poisson(grid,seed)`/`standardNormal(grid,seed)`
  dispatch to CPU or GPU by grid type like the other pointwise ops.
- **Verified on the M1 GPU** (test jar, `conrad_ext` off): Poisson mean==var==λ
  from 0.5→5000 (Knuth + PTRS), integer output, per-element λ, reproducible;
  **chi-square goodness-of-fit vs the analytic PMF passes** (λ=5 and λ=50);
  standard normal mean 0, std 1, ~0 skew/kurtosis.
- **Next — Layer 2:** GPU polychromatic absorption model + photon-counting /
  noisy detector classes (grid-resident), then one combined jar re-release.

## 2026-07-08 — CONRAD OpenCL fixes upstreamed to akmaier/CONRAD + jar re-release prep

- Studied CONRAD's GPU stack to answer "can we run it all on GPU with OpenCL
  grids?" Findings: `OpenCLGrid2D/3D` + `OpenCLGridOperators` (GPU `multiplyBy`,
  `addBy`, `sum`, `pow`, …) + GPU-resident `fastProjectRayDrivenCL` /
  `fastBackprojectPixelDrivenCL` all exist, and `SimulateXRayDetector` (render
  material path-lengths once, apply any detector) is exactly our architecture —
  so an all-GPU, OpenCL-grid-resident pipeline is supported. **But** two GPU ops
  needed for the polychromatic conversion are broken: `OpenCLGridOperators.exp()`
  called a non-existent kernel `"expontial"` (no natural-exp kernel at all) and
  `log()` called a missing `"logarithm"` — both threw `CLInvalidKernelNameException`.
  (No GPU Poisson anywhere; noise stays a cheap CPU step.)
- Per AM: **fix upstream in CONRAD, not as a `conrad_ext` shadow.** Added
  `exponential` + natural-`logarithm` kernels to `PointwiseOperators.cl` and fixed
  the `exp()` call site; **ported the `FanBeamBackprojector2D` pixel-spacing fix**
  (previously only in `conrad_ext`) into the repo too. Pushed to
  **akmaier/CONRAD `master` (`71bb20c`)**.
- Built a patched `conrad_1.1.0.jar` by recompiling the two changed classes
  against the released jar and splicing the new `.class` + `.cl` in with `jar uf`
  (`scripts/rebuild_conrad_jar.sh`). **Verified** against the patched jar with the
  `conrad_ext` shadow disabled + OpenCL on: GPU `exp(log(x))=x` to 1e-7; 0.5 mm
  recon with monotonic iron ΔHU (c20 → +7.14, identical to the `conrad_ext` result).
- Re-release prep: swapped the patched jar into the canonical FAU zip
  (`publish/conrad/CONRAD_1.1.0.zip`; jar sha `c4aa2689bfef`, was `38fe24e116da`).
  **Pending (needs FAU access):** host it at
  `https://www5.cs.fau.de/fileadmin/user_upload/CONRAD_1.1.0.zip`. See
  `publish/conrad/RELEASE.md`. Once the re-hosted jar ships via pyconrad,
  `conrad_ext` can retire (keep it until then — our env still uses the old jar).

## 2026-07-08 — M6 PCD bug fixed + full factorial re-run (GPU 0.5 mm)

- **Fixed the PCD detector bug** (was CNR ≈ 0/negative, unusable). Two root
  causes in `src/run_factorial.py`:
  1. `_pcd_weights` weighted bins by **iron contrast only**, so the low-energy
     bin (10–37.5 keV, highest contrast) got the largest weight — but a 90 kVp
     spectrum has few photons there and ~11 cm of tissue removes almost all of
     them, so the code up-weighted the noisiest bin.
  2. Bins were combined as `Σ w_b·(−log C_b)` — a **per-bin log** that diverges
     when the starved low bin's Poisson counts hit 0 (−log(0→ε) ≈ 22), which FBP
     smears into image-wide streaks (per-realization ΔHU swung ±10 HU).
  Fix: combine **counts first, then a single log** —
  `p = −log(Σ w_b C_b / Σ w_b C_air,b)` — which is dominated by the populated
  high bin and never blows up, with matched-filter weights `w_b ∝ S_b/V_b`
  (`S_b = Σ_bin N_t·c`, `V_b = Σ_bin N_t`; detected-count-weighted mean iron
  contrast). Verified against a **detector-level Monte Carlo**: the estimator
  reaches 1.29× EID, essentially the 1.34× ideal-observer ceiling from
  `spectral.py`.
- **Re-ran the full factorial** on the GPU 0.5 mm pipeline, 30 noise
  realizations/cell (was CPU/1 mm before), persisted to `results/factorial/`
  (`factorial.csv` + `factorial.json`; added a `save_results` writer). Results:
  - Detection threshold ≈ **0.54 mg Fe/ml** (Rose 3 and 5) — the realistic 6 mg
    dose — **identical for all four** EID/PCD × BH off/on cells.
  - Realistic dose (c_Fe 0.54): EID ΔHU 3.4 / CNR 5.4; PCD ΔHU 4.7 / CNR 5.1.
  - 2× dose (c_Fe 1.09): EID CNR 8.8; PCD CNR 11.5 = **1.31× EID** (PCD's ideal
    spectral gain shows only once the low bin isn't photon-starved).
  - Beam-hardening correction: negligible effect on iron detectability.
- Bumped `EVAL.noise_realizations` 10 → 30 for stable CNR. Refreshed README key
  findings + HANDOFF headline (were still carrying old CPU/1 mm numbers).
- **Still TODO:** filter/kVp sweep in the factorial (accumulators already accept
  `kvp`/`filters`); RabbitCT geometry extraction; 3D rabbit case.

## 2026-07-08 — M3 COMPLETE ✅ (real CONRAD spectrum + refined optimization)

- `src/spectrum.py`: wraps `PolychromaticXRaySpectrum`. The no-arg ctor = CONRAD
  **standard spectrum = 90 kVp** (not 120 as placeholder), 10–150 keV @0.5,
  mean 55.4 keV, W anode, peak flux at the 59.5 keV characteristic line.
  Parameterized `conrad_spectrum(kvp)` + `apply_filters` for the sweep.
- Rewired `spectral.py` to use the **real spectra**. Refined results:
  60 kVp = 1.34× / 80 = 1.09× / 120 = 0.85× the 90 kVp ideal CNR; hardening
  filters still hurt (Sn 0.58×). On the softer 90 kVp, EID already gets 72% of
  ideal, so optimal PCD weighting = **1.35× EID** (3-bin @ **37.5/50 keV**,
  97% of ideal; 2-bin @ 47.5). Updated config (kVp=90, thresholds) + SPEC §5.8.
- Dashboard spectrum/pcd_bins figures now show the real CONRAD spectrum.

**Next (M4):** 500-proj cone-beam forward projection (polychromatic, EID +
multi-bin PCD) using the phantom component volumes + real spectrum.

---

## 2026-07-08 — CL fan backprojector FIXED in-place + 0.5 mm GPU recon

- Subagent root-caused + fixed the OpenCL fan backprojector (kernel never got
  pixel spacing + `//FIXME` sign/guard hacks). Verified: matches CPU ~1e-5,
  spacing knob exact (0.5 mm marker ratio 1.998), ~850× faster.
- Per user ("we have versioning, no new class"): folded the fix into the
  **original** `edu.stanford.rsl.tutorial.fan.FanBeamBackprojector2D` + its
  `FanBeamBackProjectorPixel.cl` (added `setSpacing`, spacing kernel args, and
  the corrected 2DCL detector geometry; 1DCL untouched). Source lives in
  `conrad_ext/` (git-tracked, contributable); `scripts/build_conrad_ext.sh`
  compiles it against the 1.1.0 jar; `conrad_backend` puts `conrad_ext/out` on
  the classpath (CONRAD_DEV_DIRS) so it **shadows** the jar class. Verified the
  patched class loads from `conrad_ext/out/` and has `setSpacing`.
- `conrad_ct.fbp` now uses the GPU backprojector with `setSpacing(voxel_mm)`;
  `config.RECON_VOXEL_MM = 0.5`. **Full native pipeline at 0.5 mm on GPU** →
  correct monotonic ΔHU (c20 → +7.14, matching the CPU/1 mm result). This solved
  BOTH the CL-backprojector fix and the 0.5 mm voxel request; the separate
  "config-aware reconstruction" route is no longer needed.
- Dashboard evidence: clbp_fixed_disk.png, clbp_spacing_marker.png.

**Geometry plan (user):** rabbit case = full 3D cone-beam at the real RabbitCT
C-arm geometry (496 proj matrices); disc phantom study = 2D fan matching the same
C-arm (SID/SDD/detector/angular range). Feasible; next.

---

## 2026-07-08 — RabbitCT PUBLISHED on Zenodo (DOI 10.5281/zenodo.21267885, CC BY 4.0)

Uploaded the RabbitCT -v2 dataset from lme31 (python3 uploader; lme31 lacks
curl/jq) to a Zenodo **draft** (id 21267885, PRODUCTION, not published):
rabbitct_512-v2.rctd (2.5 GB) + rabbitct_1024-v2.rctd (4.5 GB) + reference_256.vol
+ dev kit + README. Metadata: 4 authors (Rohkohl/Keck/Hofmann/Hornegger),
CC BY 4.0, isDocumentedBy doi:10.1118/1.3180956, post-mortem ethics note.
Decision: publish v2 only (originals not added). **PUBLISHED** by the user:
DOI **10.5281/zenodo.21267885** (concept 10.5281/zenodo.21267884),
record https://zenodo.org/records/21267885. Documented in README.md. Token stays
in git-ignored publish/rabbitct/zenodo_secrets.env.

---

## 2026-07-08 — RabbitCT data (real rabbit CT) located + fetched from LME

- Sven's tip: rabbit CT is the **RabbitCT benchmark** (Erlangen group), on lme31 at
  `/disks/.../Forschung/Software/RabbitCT`. Reached it via two-hop SSH
  (maier@cluster.i5.cs.fau.de → maier@lme31.cs.fau.de; key-based, works).
- Contents: `download/reference_256.vol` (67 MB = reconstructed **rabbit volume**,
  256³ @ 1 mm — the real anatomy we can use as the realistic phantom);
  `rabbitct_512/1024.rctd` (2.7/4.8 GB benchmark projections — not fetched);
  `rabbitct_develop.zip` (format spec + geometry). CONRAD has RabbitCT support
  under the **"LolaBunny"** codename (`LolaBunnyBackprojector`).
- Fetched `reference_256.vol` + dev kit into `data/rabbitct/` (git-ignored).
- RabbitCT is a real **C-arm geometry** (496 projections, 0.5 mm-voxel recon) →
  ideal for the "set the geometry correctly" Conrad config.
- Two subagents dispatched: fix the CL fan backprojector (`// TODO: Spacing`); and
  config-aware reconstruction (rate-limited on first try, to re-run).

---

## 2026-07-08 — Investigation: CL backprojector + voxel spacing (user questions)

**CL backprojector — root cause found.** `FanBeamBackprojector2D.backprojectPixelDrivenCL`
never passes the image pixel spacing to its kernel — literally `// TODO: Spacing :)`
in the source. So its pixel↔world mapping is inconsistent with the CPU path:
looks fine on a uniform disk (differs by a scale) but **corrupts sub-HU
measurements** (per-insert ΔHU came out 196 HU, non-monotonic). Reverted: forward
projection uses CL (validated), **backprojection stays CPU** (correct). GPU thus
does NOT accelerate the reconstruction step (the factorial bottleneck) as-is.

**Voxel spacing — where it lives.** CONRAD's *global Configuration* DOES hold
voxel spacing (`Configuration.java` default `setVoxelSpacingX(1.0)`; `Trajectory`
has get/setVoxelSpacingX/Y/Z + getReconDimension*), and the MAIN reconstruction
pipeline reads it. BUT the **tutorial** `FanBeamBackprojector2D` we use ignores
the config and hardcodes 1 mm via `Translation(-imgSize/2,-imgSize/2,-1)`. So the
observed 1.0 mm is the tutorial hardcode (coincidentally == the config default).
To get 0.5 mm properly: set the global Configuration + use a config-aware
reconstructor, or a self-contained geometry rescale.

**Factorial (24 cells, CPU backproj, 17 min):** EID CNR rises to ~3.9 at the top
dose (1.09 mg Fe/ml) — just below Rose 5 → SPIONs **borderline-undetectable**
(the expected headline). PCD path is currently BROKEN (CNR≈0/neg) — a bug in my
per-bin combination to fix.

---

## 2026-07-08 — CONRAD-native pipeline COMPLETE end-to-end (src/conrad_project.py)

Full native chain working: CONRAD AnalyticPhantom → PriorityRayTracer fan-beam
**base-material sinograms** (9 materials, 500×1024, ~17 s, 38 µs/ray) →
polychromatic combine (per-material CONRAD attenuation, real 90 kVp spectrum,
EID + multi-bin PCD, Poisson noise) → CONRAD fan-beam FBP → per-insert ROI.
- Base-material path lengths exact (water 160 mm, each insert 25 mm, bone 25 mm).
- Fixed a ROI-scale bug: FanBeamBackprojector2D writes a Grid2D at **1.0 mm/px**
  (default), not geo["spacing"]; calibrated via the bone marker.
- **Result (noise-free EID):** high-c ΔHU monotonic and matches the earlier
  numpy pipeline (c20 → +6.9 HU vs +7.6). BUT low-c inserts sit at a spurious
  ~−2.5 HU baseline from **bone-insert streak artifacts** — i.e. the sub-HU iron
  signal is below the reconstruction-artifact floor. Real finding + an M6
  measurement-refinement TODO (same-radius annular reference / more views /
  apodization / bone-free contrast sub-scan). Figure → docs/assets/recon_grid.png.

**Next (M6):** refine per-insert detectability (kill streak bias), then run the
full factorial (detectors × BH × noise × filters) + dashboard.

---

## 2026-07-08 — CONRAD-native phantom built (src/conrad_phantom.py)

- Real CONRAD `AnalyticPhantom` (PrioritizableScene) built at runtime via pyconrad:
  water body Cylinder + 7 SPION insert Cylinders on a circle + bone insert (9
  PhysicalObjects). Inserts added after the body → PriorityRayTracer overrides.
- Registered magnetite SPION `Mixture` materials (H2O + Fe3O4 per c_Fe) via the
  CreateCustomMaterial pattern (WeightedAtomicComposition + MaterialsDB.put).
- **Cross-validated:** CONRAD material μ@60keV rises 0.20584 (c=0=water) → 0.20723
  (c=20), matching the independent materials.py oxide model to 5 decimals.
- Next: src/conrad_project.py (PriorityRayTracer fan-beam base-material sinograms).

---

## 2026-07-08 — PATH A IMPLEMENTED ✅ (CONRAD OpenCL on the M1 GPU)

Wired Path A (user go-ahead). Now working:
- `scripts/install_opencl.sh`: fetches jogamp 2.6.0 (universal arm64 natives)
  into the pyconrad bundle dir (globbed onto the classpath before conrad_1.1.0.jar
  → shadows the stale 2.3.2 JOCL) and builds the `.jogamp/libOpenCL.dylib`
  reexport shim for Apple's OpenCL.framework.
- `conrad_backend.setup()` sets DYLD_LIBRARY_PATH to `.jogamp` before JVM start;
  `conrad_backend.opencl_available()` probes `OpenCLUtil.getStaticContext()`.
- **Verified:** OpenCL init → Apple M1 Max, 32 CUs, GPU. CONRAD CL fan **forward
  projector** matches CPU to **0.03%** and is **~4000× faster** (88.4 s → 0.02 s
  for a 512² × 500-view projection). `conrad_ct.project(...)` uses it with CPU
  fallback.
- **Caveat:** `backprojectPixelDrivenCL` reconstructs incorrectly here (output
  inverted vs CPU — a convention/setup mismatch); FBP stays on the validated CPU
  backprojector for now (TODO(opencl) to fix). Forward projection was the
  bottleneck, so the win is captured.

---

## 2026-07-08 — OpenCL PATH A PROVEN VIABLE (subagent investigation)

Correction to the note below: the true root cause is a **pure-Java gluegen guard**
(`RuntimeException: Please port CPU detection to your platform (mac os x/aarch64)`)
that throws in a static initializer before any native loads — the 2015 jogamp has
no aarch64 CPU-detection case. (Wrong-arch natives were only a secondary issue.)

**Path A proven working on this M1 Max** (subagent tested end-to-end):
- jogamp **2.6.0** ships universal (x86_64+arm64) gluegen/JOCL natives.
- Apple's OpenCL.framework isn't on JOCL's macOS search list → add a 1-line
  reexport shim `libOpenCL.dylib` (`clang -dynamiclib -Wl,-reexport_framework,
  OpenCL -framework OpenCL`) discoverable via `DYLD_LIBRARY_PATH`.
- With that: CONRAD's own `OpenCLUtil.getStaticContext()` runs on the 32-CU GPU
  and a CONRAD OpenCL kernel executes correctly. **No CONRAD recompile** (JOCL
  high-level API binary-compatible 2.3.2→2.6.0). Non-invasive: prepend the 2.6.0
  jars + set DYLD_LIBRARY_PATH in conrad_backend.setup().
- GPU projectors are drop-in: `projectRayDrivenCL` / `backprojectPixelDrivenCL`.

**Cost:** Path A ≈ ½ day, low risk, ~10–30× (2D) / ~50–100× (3D) payoff. Path B
(x86_64 Rosetta) worse. Path C (CPU) fine for THIS study, blocks 3D.
**Plan:** CPU path stays the baseline for the SPION study; wire Path A as an
optional GPU accelerator (pending user go-ahead) since it pays off long-term.

---

## 2026-07-08 — OpenCL check: not available on this Mac (arch, not hardware)

User asked whether CONRAD's OpenCL runs here (it's OpenCL, not CUDA). Tested:
- Hardware present: Apple M1 Max, 32 GPU cores, Apple OpenCL 1.2 exists.
- But CONRAD OpenCL **fails to initialize**: JOCL Java classes load, native
  binding does not — `NoClassDefFoundError: Could not initialize class
  com.jogamp.opencl.llb.impl.CLAbstractImpl` / `ExceptionInInitializerError`.
  Cause: pyconrad bundles jogamp **2.3.2 (2015)** → gluegen/JOCL natives are
  x86_64-only, incompatible with the arm64 JVM. Not a hardware problem.
- Fix options (not pursued): arm64 JOCL 2.4+ natives (API-mismatch risk vs CONRAD
  1.1.0) or a full x86_64/Rosetta stack. Not worth it for a 2D study.
- **Decision:** keep the proven CPU path (PriorityRayTracer base-material
  sinograms + CPU fan FBP). Corrected the earlier "no CUDA" wording in conrad_ct.

---

## 2026-07-08 — CONRAD-native path proven (analytic phantom + material projector)

User steered to: use a CONRAD Phantom + CONRAD's material projectors; checked out
the CONRAD source. Findings & proofs:
- **Cloned CONRAD source** to `~/Documents/CONRAD` (github.com/akmaier/CONRAD) for
  reference + to host the new phantom. (No Maven/Gradle — Eclipse project; a single
  new class can be compiled against the pyconrad jar instead of a full rebuild.)
- **Material projector = `MaterialPathLengthDetector`** (records per-material path
  length into a `MultiChannelGrid2D` = base-material sinograms) +
  `XRayDetector.accumulatePathLenghtForEachMaterial(segments)`.
- **Analytic ray tracer = `PriorityRayTracer`** (NOT WatertightRayTracer, which is
  STL/mesh-only). `setScene(phantom)` + `castRay(StraightLine)` → material segments.
- **PROVEN** via pyconrad: built an AnalyticPhantom (water Cylinder + iron insert
  Cylinder), cast a ray → exact base-material path lengths (iron 25.00 mm = insert
  dia; water 135.00 mm = body − insert). CPU, headless, no Configuration.
- **`CirclesPhantom` is the template** for the new phantom (Cylinder body +
  PhysicalObject inserts on a circle + MaterialsDB materials + add()).

**New phantom location:** `CONRAD/src/edu/stanford/rsl/conrad/phantom/SPIONInsertPhantom.java`
(contributable), built identically at runtime via pyconrad for now.

**Build plan (next):** src/conrad_phantom.py — runtime AnalyticPhantom (ED inserts
+ registered magnetite SPION materials); src/conrad_project.py — PriorityRayTracer
fan-beam base-material sinograms; then polychromatic EID/PCD + CONRAD fan FBP.
Supersedes the numpy-geometry hybrid.

---

## 2026-07-08 — COURSE CORRECTION (user review): use CONRAD for recon; ED phantom

User caught four issues in the M4/M5 review — all valid:
1. **Not using CONRAD for projection/recon** — M4/M5 were custom numpy Radon/FBP
   (CONRAD only gave attenuation+spectrum). FIXED: proved CONRAD's **CPU**
   fan-beam works on this Mac — `FanBeamProjector2D.projectRayDriven` +
   `FanBeamBackprojector2D.backprojectPixelDriven` + numpy ramp/cosine weighting
   (`src/conrad_ct.py`). Round-trip on a disk reconstructs sharp & flat. No GPU.
2. **Noise not visible** — it IS simulated (~17 HU), but the montage window
   spanned the bone spike (4× body) so it was invisible; will window to soft tissue.
3. **Phantom design** — switching to an **ED-phantom style**: concentration inserts
   on a circle at equal radius (all imaged in one scan, cupping equal, directly
   comparable), each insert a 2.5 cm (8 cm³-equiv) disk. Replaces the single tumor.
4. **Persist sinograms** — will save sinogram arrays to results/ and show them.

Rework in progress: phantom.py (ED), simulate.py/reconstruct.py -> conrad_ct.

---

## 2026-07-08 — M2 COMPLETE ✅ (round geometric phantom)

- `src/phantom.py`: analytic geometry (soft-tissue cylinder + bone rod + tumor
  sphere) voxelized into 3 component volumes (soft frac, bone frac, iron mg/ml)
  so M4 can integrate each material's path length independently for polychromatic
  projection. Iron modeled as magnetite downstream.
- Tumor models: **homogeneous** and **vessel** (150 µm vessels @10%, mass
  conserved via per-voxel Binomial partial-volume sampling).
- Verified: homogeneous tumor mean = 1.086 mg/ml (= tumor_iron_conc(20), std≈0);
  vessel mean = 1.083 (mass conserved), std = 0.274 (~25% structural noise),
  max 2.31 (local enrichment). Dashboard `phantom_axial.png` refreshed.

**Next (M3):** extract the real CONRAD standard polychromatic spectrum; then M4
projector.

---

## 2026-07-08 — CORRECTION: ROBY is a RAT; no digital rabbit exists

- Verified via literature: RADAR realistic animal series (Keenan/Stabin/Segars,
  JNM 2010;51:471) = **MOBY (mouse) + ROBY (RAT)** only. **No rabbit** in the
  RADAR/Segars/XCAT family. My earlier "ROBY = rabbit" was wrong.
- MOBY/ROBY are licensed via **Duke CVIT** (paul.segars@duke.edu,
  https://cvit.duke.edu/resource/moby-roby-phantoms/), NURBS + voxel — NOT from
  doseinfo-radar.com (that only describes them; cert broken).
- A realistic **rabbit** needs a **real rabbit CT** (one-off in the literature,
  e.g. PubMed 25971772 / 34487228). SEON co-authors may have one.
- Corrected SPEC §5.1 and config: arm (B) background is now **pending user
  decision** (real rabbit CT [preferred] vs ROBY rat [rat-scale, licensed] vs
  drop). Batch runs on (A) round geometric meanwhile. Config
  `PHANTOM_BACKGROUNDS=['round']`, `REALISTIC` paths for rabbit-CT/ROBY.

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

- **ROBY format check (introspected the CONRAD jar):** ROBY = same Segars NURBS
  spline format as XCAT. CONRAD reads it natively via `XCatScene` (loads splines,
  spline→material LUT, tessellate) + `AnalyticPhantomProjector`. BUT the loader is
  XCAT-organ-specific → ROBY needs either (1) a ROBY spline→material LUT adapter
  or (2) the STL/mesh route via `AsciiSTLMeshPhantom`. Preferred: analytic NURBS.
  **Q for user:** do you have ROBY as `.nrb` NURBS splines (→ analytic) or only
  voxelized output (→ voxel adapter)? Determines which loader I build.

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
