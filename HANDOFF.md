# HANDOFF — SPIONvsXRay (updated 2026-07-09)

For the next agent. Read this, then `SPEC.md` (the plan, current), `README.md`
(state), `DEVLOG.md` (chronological). The user is **akmaier** (Andreas Maier, the
CONRAD author) — works on `main`, both `git@github.com:akmaier/SPIONvsXRay.git` and
`https://github.com/akmaier/CONRAD.git`.

## 1. What this is
Simulation study: **at what tumor iron concentration do PAA-coated magnetite SPIONs
become visible in a C-arm CT**, and does realistic melanoma uptake reach it?
Fully **CONRAD-native** via `pyconrad`. Reference: Heinen et al., *Ultrason.
Sonochem.* 130 (2026) 107876 (`1-s2.0-S1350417726001410-main.pdf` + supplement
`...-mmc1.pdf`, both in repo root). SPIE Medical Imaging paper (`paper/spie_manuscript.tex`).

## 2. Reconstruction — rebuilt CONRAD-native this session (DONE, committed + upstreamed)
2D **fan-beam** FBP is the detectability reconstructor. Geometry (`src/conrad_ct.py`
`fan_geometry`): CONRAD convention — `focal = SID = 750 mm`, virtual detector at the
isocenter, **`deltaT = 1 mm` (DETECTOR_DT_MM), decoupled from the recon voxel** (the
ramp kernels are only deltaS-exact at deltaT=1; the grid projector assumes 1 px = 1 mm).
No magnification/1.15 factors. Recon 512×512 @ 0.5 mm, short scan (Parker 180°+2γ),
500 views. Chain: forward → ParkerWeights → CosineFilter → SheppLoganKernel →
distance-weighted FanBeamBackprojector2D, **water precorrection always applied**.

Forward projection (`src/conrad_project.py` `project_base_materials`, `method=`):
`aa` (anti-aliased grid + `FanBeamProjector2D.projectRayDrivenCL`, DEFAULT),
`analytic` (CONRAD `FanBeamAnalyticProjector2D`), `disk_chords` (numpy oracle).

### CONRAD upstream contributions (akmaier/CONRAD master, all pushed)
- `731d345` reinstate fan 1/U² distance weighting (FanBeamBackprojector2D CPU+CL) — fixes cupping.
- `6817da2` SheppLoganKernel deltaS²→deltaS fix; `3ad5157` same for RamLakKernelLinux.
- `7152c33` **FanBeamAnalyticProjector2D** (+ example) — analytic per-material fan projector (ray tracer → MultiChannelGrid2D).
- `9018413` **WaterPrecorrectionTool** — single-material water precorrection (IndividualImageFilteringTool); matches SPION's numpy poly to 8e-16.
- All mirrored in `conrad_ext/` (built by `scripts/build_conrad_ext.sh`); SPION mirrors: `96db04e`, `d0955a8`, conrad_ext copies committed alongside.

### 3D cone-beam FDK — CLOSED, verified
`scripts/cone_resolution_study.py` (`250d002`): DC = **1.000 scale-invariant** across
voxel {4,2,1,0.5} mm × detector {2,1,0.5} mm, both kernels; only axial FDK cone
shading; fixed SheppLogan == RamLak. **Key fix:** `ConeBeamProjector.projectRayDrivenCL(Grid3D)`
returns all-zero on JOCL/Apple (host buffer never uploaded) — use
`fastProjectRayDrivenCL(OpenCLGrid3D, OpenCLGrid3D)` after
`gridCL.getDelegate().prepareForDeviceOperation()`. 0.25 mm detector cells FAIL with a
32-bit int-overflow ("Negative capacity") in CONRAD's OpenCL grid alloc — deferred.

## 3. Materials (corrected this session)
- Body + tumor host = **ICRU soft tissue** (`config.BODY_MATERIAL = "body"`, ρ=1.0),
  NOT water (fixed a drift).
- SPION material (`src/conrad_phantom.py` `register_spion_materials`): magnetite
  (Fe₃O₄) core + PAA coating, built in the **mass basis**. GOTCHA fixed: CONRAD
  `WeightedAtomicComposition.getCompositionTable()` returns per-element **mass** but
  `WAC.add()` expects **moles** — mixing them zeroed the iron. Build in mass basis.
- **Full particle simulated:** `c_NP = c_Fe/(0.724·(1−φ))`; PAA (C₃H₄O₂) registered as
  a CONRAD C/H/O material. Coating is low-Z → negligible for μ (~+3.6 % ΔHU at φ=0.15).
- **Coating fraction φ — decided value, NOT yet coded per-formulation.** The
  supplement TGA (Table A.1, `...-mmc1.pdf` in repo root) fixes φ per formulation:
  SPION I (12 nm) inorganic residual 87.7 % → magnetite 84.8 % → **φ ≈ 15 %**;
  SPION II (8 nm) 66.4 % → 64.2 % → **φ ≈ 36 %** (after the +3.4 % magnetite→hematite
  oxidation, article Eq A.1). **BUT `config.py` still hardcodes a single
  `PAA_MASS_FRAC = 0.15`** (exact for SPION I) and the phantom applies it to every
  insert — so per-formulation φ is a REMAINING task, folded into the Study A
  per-configuration sampling (§6/§7). README/SPEC document the per-formulation values
  as the model; the code has not been switched yet.
- Iron μ verified: mono ΔHU at c10 (0.543 mg Fe/ml) ≈ **+3.6 HU** (matches SPEC §5.2
  ~3 HU); `SPION_c0 == soft-tissue body` exact.

## 4. Detectability study (`src/run_detectability.py`)
Factors: **detector** (EID/PCD) × **bone** (absent/present) × **dose** (low 20 000 /
high 100 000 ph/px) × **distribution** (homogeneous/vessel) × **iron concentration**;
30 noise reps. Always-on (NOT factors): water precorrection, short scan, 90 kVp
spectrum. Spectrum optimization (kVp/filter) is the SEPARATE `run_spectral_sweep.py`.
Metrics: per-insert iron ΔHU (c0-corrected, local annulus), CNR, Rose-3/5 thresholds.
**Terminology: "detectability study" / "factor" — NEVER "factorial"** (the user objects).

## 5. Concentration / two-experiment model (IN PROGRESS)
Central question: at what tumor mg Fe/ml are SPIONs CT-visible, and does real melanoma
uptake reach it? Anchor concentrations to the article, not an arbitrary dose.
- **Study A — homogeneous (cellular uptake):** iron internalized by tumor cells →
  uniform. Sample the article's measured **cellular loading (Fig 5), all configs incl.
  fresh:** SPION I-113 **8.23 / 3.78** pg Fe/cell (0 h / 24 h), SPION II-115 **3.86 /
  1.52**, II-98 **3.60 / 1.37**, II-76 **2.51 / 0.85**. Convert pg Fe/cell → tumor mg
  Fe/ml via tumor cell density; each config uses its formulation's coating (I 15 %, II 36 %).
- **Study B — vascular / "fresh" (direct injection):** freshly-injected SPIONs still in
  the **blood carrier inside vessels** (150 µm vessels, ~10 % tumor volume, 10× local),
  at the injection concentration — heterogeneous, pre-uptake. Reference the Genç/Lyer
  flow-accumulation paper (suspension ~0.84 mg Fe/mL).

## 6. OPEN DECISIONS (blocking the concentration implementation — ASK THE USER)
1. **Study A cell density** for pg Fe/cell → mg Fe/ml (e.g. ~1–3×10⁸ cells/cm³; 8.23 pg
   × 3e8 → ~2.5 mg Fe/ml).
2. **Study B injection concentration** — blood-carrier level (article suspension 1–10 mg
   Fe/ml, or a specific injected dose).

## 7. REMAINING TASKS
1. Get the two open decisions, then implement in `config.py`/`run_detectability.py`:
   Study A multi-config sampling (8 configs, per-formulation coating) and Study B
   blood-vessel injection model (blood carrier in vessels, 10 %/10×, injection dilution).
2. Regenerate paper figures (`make_paper_figures.py`) + dashboard assets from the
   corrected `results/detectability.json`; refresh README headline numbers.
3. Verify the full-particle ΔHU per formulation (SPION II's ~36 % coating adds more
   low-Z mass — confirm ΔHU stays iron-dominated).
4. (Deferred/optional) 0.25 mm detector int-overflow: upstream `long`/size_t buffer fix.

## 8. Key files
`src/conrad_ct.py` (geometry+fbp), `src/conrad_project.py` (forward + EID/PCD +
measure_inserts), `src/conrad_phantom.py` (phantom + materials), `src/config.py`
(factors, FE_FRACTION=0.724, PAA_MASS_FRAC, dose levels), `src/run_detectability.py`,
`src/run_spectral_sweep.py` (separate), `conrad_ext/` (shadowed CONRAD, built by
`scripts/build_conrad_ext.sh`), `scripts/cone_resolution_study.py` (3D). Docs:
`SPEC.md`, `README.md`, `docs/index.html`.

## 9. Working style the user expects (IMPORTANT — this is why the session went long)
- **CONRAD API for ALL recon steps** — no hand-rolled ramp/geometry/projection.
- **Show images AND numbers**, not diffs or prose. The user is a 20-yr CT expert who
  reads results himself and distrusts self-assessed "looks good."
- **Step-by-step, verify together before code actions.** Look at results, then act.
- **Commit + push regularly** (both repos). Pushing to `main` is the established
  workflow (the user directs it), though subagents trigger a generic security warning.
- **Use subagents to keep moving** (don't block on background experiments); the user
  reviews the results. Spawn read-only investigations and verified implementations.
- Be honest about failures (report OOM/overflow/regressions, don't hide zeros).

## 10. Gotchas
- `deltaT = 1 mm` is REQUIRED (ramp kernels deltaS-exact only at 1; grid projector 1 px=1 mm).
- CONRAD material WAC: `getCompositionTable()` = MASS, `WAC.add()` = MOLES — build in mass basis.
- 3D: `fastProjectRayDrivenCL` + `prepareForDeviceOperation`, not `projectRayDrivenCL`.
- Iron audit: the study reports **iron** mg Fe/ml (contrast driver); particle mg SPION/ml
  is derived (`×1/(0.724(1−φ))`). The article's 1–10 mg/ml is **iron** (ICP-MAES).
- `tnan-2025-ra-0137-File003.docx` in root is the WRONG paper (Genç/Lyer gold-SPION) —
  relevant only as Study B flow context; the Heinen supplement is `...-mmc1.pdf`.
