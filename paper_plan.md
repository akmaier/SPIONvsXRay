# Paper Plan — SPIE Medical Imaging

Target conference paper based on the SPION CT-visibility simulation study
(`SPEC.md`). This document holds the venue/logistics, the paper structure, and a
**writing-style analysis of Andreas Maier** to keep the manuscript in his voice.

---

## 1. Target venue

- **Conference:** SPIE Medical Imaging (annual, San Diego, ~February).
  <https://spie.org/conferences-and-exhibitions/medical-imaging>
- **Track:** *Physics of Medical Imaging* (confirmed) — the natural home for CT
  simulation, spectral / photon-counting detectors, and detectability studies.
- **Format:** SPIE Proceedings, single-column, letter paper, no page numbers;
  built with the official `spie.cls` v3.4. Template downloaded and verified to
  compile under `paper/template/` (`spie_template.tex` = the official style guide;
  `spie.cls`, `spiebib.bst` = class + bibliography style).
- **Length:** SPIE proceedings are typically **6–10 pages**; aim for ~8.

**Deadlines — CONFIRM on the SPIE MI 2027 call-for-papers** (dates move yearly;
the pattern is): abstract submission ~**August 2026**, accepted-author
notification ~October 2026, **manuscript due ~January 2027**, conference
**February 2027**. → Set the internal schedule off the confirmed manuscript date.

## 2. Working titles (pick/refine)

1. *Are magnetic nanoparticles visible in CT? A simulation study of SPION
   detectability with energy-integrating and photon-counting detectors*
2. *Simulated CT visibility of superparamagnetic iron-oxide nanoparticles: dose,
   detector, and beam-hardening effects*
3. *How much iron does it take? Detectability limits of SPIONs in C-arm CT*

Recommended: **#1** — poses the question (Maier-style hook) and names the method.

## 3. Authors & affiliations (confirmed)

- **Andreas Maier** (corresponding) — Pattern Recognition Lab, FAU
  Erlangen-Nürnberg — simulation & methodology.
- **Stefan Lyer** — SEON, Section of Experimental Oncology and Nanomedicine,
  Dept. of Otorhinolaryngology, Universitätsklinikum Erlangen / FAU.
- **Lukas Heinen** — SEON, Dept. of Otorhinolaryngology, Universitätsklinikum
  Erlangen / FAU — SPION synthesis & characterization (reference article).
- **Rainer Tietze** — SEON, Dept. of Otorhinolaryngology, Universitätsklinikum
  Erlangen / FAU — nanoparticle characterization.

(Exact affiliation strings/order to be finalized on the manuscript; SEON authors
provide the material and dose grounding from the reference article.)

## 4. Core contribution (one-sentence thesis)

*Using a physically grounded CONRAD simulation of a rabbit-scale phantom, we
quantify whether biologically realistic SPION iron loads are detectable in C-arm
CT, and show how detector type (energy-integrating vs. multi-bin photon-counting)
and beam-hardening correction shift that detection limit.*

## 5. Paper structure (SPIE proceedings, ~8 pp.)

### 5.1 Abstract (~200 words)
Motivation (SPIONs for MDT/hyperthermia/cell tracking) → question (CT
visibility?) → method (CONRAD polychromatic sim, delivered-mass dose model,
factorial: dose × detector × beam-hardening, ΔHU + CNR) → headline result
(detection threshold in mg Fe/ml; EID vs PCD) → takeaway.

### 5.2 Introduction (~1–1.5 pp.)
- Why SPIONs matter clinically; why knowing their CT appearance matters
  (dose planning, incidental visibility, image-guided delivery).
- The problem: at biological loading the iron mass is small (~≤0.5 mg Fe/ml —
  ~10× below iodine enhancement); is it above the noise?
- Gap: little quantitative simulation of SPION CT detectability across modern
  detectors.
- **Contributions** (bulleted): (i) physically grounded delivered-mass dose
  model tied to a real SPION formulation; (ii) factorial detectability study
  EID vs multi-bin PCD with/without beam-hardening correction; (iii) reported
  detection thresholds (ΔHU, CNR/Rose).

### 5.3 Materials and Methods (~2.5–3 pp.)
- **2.1 Nanoparticle & dose model** — magnetite core (72.4% Fe); delivered-mass
  model anchored at 6 mg SPIONs for the 10 mg/ml formulation over an 8 cm³ tumor;
  `c_Fe = 0.0543·c_form`. Reproduce the conversion table from `SPEC.md` §5.2.
- **2.2 Digital phantom** — rabbit-scale ICRU soft-tissue cylinder + cortical-bone
  insert (beam-hardening source) + iron-loaded 8 cm³ tumor; 20 cm FOV.
- **2.3 X-ray simulation** — CONRAD standard polychromatic spectrum; standard
  C-arm cone-beam geometry, 500 projections; **70 000 photons/pixel** Poisson
  noise; two detector models — energy-integrating (EID) and energy-resolved
  multi-bin photon-counting (PCD).
- **2.4 Reconstruction** — cone-beam FDK; water beam-hardening correction on/off.
- **2.5 Evaluation** — factorial design (7 doses × 2 detectors × 2 BH × 10 noise
  reps ≈ 308 volumes); metrics ΔHU and CNR; detection threshold via
  Rose (CNR ≥ 3–5).

### 5.4 Results (~2 pp.)
- Example reconstructions (dose × detector grid).
- ΔHU and CNR vs. iron dose curves; EID vs PCD; BH on/off.
- Detection-threshold table (lowest visible mg Fe/ml per condition).
- PCD energy-bin / optimal-weighting benefit (if pursued).

### 5.5 Discussion (~1 pp.)
- Interpretation: is the realistic 6 mg dose visible? Which factor helps most?
- Comparison to iodine CT; spectral advantage of PCD for low-Z-contrast iron.
- **Limitations (state candidly):** simulation-only, idealized detectors,
  **round geometric (cylindrical) phantom rather than real rabbit anatomy** — no
  standard digital rabbit phantom exists (ROBY = rat), so beam hardening/scatter
  are idealized; homogeneous tumor uptake (vs the vessel-model Study B); single
  base geometry; no scatter/motion.

### 5.6 Conclusion (~0.3 pp.)
Answer the title question; state the threshold; note the detector implication;
one line of future work (experimental validation, spectral optimization).

### 5.7 Acknowledgments / References
CONRAD attribution; funding; `spiebib.bst` numeric style.

## 6. Figures & tables (plan for ~5 figs, 2 tables)

- **F1** Study-overview schematic: dose model → phantom → 500-proj scan → FDK →
  ROI metrics.
- **F2** Digital phantom (axial slice + labeled inserts; optional 3D).
- **F3** Reconstruction gallery: rows = dose, cols = {EID, PCD, BH on/off}.
- **F4** ΔHU vs. iron concentration (curves per detector/BH).
- **F5** CNR vs. dose with Rose threshold line; detection limits marked.
- **F6** (optional) PCD per-bin contrast / optimal energy weighting.
- **T1** Simulation parameters (spectrum, geometry, dose, voxel/detector).
- **T2** Detection thresholds (mg Fe/ml) per detector × BH.

## 7. Result-dependent claims (fill after simulation)

Leave placeholders; do **not** assert numbers until produced. Expected shape:
SPIONs at realistic loading are borderline in EID, improved by multi-bin PCD;
report the exact thresholds. Keep the paper honest if the result is "not
detectable at realistic dose" — that is itself a publishable, useful finding.

## 8. Template usage notes

- Edit a copy of `paper/template/spie_template.tex` → `paper/spie_manuscript.tex`
  (don't overwrite the pristine template). Keep `spie.cls`, `spiebib.bst` beside it.
- `\documentclass[]{spie}` (US letter) or `[a4paper]`; single spacing;
  `\title`, `\author[a]{}`, `\affil[a]{}`, `\authorinfo{}`, `abstract`,
  `\keywords{}`, then sections.
- Build: `latexmk -pdf spie_manuscript.tex` (verified working locally).

---

# 9. Writing-Style Analysis — Andreas Maier

Derived from *A Gentle Introduction to Deep Learning in Medical Image Processing*
(Maier et al., Z. Med. Phys. 29 (2019) 86–101; arXiv:1810.05401). The goal is to
write this SPIE paper in the same recognizable voice. Each trait below is paired
with a concrete rule for our manuscript.

### 9.1 Observed traits (with evidence)

1. **Reader-centric, didactic tone.** Explicitly lowers the barrier to entry:
   *"Readers of this article do not have to be closely acquainted with deep
   learning at its terminology."* Teaches rather than lectures.
2. **Intuition before formalism.** Concepts are motivated visually/biologically
   first, math second: the neuron is *"inspired by biological neurons"* and
   *"computes a weighted sum of its inputs"* before any equation appears.
3. **Analogies to ground abstraction** — apples vs. pears for feature extraction;
   a single neuron *"can already be interpreted as a classifier."*
4. **Inclusive, guiding "we".** *"We first discuss general reasons…"*,
   *"we still need to discuss how its parameters are actually determined."* The
   reader is walked through, step by step.
5. **Strong signposting / metacognition.** Transitions state where the argument
   is: *"Having gained basic insights into neural networks and their basic
   topology, we still need to discuss…"* Each move is announced.
6. **Motivation-first (why before how).** A limitation is named, its consequence
   explained, then the solution introduced, then grounded in theory
   (e.g., single neurons can't solve XOR → add layers → universal approximation).
7. **Varied sentence rhythm.** Short declaratives land key points
   (*"A single neuron itself can already be interpreted as a classifier"*);
   longer sentences build multi-clause explanations. Readable, not choppy.
8. **Notation introduced verbally.** Symbols are defined in words around the
   equation: *"the classifier has to predict the correct class y, which is
   typically estimated by a function ŷ = f̂(x)."*
9. **Figures woven into the narrative**, introduced before shown, cited inline:
   *"cf. Fig. 1"*, *"as summarized graphically in Fig. 3."*
10. **Honest, balanced framing.** Enthusiasm plus candor about limits:
    *"we also clearly indicate potential weaknesses of the current technology."*
    Hedges appropriately (*"we believe it is worthwhile…"*).
11. **Big-picture, forward-looking close.** Ends by situating the work in a
    fast-moving field and pointing ahead, rather than merely summarizing.

### 9.2 Distilled "Maier voice" rules for this paper

- **Open with the question, not the machinery.** Lead the abstract and intro with
  *why SPION CT visibility matters* and pose it as a crisp question
  ("Are they visible?") — mirror the title's hook.
- **Use "we" throughout**; guide the reader move by move.
- **Explain before you formalize.** Precede the `c_Fe = 0.0543·c_form` equation
  with one intuitive sentence ("we spread a fixed particle mass through the tumor
  and ask how much iron that leaves per millilitre").
- **Signpost every section transition** with a metacognitive sentence ("Having
  fixed the dose model, we now describe the phantom…").
- **Lead each paragraph with a topic sentence** stating its single point.
- **Vary sentence length deliberately** — short sentence for each key result,
  longer for mechanism.
- **Introduce every figure in prose before it appears**; cite inline as
  "cf. Fig. X"; make each figure carry one message.
- **Be candid about limitations** in the discussion — simulation-only, idealized
  detectors, homogeneous uptake. Maier's style treats honesty as credibility.
- **Motivation-first paragraph architecture:** limitation → consequence →
  our approach → evidence.
- **Close with the big picture**: what the detection limit means for
  image-guided nanoparticle therapy, and the path to experimental validation.
- **Accessible but rigorous:** define terms (SPION, EID, PCD, CNR, Rose
  criterion) on first use; keep equations minimal and verbally framed.

### 9.3 Anti-patterns to avoid (contrary to the style)

- Front-loading dense notation before motivation.
- Passive, impersonal "it was found that…" (prefer "we found").
- Figures dumped without narrative setup.
- Overclaiming — especially asserting SPIONs are "clearly visible" if the data
  say otherwise.

---

## 9b. Reproducibility & audit dashboard (supports the paper)

The project ships a public **GitHub Pages audit dashboard** (`docs/`,
<https://akmaier.github.io/SPIONvsXRay/>) that shows every intermediate result —
spectrum, material attenuation, phantom, projections, reconstructions,
detectability curves — plus **annotated code snippets for auditing**. For the
manuscript this gives us: (i) a citable reproducibility/data-availability
statement, (ii) ready-made figure sources, and (iii) reviewer-facing
transparency. The project is **MIT-licensed** (`LICENSE`), so code and figures
can be reused/cited. Add a "Code and data availability" line to the paper linking
the dashboard and repository.

## 10. Next actions

1. Confirm SPIE MI 2027 dates + track; confirm co-authors.
2. Run the simulation study (`SPEC.md` M0–M7) → produce F3–F6 + tables.
3. Draft `paper/spie_manuscript.tex` from the template, applying §9.2.
4. Internal review against the style rules before submission.
