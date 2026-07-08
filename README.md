# SPIONvsXRay

**Can superparamagnetic iron-oxide nanoparticles (SPIONs) be seen in CT?**

This project simulates the X-ray absorption and cone-beam CT reconstruction of
PAA-coated iron-oxide nanoparticle suspensions at varying iron concentrations,
to determine whether — and at what concentration — they become visible in a
standard C-arm CT. The simulation uses the [CONRAD](https://www5.cs.fau.de/conrad/)
framework via [`pyconrad`](https://pypi.org/project/pyconrad/).

The nanoparticles are those described in Heinen et al., *Ultrasonics
Sonochemistry* 130 (2026) 107876 (PAA-coated SPION clusters, iron-oxide cores
~8–12 nm, used at ~1–10 mg Fe/ml). *The article PDF is git-ignored (copyright).*

## Idea in one picture

```
SPION suspensions          multi-vial phantom        C-arm cone-beam scan
@ 0,1,2,5,10 mg Fe/ml  ──►  (water + inserts)   ──►  500 projections, std spectrum
                                                             │
                                          HU vs. concentration, CNR  ◄── FDK reco
```

Because X-ray attenuation of a dilute suspension is dominated by its **iron
content**, each concentration is modeled as *water + Fe* at the matching mass
fraction (CONRAD custom-material "mixture" recipe), imaged polychromatically and
reconstructed, then the reconstructed contrast (Hounsfield Units) and
contrast-to-noise ratio are measured per insert.

## Method summary

1. **Materials** — register a CONRAD material per iron concentration.
2. **Phantom** — a water cylinder with rod inserts, one concentration each.
3. **Spectrum** — CONRAD standard polychromatic X-ray spectrum.
4. **Acquisition** — 500 projections in the default C-arm cone-beam geometry.
5. **Reconstruction** — cone-beam FDK.
6. **Analysis** — HU-vs-concentration curve, CNR, detectability threshold.

See [`SPEC.md`](SPEC.md) for the full plan and milestones, and
[`DEVLOG.md`](DEVLOG.md) for progress.

## Status

🚧 Planning complete; environment setup next. This repo currently contains the
plan (`SPEC.md`), this README, and the development log (`DEVLOG.md`). Code lands
under `src/` per the milestones in the SPEC.

## Getting started (planned)

CONRAD requires **Java 8** (already present) and a pyconrad-compatible Python
(3.9/3.10 recommended — the system Python 3.14 is too new for JPype):

```bash
python3.10 -m venv .venv && source .venv/bin/activate
pip install pyconrad
python -c "import pyconrad; pyconrad.setup_pyconrad()"   # downloads CONRAD.jar, starts JVM
```

## Missing information (to finalize the experiment)

The plan uses reasonable defaults, but the following would change results and
should be confirmed:

1. **Iron-oxide chemistry & density** — magnetite Fe₃O₄ vs maghemite γ-Fe₂O₃
   (affects particle-mass ↔ Fe-mass conversion; X-ray model uses Fe directly).
2. **Concentrations & units** — exact list to sweep, and **mg Fe/ml** vs
   mg-particle/ml. In-vitro/clinical suspension range, or cellular loading
   (pg Fe/cell) too?
3. **Standard spectrum specifics** — the CONRAD default kVp / tube voltage,
   filtration, and energy sampling; whether to model a specific C-arm tube.
4. **Dose & noise model** — model Poisson photon noise (and at what dose/mAs), or
   report the noise-free contrast ceiling? This sets the detectability limit.
5. **Phantom & object scale** — vial/insert sizes, background material
   (water vs blood vs soft tissue), number of inserts, and physical object size
   (mouse-scale in-vivo vs a benchtop vial phantom).
6. **C-arm geometry** — use CONRAD defaults, or specify SID/SDD, detector
   pixel size & count, angular range (short-scan ~200° vs full 360°) and
   reconstruction voxel size?
7. **Detectability criterion** — what counts as "visible" (e.g. Rose criterion
   CNR ≥ 3–5, or a minimum HU difference)?
8. **Detector model** — energy-integrating vs photon-counting; beam-hardening
   correction on/off.

## Repository layout

```
SPEC.md      Full experimental plan and milestones
README.md    This file
DEVLOG.md    Running development log
src/         Simulation & analysis code (to be added)
results/     Generated outputs (git-ignored)
```

## License / attribution

Simulation built on CONRAD (FAU Pattern Recognition Lab). The reference article
PDF is not redistributed (git-ignored).
