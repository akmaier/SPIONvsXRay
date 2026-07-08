# CONRAD OpenCL contributions — re-release via pyconrad

Four OpenCL contributions were upstreamed to **akmaier/CONRAD** `master`
(commits `71bb20c`, `631fbed`, `fdf5365`) and baked into a patched `conrad_1.1.0.jar`:

1. **GPU pointwise `exp`/`log`** — `OpenCLGridOperators.exp()` called a kernel
   `"expontial"` that did not exist (no natural-exp kernel at all) and `log()`
   called a missing `"logarithm"`; both threw `CLInvalidKernelNameException`.
   Added `exponential` + `logarithm` kernels to `PointwiseOperators.cl` and fixed
   the `exp()` call site.
2. **Fan backprojector pixel spacing** — `FanBeamBackprojector2D`'s CL kernel never
   received the reconstruction pixel spacing (`// TODO: Spacing`) and used
   placeholder detector geometry, so it only reconstructed correctly at 1 mm/px.
   Added `setSpacing()`, passed spacing to the kernel, corrected the 2D detector
   geometry (matches CPU to ~1e-5, ~850× faster).
3. **GPU noise generators** — CONRAD had no GPU RNG. Added a Philox4x32-10
   counter-based RNG in `PointwiseOperators.cl` with `poisson` (Knuth / Hörmann
   PTRS, exact) and `standardNormal` (Box–Muller) kernels; CPU reference impls in
   `NumericGridOperator`; `NumericPointwiseOperators.poisson/standardNormal(grid,seed)`.
4. **GPU spectral detectors** — `OpenCLEnergyIntegratingDetector` and
   `OpenCLPhotonCountingDetector` (off a shared `OpenCLSpectralDetector` + fused
   `SpectralDetector.cl` kernels): integrate the Beer–Lambert law over energy on a
   per-material path-length stack, with **per-energy** Poisson noise (correct
   energy-integrated variance ΣE²·n, not Poisson(ΣE·n)).

Verified against the patched jar with the `conrad_ext` shadow disabled and OpenCL
active: `exp(log(x))=x` to 1e-7; 0.5 mm recon with monotonic iron ΔHU (c20 → +7.14);
GPU Poisson mean==var==λ (0.5–5000) with chi-square GoF; EID/PCD match a
CONRAD-sourced reference to 1e-7 / exactly, and the noisy EID variance equals ΣE²·n.

## Artifacts (git-ignored — large binaries)
- `publish/conrad/conrad_1.1.0.jar` — patched jar (sha1 `6d41971c3f5d…`; original was `38fe24e116da…`).
- `publish/conrad/CONRAD_1.1.0.zip` — the canonical FAU release zip with **only**
  `CONRAD 1.1.0/conrad_1.1.0.jar` swapped for the patched jar (structure preserved,
  `__MACOSX` cruft stripped). This is the file to host.

The jar was produced by recompiling the changed classes against the released jar
and splicing the new `.class` + `.cl` resources in with `jar uf` (no full
1205-file rebuild). To reproduce: `scripts/rebuild_conrad_jar.sh` (see below) or
Eclipse *Export → Runnable JAR*.

## To publish (requires FAU access — Claude cannot reach the fileadmin)
pyconrad downloads the jar from a zip on the FAU web server
(`pyconrad/_download_conrad.py`):

```
https://www5.cs.fau.de/fileadmin/user_upload/CONRAD_1.1.0.zip
```

**Same-version re-host (simplest):** overwrite that URL with
`publish/conrad/CONRAD_1.1.0.zip`. Any fresh `pip install pyconrad` (or a call to
`pyconrad._download_conrad.download_conrad()` into a clean dir) then pulls the
fixed jar. Existing installs keep their cached jar until re-downloaded.

**Versioned release (cleaner, optional):** host as `CONRAD_1.1.1.zip`, bump the
URL in `pyconrad/_download_conrad.py`, and release pyconrad to PyPI.

## After the fixed jar is live
The local `conrad_ext/` shadow (patched `FanBeamBackprojector2D`) becomes
redundant — the fix is in the jar. Until the re-hosted jar is picked up, keep
`conrad_ext` active (our env still uses the currently-hosted, unpatched jar).
Retire it by removing the `CONRAD_DEV_DIRS` shadow in `src/conrad_backend.py`
once `pyconrad` ships the patched jar.
