# RabbitCT — C-arm Cone-Beam CT Benchmark Dataset (post-mortem rabbit)

This record re-publishes the **RabbitCT** benchmark dataset, an open platform for
comparing and ranking 3-D cone-beam CT **backprojection** implementations on a
single high-resolution C-arm CT scan of a **post-mortem rabbit**. The original
challenge site is no longer online; this archival copy is released (with the
authors' agreement) under **CC BY 4.0** to keep the benchmark available.

## Please cite

Rohkohl C, Keck B, Hofmann HG, Hornegger J. *Technical Note: RabbitCT — an open
platform for benchmarking 3D cone-beam reconstruction algorithms.* **Medical
Physics** 36(9):3940–3944, 2009. DOI: 10.1118/1.3180956

## Ethics

The scan is a **post-mortem** acquisition; no live-animal procedure was performed.
See the Medical Physics paper above for details.

## Contents

Both the **original (2011)** and the **corrected `-v2` (2012)** versions of the
512³ and 1024³ datasets are included for a complete historical record.

| File | Size | Description |
|------|------|-------------|
| `rabbitct_512-v2.rctd` / `rabbitct_512.rctd`   | ~2.5 GB | 496 projections + calibrated 3×4 projection matrices; 512³ reconstruction target (voxel 0.5 mm, 256 mm³ volume). `-v2` = corrected 2012 release. |
| `rabbitct_1024-v2.rctd` / `rabbitct_1024.rctd` | ~4.5 GB | Same scan at the 1024³ target (voxel 0.25 mm). `-v2` = corrected 2012 release. |
| `reference_256.vol`  | 64 MB   | Reference reconstruction, 256³ float32 (256 mm³, 1.0 mm voxels). |
| `rabbitct_develop-v2.zip` / `rabbitct_develop.zip` | small | Developer kit: `rabbitct.h` (data struct), CMake project, `LolaBunny` example backprojection module. |

## `.rctd` format (summary)

Binary container per the developer kit (`rabbitct.h`, struct `RabbitCtGlobalData`):
projection image dimensions `S_x`×`S_y`, per-projection **3×4 projection matrices**
`A_n` (world→detector, in mm), the projection images `I_n`, and the isotropic
reconstruction voxel size `R_L`. A backprojection module implements
`RCTAlgorithmBackprojection` and is driven by the RabbitCTRunner. The reconstruction
volume is a centered 256 mm cube; voxel size = 256 mm / N (N = 256/512/1024).

## Geometry

Circular short-scan C-arm cone-beam trajectory, 496 projections, calibrated per
projection via the `A_n` matrices (no nominal SID/SDD needed — the matrices are
the geometry). CONRAD provides a compatible backprojector under the codename
`LolaBunny` (`edu.stanford.rsl.conrad.reconstruction.LolaBunnyBackprojector`).

## License

Creative Commons Attribution 4.0 International (CC BY 4.0). Attribution: cite the
Medical Physics paper above.
