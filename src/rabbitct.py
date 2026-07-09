"""RabbitCT dataset reader — projection geometry + reference volume.

The RabbitCT benchmark ships a binary `.rctd` file. There is no public parser in
the distributed sources (the loader lived in the challenge's C++ framework), so we
reverse-engineered the layout from the header struct (data/rabbitct/develop/.../
rabbitct.h) and the binary:

  bytes  0.. 3  uint32  version (=2)
  bytes  4.. 7  uint32  S_x  (projection width  = 1248)
  bytes  8..11  uint32  S_y  (projection height = 960)
  bytes 12..15  uint32  numProjections (= 496)
  bytes 16..19  float32 R_L  (isotropic voxel size, mm)
  bytes 20..23  float32 O_L  (0-index world position)
  then per projection:
     96 bytes   3x4 double projection matrix A_n, COLUMN-major
     S_x*S_y*4  float32 projection image I_n

The matrix maps world[mm] homogeneous -> detector[px] homogeneous. Decomposing it
recovers a real Siemens C-arm short scan: SID ~ 745 mm, SDD ~ 1196 mm (0.308 mm
pixel pitch), 198 deg over 496 views, ~0.40 deg/step, rotation about z.

Cite: Rohkohl, Keck, Hofmann, Hornegger, "Technical Note: RabbitCT ...",
Med. Phys. 36(9):3940-3944 (2009).
"""
from __future__ import annotations
import struct
import numpy as np

HEADER = 24
PIXEL_PITCH_MM = 0.308      # RabbitCT flat detector pitch


def read_header(path):
    with open(path, "rb") as f:
        h = f.read(HEADER)
    version, S_x, S_y, nP = struct.unpack("<4I", h[:16])
    R_L, O_L = struct.unpack("<2f", h[16:24])
    return dict(version=version, S_x=S_x, S_y=S_y, num_proj=nP, R_L=R_L, O_L=O_L)


def read_matrices(path):
    """Return the (num_proj, 3, 4) world[mm]->detector[px] projection matrices."""
    hd = read_header(path)
    S_x, S_y, nP = hd["S_x"], hd["S_y"], hd["num_proj"]
    stride = 96 + S_x * S_y * 4
    mats = np.zeros((nP, 3, 4))
    for k in range(nP):
        d = np.fromfile(path, dtype="<f8", count=12, offset=HEADER + k * stride)
        mats[k] = d.reshape(3, 4, order="F")
    return mats, hd


def _rq(M):
    P = np.flipud(M)
    Q, R = np.linalg.qr(P.T)
    R = np.flipud(R.T); Q = Q.T
    R = R[:, ::-1]; Q = Q[::-1, :]
    D = np.diag(np.sign(np.diag(R)))
    return R @ D, D @ Q


def geometry(path):
    """Decompose the matrices into a C-arm geometry summary."""
    mats, hd = read_matrices(path)
    sids, angs, cz = [], [], []
    for P in mats:
        M = P[:, :3]
        C = -np.linalg.solve(M, P[:, 3])
        sids.append(np.linalg.norm(C))
        angs.append(np.degrees(np.arctan2(C[1], C[0])))
        cz.append(C[2])
    angs = np.degrees(np.unwrap(np.radians(angs)))
    K, _ = _rq(mats[0][:, :3]); K = K / K[2, 2]
    f_px = 0.5 * (abs(K[0, 0]) + abs(K[1, 1]))
    return dict(num_proj=len(mats),
                detector_px=(hd["S_x"], hd["S_y"]),
                pixel_pitch_mm=PIXEL_PITCH_MM,
                voxel_mm=hd["R_L"],
                sid_mm=float(np.mean(sids)), sid_std=float(np.std(sids)),
                sdd_mm=float(f_px * PIXEL_PITCH_MM),
                focal_px=float(f_px),
                principal_point_px=(float(K[0, 2]), float(K[1, 2])),
                angle_start_deg=float(angs[0]), angle_end_deg=float(angs[-1]),
                angle_range_deg=float(angs[-1] - angs[0]),
                angle_step_deg=float(np.mean(np.diff(angs))),
                source_z_mm=float(np.mean(cz)))


def load_reference_volume(path, L=256):
    """Load the reference reconstruction volume (L^3 float32)."""
    vol = np.fromfile(path, dtype="<f4", count=L * L * L)
    return vol.reshape(L, L, L)


if __name__ == "__main__":
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import conrad_backend

    root = conrad_backend.REPO_ROOT
    rctd = str(root / "data" / "rabbitct" / "rabbitct_512-v2.rctd")
    volf = str(root / "data" / "rabbitct" / "reference_256.vol")

    g = geometry(rctd)
    print("RabbitCT C-arm geometry:")
    for k, v in g.items():
        print(f"  {k}: {v}")

    outdir = str(root / "results" / "rabbit")
    os.makedirs(outdir, exist_ok=True)
    if os.path.isfile(volf):
        vol = load_reference_volume(volf, 256)
        print(f"reference volume: {vol.shape}, range [{vol.min():.1f},{vol.max():.1f}]")
        fig, ax = plt.subplots(1, 3, figsize=(11, 3.8))
        for a, (idx, name) in zip(ax, [(128, "axial"), (128, "coronal"), (128, "sagittal")]):
            sl = {"axial": vol[idx], "coronal": vol[:, idx], "sagittal": vol[:, :, idx]}[name]
            lo, hi = np.percentile(sl, [40, 99.5])
            a.imshow(sl, cmap="gray", vmin=lo, vmax=hi); a.set_title(f"rabbit {name}"); a.axis("off")
        fig.tight_layout(); fig.savefig(f"{outdir}/rabbit_slices.png", dpi=130); plt.close(fig)
        print("wrote", f"{outdir}/rabbit_slices.png")
