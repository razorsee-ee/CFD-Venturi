"""Turn a cut plane of the full-scrubber solution into inlet data for case-throat.

Input : a 'foam'-format cut plane written by the cutPlaneSurface function
        (faceCentres, vectorField/U, scalarField/{T,k,omega})
Output: <case>/constant/inletProfile/points and .../0/{U,T,k,omega}. Allrun copies
        it to constant/boundaryData/inlet (read by timeVaryingMappedFixedValue)
        and scales U to the exact mass flow (tools/scale_inlet_flow.py).

Usage:  python make_inlet_data.py <cutPlane dir> <case dir>
"""
import math
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "geometry"))
import params as geo  # noqa: E402

N_RING = 128   # extra points on the wall circle: U = 0 there, so the boundary layer is not smeared to the wall

src, case = Path(sys.argv[1]), Path(sys.argv[2])
dst = case / "constant" / "inletProfile"
(dst / "0").mkdir(parents=True, exist_ok=True)


def body(path):
    """Text of an OpenFOAM list file without any header before the size line."""
    txt = path.read_text()
    m = re.search(r"^\s*(\d+)\s*\n\s*\(", txt, re.M)
    return txt[m.start():]


VEC = r"\(([-\d.eE+]+) ([-\d.eE+]+) ([-\d.eE+]+)\)"
pts = np.array(re.findall(VEC, body(src / "faceCentres")), float)
x0 = pts[0, 0]
R = geo.wall_radius(x0 * 1e3) * 1e-3
ang = np.linspace(0, 2 * math.pi, N_RING, endpoint=False)
ring = np.c_[np.full(N_RING, x0), R * np.cos(ang), R * np.sin(ang)]
r_src = np.hypot(pts[:, 1], pts[:, 2])
near = [int(np.argmin((pts[:, 1] - p[1])**2 + (pts[:, 2] - p[2])**2)) for p in ring]


def write(path, header, rows):
    path.write_text(f"// {header}\n{len(rows)}\n(\n" + "\n".join(rows) + "\n)\n")


write(dst / "points", "cut-plane face centres + wall ring, from tools/make_inlet_data.py",
      [f"({p[0]:.6g} {p[1]:.6g} {p[2]:.6g})" for p in np.r_[pts, ring]])
U = np.array(re.findall(VEC, body(src / "vectorField" / "U")), float)
write(dst / "0" / "U", "U: cut plane + 0 on the wall ring",
      [f"({u[0]:.6g} {u[1]:.6g} {u[2]:.6g})" for u in np.r_[U, np.zeros((N_RING, 3))]])
for name in ("T", "k", "omega"):
    v = np.array(re.findall(r"^\s*([-\d.eE+]+)\s*$", body(src / "scalarField" / name), re.M)[1:], float)
    write(dst / "0" / name, f"{name}: cut plane + nearest values on the wall ring",
          [f"{x:.6g}" for x in np.r_[v, v[near]]])

# report what the inlet carries
r = np.hypot(pts[:, 1], pts[:, 2])
ut = (pts[:, 1] * U[:, 2] - pts[:, 2] * U[:, 1]) / np.maximum(r, 1e-9)
print(f"{len(pts)} points at x = {x0:.3f} m, radius up to {r.max():.3f} m, + {N_RING} wall points at r = {R:.3f} m")
print(f"axial velocity {U[:, 0].mean():.1f} m/s (mean of faces), tangential velocity "
      f"{np.abs(ut).max():.1f} m/s (max)")
print(f"written {dst}")
