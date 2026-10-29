"""Turn a cut plane of the full-scrubber solution into inlet data for case-throat.

Input : a 'foam'-format cut plane written by the cutPlaneSurface function
        (faceCentres, vectorField/U, scalarField/{T,k,omega})
Output: <case>/constant/inletProfile/points and .../0/{T,k,omega}: mapped as they are
        (timeVaryingMappedFixedValue; Allrun copies them to constant/boundaryData/inlet);
        <case>/constant/inletSwirl: tangential velocity W(r), circumferential average,
        for the swirlFlowRateInletVelocity inlet. That keeps the exact mass flow
        whatever the inlet pressure (a mapped velocity would not: the inlet
        density rises with the spray pressure drop).

Usage:  python make_inlet_data.py <cutPlane dir> <case dir>
"""
import math
import re
import sys
from pathlib import Path

import numpy as np

# Windows consoles/files default to cp1252: the reports contain Cyrillic and m²
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "geometry"))
import params as geo  # noqa: E402

N_RING = 128   # extra points on the wall circle: U = 0 there, so the boundary layer is not smeared to the wall

src, case = Path(sys.argv[1]), Path(sys.argv[2])
dst = case / "constant" / "inletProfile"
(dst / "0").mkdir(parents=True, exist_ok=True)


def body(path):
    """Text of an OpenFOAM list file without any header before the size line."""
    txt = path.read_text(encoding="utf-8", errors="ignore")
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
    path.write_text(f"// {header}\n{len(rows)}\n(\n" + "\n".join(rows) + "\n)\n", encoding="utf-8")


write(dst / "points", "cut-plane face centres + wall ring, from tools/make_inlet_data.py",
      [f"({p[0]:.6g} {p[1]:.6g} {p[2]:.6g})" for p in np.r_[pts, ring]])
U = np.array(re.findall(VEC, body(src / "vectorField" / "U")), float)
(dst / "0" / "U").unlink(missing_ok=True)
for name in ("T", "k", "omega"):
    v = np.array(re.findall(r"^\s*([-\d.eE+]+)\s*$", body(src / "scalarField" / name), re.M)[1:], float)
    write(dst / "0" / name, f"{name}: cut plane + nearest values on the wall ring",
          [f"{x:.6g}" for x in np.r_[v, v[near]]])

# tangential velocity about +x, circumferential average in radial bins, 0 at the wall
r = np.hypot(pts[:, 1], pts[:, 2])
ut = (pts[:, 1] * U[:, 2] - pts[:, 2] * U[:, 1]) / np.maximum(r, 1e-9)
edges = np.linspace(0, R, 21)
rows = ["(0 0)"]
for a, b in zip(edges[:-1], edges[1:]):
    sel = (r >= a) & (r < b)
    if sel.any():
        rows.append(f"({r[sel].mean():.5g} {ut[sel].mean():.5g})")
rows.append(f"({R:.5g} 0)")
(case / "constant" / "inletSwirl").write_text(
    "// tangential velocity W(r) [m/s] at the inlet plane, from tools/make_inlet_data.py\n"
    "tangentialVelocity\n{\n    type    function1;\n    value2  table\n    (\n        "
    + "\n        ".join(rows) + "\n    );\n}\n", encoding="utf-8")

print(f"{len(pts)} points at x = {x0:.3f} m, radius up to {r.max():.3f} m, + {N_RING} wall points at r = {R:.3f} m")
print(f"axial velocity {U[:, 0].mean():.1f} m/s (mean of faces), tangential velocity "
      f"{np.abs(ut).max():.1f} m/s (max)")
print(f"written {dst}")
