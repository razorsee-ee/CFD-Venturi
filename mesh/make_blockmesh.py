"""Write <case>/system/blockMeshDict: structured O-grid that follows the Venturi wall.

Five blocks per axial segment (central square + four wall blocks). The wall
patch is exact (cones between stations), cells are graded towards the wall
for wall functions and refined axially around the nozzle and the throat.
The nozzle itself is cut out afterwards by snappyHexMesh.

Usage:  python3 make_blockmesh.py [scale] [--zone full|throat]
        scale < 1 coarsens the mesh (e.g. 0.5 for a quick test), > 1 refines it.
        --zone full   (default) whole scrubber + 2 m outlet extension -> case/
        --zone throat nozzle, end of the convergent, throat and 0.6 m of the
                      diffuser + 0.3 m straight extension -> case-throat/
        --zone rim    same domain -> case-rim/ (spray from the disk rim, breakup)
"""
import argparse
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "geometry"))
import params as p  # noqa: E402

_ap = argparse.ArgumentParser()
_ap.add_argument("scale", nargs="?", type=float, default=1.0)
_ap.add_argument("--zone", choices=("full", "throat", "rim"), default="full")
ARGS = _ap.parse_args()
SCALE = ARGS.scale
N_T = 24                # cells per quarter circumference (and across the central square)
N_R = 22                # cells from the central square to the wall
GRADE_R = 8.0           # radial expansion ratio, wall -> centre (first cell ~3 mm at the inlet)
SQUARE = 0.55           # central square corner at SQUARE*R on the diagonal
BULGE = 0.45            # central square edge mid-point at BULGE*R

X_REF0 = p.X_STEM - 150.0  # start of the fine axial zone around the nozzle

# per zone: case folder, x where the real wall ends (constant radius after it:
# numerical outlet extension) and axial stations (mm) with the target axial
# cell size there (mm)
ZONES = {
    "full": dict(
        case="case",
        wall_end=p.X_END,
        stations=[
            (0.0, 60.0),            # coarse at the inlet ...
            (p.X_CONV, 20.0),
            (X_REF0, 6.0),
            (p.X_THROAT, 6.0),
            (p.X_DIFF, 8.0),
            (p.X_END, 30.0),
            (p.X_END + 2000.0, 100.0),  # ... and on the 2 m outlet extension
        ]),
    "throat": dict(
        case="case-throat",
        wall_end=p.X_DIFF + 600.0,
        stations=[
            (1750.0, 12.0),         # inlet: profiles mapped from the full case
            (X_REF0, 4.0),
            (p.X_THROAT, 4.0),
            (p.X_DIFF, 5.0),
            (p.X_DIFF + 600.0, 12.0),
            (p.X_DIFF + 900.0, 25.0),   # 0.3 m straight extension
        ]),
}
ZONES["rim"] = dict(ZONES["throat"], case="case-rim")   # same domain, spray from the disk rim
ZONE = ZONES[ARGS.zone]
STATIONS = ZONE["stations"]


def radius(x):
    return p.wall_radius(min(x, ZONE["wall_end"]))


def axial_cells(x0, d0, x1, d1):
    """Number of cells and expansion ratio for a geometric spacing d0 -> d1 over x0..x1."""
    L = x1 - x0
    d0, d1 = d0 / SCALE, d1 / SCALE
    if abs(d1 - d0) < 1e-9:
        return max(1, round(L / d0)), 1.0
    n = L * math.log(d1 / d0) / (d1 - d0)
    return max(1, round(n)), d1 / d0


def main():
    nt = max(4, round(N_T * SCALE))
    nr = max(4, round(N_R * SCALE))
    mm = 1e-3
    verts, edges, blocks = [], [], []

    # 8 vertices per station: inner square i0..i3, outer o0..o3 (CCW in the y-z plane)
    for x, _ in STATIONS:
        R = radius(x)
        s = SQUARE * R / math.sqrt(2)
        for (cy, cz) in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
            verts.append((x, cy * s, cz * s))
        for ang in (225, 315, 45, 135):
            a = math.radians(ang)
            verts.append((x, R * math.cos(a), R * math.sin(a)))
        k = 8 * (len(verts) // 8 - 1)
        # arcs: outer circle and slightly convex inner square, mid-points at 270, 0, 90, 180 deg
        for j, ang in enumerate((270, 0, 90, 180)):
            a = math.radians(ang)
            edges.append((k + 4 + j, k + 4 + (j + 1) % 4, (x, R * math.cos(a), R * math.sin(a))))
            edges.append((k + j, k + (j + 1) % 4, (x, BULGE * R * math.cos(a), BULGE * R * math.sin(a))))

    for st in range(len(STATIONS) - 1):
        (x0, d0), (x1, d1) = STATIONS[st], STATIONS[st + 1]
        nx, gx = axial_cells(x0, d0, x1, d1)
        a, b = 8 * st, 8 * (st + 1)
        i = [a + j for j in range(4)]
        o = [a + 4 + j for j in range(4)]
        blocks.append((i + [b + j for j in range(4)], (nt, nt, nx), (1, 1, gx)))
        for j in range(4):
            j1 = (j + 1) % 4
            bot = [o[j], o[j1], i[j1], i[j]]
            top = [v + 8 for v in bot]
            blocks.append((bot + top, (nt, nr, nx), (1, GRADE_R, gx)))

    nst = len(STATIONS)
    inlet = [[0, 3, 2, 1]] + [[4 + j, j, (j + 1) % 4, 4 + (j + 1) % 4] for j in range(4)]
    last = 8 * (nst - 1)
    outlet = [[last + j for j in range(4)]] + \
        [[last + 4 + j, last + 4 + (j + 1) % 4, last + (j + 1) % 4, last + j] for j in range(4)]
    wall = [[8 * s + 4 + j, 8 * (s + 1) + 4 + j, 8 * (s + 1) + 4 + (j + 1) % 4, 8 * s + 4 + (j + 1) % 4]
            for s in range(nst - 1) for j in range(4)]

    f = lambda v: f"({v[0] * mm:.6g} {v[1] * mm:.6g} {v[2] * mm:.6g})"  # noqa: E731
    out = [HEADER, f"// generated by mesh/make_blockmesh.py, scale = {SCALE}, zone = {ARGS.zone}",
           f"// domain (mm): x0 = {STATIONS[0][0]:g}, wall_end = {ZONE['wall_end']:g}, x_end = {STATIONS[-1][0]:g}\n",
           "convertToMeters 1;\n", "vertices\n("]
    out += [f"    {f(v)} // {n}" for n, v in enumerate(verts)]
    out += [");\n", "edges\n("]
    out += [f"    arc {e[0]} {e[1]} {f(e[2])}" for e in edges]
    out += [");\n", "blocks\n("]
    ncells = 0
    for vs, (n1, n2, n3), (g1, g2, g3) in blocks:
        out.append(f"    hex ({' '.join(map(str, vs))}) ({n1} {n2} {n3}) simpleGrading ({g1:g} {g2:g} {g3:.5g})")
        ncells += n1 * n2 * n3
    out += [");\n", "boundary\n("]
    for name, typ, faces in (("inlet", "patch", inlet), ("outlet", "patch", outlet), ("walls", "wall", wall)):
        out += [f"    {name}\n    {{\n        type {typ};\n        faces\n        ("]
        out += [f"            ({' '.join(map(str, q))})" for q in faces]
        out += ["        );\n    }"]
    out += [");\n", "// " + "*" * 73 + " //\n"]

    target = ROOT / ZONE["case"] / "system" / "blockMeshDict"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("\n".join(out))
    print(f"{target}: {len(blocks)} blocks, {ncells} cells (before snappyHexMesh)")


HEADER = r"""/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  12
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""

if __name__ == "__main__":
    main()
