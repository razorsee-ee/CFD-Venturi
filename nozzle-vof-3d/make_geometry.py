"""STL surfaces for snappyHexMesh: the lance (with the oil orifice on its end face) and the disk.

lance.stl  regions: nozzle (side, back cap, ring around the orifice), oilInlet (orifice
           disc r < 10 mm on the end face at x = 0.082, 1 cm in front of the disk)
disk.stl   region:  disk (D65 x 8 mm, x 0.092...0.100)
Usage: python make_geometry.py
"""
import math
from pathlib import Path

N = 720                                  # segments around the axis
X_BACK, X_TIP = -0.06, 0.082             # lance: starts outside the domain (x = -0.05)
R_ORIF, R_LANCE = 0.010, 0.014
DISK_X, DISK_R = (0.092, 0.100), 0.0325
OUT = Path(__file__).resolve().parent / "constant" / "geometry"


def ring(x, r):
    return [(x, r * math.cos(2 * math.pi * k / N), r * math.sin(2 * math.pi * k / N)) for k in range(N)]


def facet(f, a, b, c):
    ux, uy, uz = (b[i] - a[i] for i in range(3))
    vx, vy, vz = (c[i] - a[i] for i in range(3))
    n = (uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx)
    m = math.sqrt(sum(q * q for q in n)) or 1.0
    f.append(f"  facet normal {n[0] / m:.6e} {n[1] / m:.6e} {n[2] / m:.6e}\n    outer loop\n"
             + "".join(f"      vertex {p[0]:.9e} {p[1]:.9e} {p[2]:.9e}\n" for p in (a, b, c))
             + "    endloop\n  endfacet\n")


def side(f, x0, x1, r):                  # outward normals
    a, b = ring(x0, r), ring(x1, r)
    for k in range(N):
        k1 = (k + 1) % N
        facet(f, a[k], a[k1], b[k1])
        facet(f, a[k], b[k1], b[k])


def disc(f, x, r, outward):              # outward = +1: normal +x
    c, p = (x, 0.0, 0.0), ring(x, r)
    for k in range(N):
        k1 = (k + 1) % N
        facet(f, c, p[k], p[k1]) if outward > 0 else facet(f, c, p[k1], p[k])


def annulus(f, x, r0, r1, outward):
    a, b = ring(x, r0), ring(x, r1)
    for k in range(N):
        k1 = (k + 1) % N
        if outward > 0:
            facet(f, a[k], b[k], b[k1]); facet(f, a[k], b[k1], a[k1])
        else:
            facet(f, a[k], b[k1], b[k]); facet(f, a[k], a[k1], b[k1])


def solid(name, facets):
    return f"solid {name}\n" + "".join(facets) + f"endsolid {name}\n"


OUT.mkdir(parents=True, exist_ok=True)
noz, oil = [], []
side(noz, X_BACK, X_TIP, R_LANCE)
disc(noz, X_BACK, R_LANCE, -1)
annulus(noz, X_TIP, R_ORIF, R_LANCE, +1)
disc(oil, X_TIP, R_ORIF, +1)
(OUT / "lance.stl").write_text(solid("nozzle", noz) + solid("oilInlet", oil))
d = []
side(d, DISK_X[0], DISK_X[1], DISK_R)
disc(d, DISK_X[0], DISK_R, -1)
disc(d, DISK_X[1], DISK_R, +1)
(OUT / "disk.stl").write_text(solid("disk", d))
print("written constant/geometry/lance.stl, disk.stl")
