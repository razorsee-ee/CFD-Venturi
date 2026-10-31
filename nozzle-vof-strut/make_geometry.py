"""STL surfaces for snappyHexMesh: lance (oil exit on its end face), disk, one strut.

lance.stl  regions: nozzle (side, back cap, ring around the exit), oilInlet (exit disc
           r < 11.2 mm on the end face at x = 0.082; the mesh starts at the stem r = 5 mm,
           so the oil leaves through the annulus D22.4/10 = the area of the D20 orifice)
disk.stl   region:  disk (D65 x 8 mm, x 0.092...0.100)
stem.stl   region:  stem - D10 rod from inside the lance to inside the disk (the disk is
           carried by it; behind the disk there is nothing)
strut.stl  region:  strut - three flat bars from the lance to the disk rim at 0, 120 and
           240 deg (the first in the plane z = 0, y > 0); ends sunk into the lance and the disk
Usage: python make_geometry.py
"""
import math
from pathlib import Path

N = 720                                  # segments around the axis
X_BACK, X_TIP = -0.06, 0.082             # lance: starts outside the domain (x = -0.05)
R_ORIF, R_LANCE = 0.0112, 0.014
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
sm = []
side(sm, X_TIP - 0.01, DISK_X[0] + 0.004, 0.005)
disc(sm, X_TIP - 0.01, 0.005, -1)
disc(sm, DISK_X[0] + 0.004, 0.005, +1)
(OUT / "stem.stl").write_text(solid("stem", sm))

# strut: flat bar, centre line A -> B in the (x, r = y) plane at z = 0
A, B = (0.070, 0.016), (0.098, 0.032)           # (x, r): starts on the lance, ends inside the disk
THICK, WIDTH = 0.005, 0.012                      # across the bar in the x-r plane, around (z)
L = math.hypot(B[0] - A[0], B[1] - A[1])
e1 = ((B[0] - A[0]) / L, (B[1] - A[1]) / L, 0.0)
e2 = (-e1[1], e1[0], 0.0)
e3 = (0.0, 0.0, 1.0)


def corner(s1, s2, s3):
    base = (A[0], A[1], 0.0)
    return tuple(base[i] + s1 * L * e1[i] + s2 * THICK / 2 * e2[i] + s3 * WIDTH / 2 * e3[i] for i in range(3))


c = {(a, b, k): corner(a, b, k) for a in (0, 1) for b in (-1, 1) for k in (-1, 1)}
quads = [((0, -1, -1), (0, 1, -1), (0, 1, 1), (0, -1, 1)),       # end at A   (normal -e1)
         ((1, -1, -1), (1, -1, 1), (1, 1, 1), (1, 1, -1)),       # end at B   (+e1)
         ((0, -1, -1), (0, -1, 1), (1, -1, 1), (1, -1, -1)),     # -e2 side
         ((0, 1, -1), (1, 1, -1), (1, 1, 1), (0, 1, 1)),         # +e2 side
         ((0, -1, -1), (1, -1, -1), (1, 1, -1), (0, 1, -1)),     # -z side
         ((0, -1, 1), (0, 1, 1), (1, 1, 1), (1, -1, 1))]         # +z side
st = []
for k in range(3):                               # three struts, 120 deg apart
    ca, sa = math.cos(2 * math.pi * k / 3), math.sin(2 * math.pi * k / 3)
    rot = {key: (p[0], p[1] * ca - p[2] * sa, p[1] * sa + p[2] * ca) for key, p in c.items()}
    for q in quads:
        p = [rot[key] for key in q]
        facet(st, p[0], p[2], p[1])             # outward normals
        facet(st, p[0], p[3], p[2])
(OUT / "strut.stl").write_text(solid("strut", st))
print("written constant/geometry/lance.stl, disk.stl, stem.stl, strut.stl")
