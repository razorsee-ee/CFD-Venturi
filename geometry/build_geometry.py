"""Build the gas-side solid of the Venturi scrubber and export it.

Outputs (geometry/out/):
  venturi_fluid.step     fluid volume in mm, for checking in FreeCAD
  inlet.stl, outlet.stl, wall.stl, nozzle.stl
                         surfaces in metres, one file per patch, for snappyHexMesh
"""
import math
from pathlib import Path

import cadquery as cq

import params as p

OUT = Path(__file__).parent / "out"


def shell():
    """Fluid volume of the bare Venturi: revolve the inner profile about the x axis."""
    prof = [(0, 0), (0, p.D_IN / 2), (p.X_CONV, p.D_IN / 2), (p.X_THROAT, p.D_THROAT / 2),
            (p.X_DIFF, p.D_THROAT / 2), (p.X_END, p.D_OUT / 2), (p.X_END, 0)]
    return cq.Workplane("XY").polyline(prof).close().revolve(360, (0, 0, 0), (1, 0, 0))


def nozzle():
    """Lance, nozzle body, struts and deflector disk as one solid."""
    top = p.wall_radius(p.X_STEM) + 50  # extend through the wall so the cut is clean
    stem = (cq.Workplane("XZ", origin=(p.X_STEM, 0, 0)).circle(p.D_STEM / 2)
            .extrude(-(top + p.STEM_BELOW_AXIS)).translate((0, -p.STEM_BELOW_AXIS, 0)))
    body = cq.Workplane("YZ", origin=(p.X_STEM, 0, 0)).circle(p.D_BODY / 2).extrude(p.L_BODY)
    disk = cq.Workplane("YZ", origin=(p.X_DISK_BACK, 0, 0)).circle(p.D_DISK / 2).extrude(p.T_DISK)
    solid = stem.union(body).union(disk)
    for k in range(p.N_STRUTS):
        a = math.radians(90 + 360 * k / p.N_STRUTS)
        a0 = cq.Vector(p.X_NOZZLE - 2, p.R_STRUT_BODY * math.cos(a), p.R_STRUT_BODY * math.sin(a))
        a1 = cq.Vector(p.X_DISK_BACK + 2, p.R_STRUT_DISK * math.cos(a), p.R_STRUT_DISK * math.sin(a))
        d = a1 - a0
        strut = cq.Solid.makeCylinder(p.D_STRUT / 2, d.Length, a0, d.normalized())
        solid = solid.union(cq.Workplane().add(strut))
    return solid


def export_patch(faces, name):
    comp = cq.Compound.makeCompound(faces).scale(1e-3)  # mm -> m
    comp.exportStl(str(OUT / f"{name}.stl"), tolerance=1e-4, angularTolerance=0.1, ascii=True)
    # name the solid so snappyHexMesh / surfaceCheck report the patch name
    path = OUT / f"{name}.stl"
    txt = path.read_text()
    if txt.startswith("solid"):
        lines = txt.splitlines()
        lines[0] = f"solid {name}"
        lines[-1] = f"endsolid {name}"
        path.write_text("\n".join(lines) + "\n")


def main():
    OUT.mkdir(exist_ok=True)
    noz = nozzle()
    fluid = shell().cut(noz)
    cq.exporters.export(fluid, str(OUT / "venturi_fluid.step"))
    cq.exporters.export(noz, str(OUT / "nozzle_solid.step"))

    tol = 1e-3
    inlet, outlet, wall, nozzle_faces = [], [], [], []
    for f in fluid.faces().vals():
        bb = f.BoundingBox()
        if bb.xmax < tol:
            inlet.append(f)
        elif bb.xmin > p.X_END - tol:
            outlet.append(f)
        elif bb.xmin > p.X_STEM - p.D_STEM and bb.xmax < p.X_DISK_FRONT + 1 and bb.ymin < p.D_DISK:
            nozzle_faces.append(f)
        else:
            wall.append(f)
    for faces, name in ((inlet, "inlet"), (outlet, "outlet"), (wall, "wall"), (nozzle_faces, "nozzle")):
        export_patch(faces, name)
        print(f"{name:7s}: {len(faces)} faces")
    print(f"fluid volume: {fluid.val().Volume() * 1e-9:.4f} m3")


if __name__ == "__main__":
    main()
