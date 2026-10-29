"""Copy constant/inletProfile -> constant/boundaryData/inlet and scale U to a mass flow.

Interpolating the mapped profile onto another mesh changes the mass flow by a
few per cent (the boundary layer at the wall). This measures the inlet mass
flow the mapped profile gives on this mesh (foamPostProcess, time 0) and scales
the velocity data so that it is exactly the target.

Usage (case directory, after the mesh and 0/ exist):
    python ../tools/scale_inlet_flow.py [target kg/s, default 21.11]
"""
import re
import shutil
import subprocess
import sys
from pathlib import Path

target = float(sys.argv[1]) if len(sys.argv) > 1 else 21.11
case = Path(".")
src, dst = case / "constant" / "inletProfile", case / "constant" / "boundaryData" / "inlet"


def mass_flow():
    out = subprocess.run(["foamPostProcess", "-solver", "fluid", "-time", "0",
                          "-func", "patchFlowRate(patch=inlet)"],
                         capture_output=True, text=True).stdout
    m = re.findall(r'sum\("?inlet"?\) of phi = ([-\d.eE+]+)', out)
    if not m:
        sys.exit("scale_inlet_flow: could not read the inlet mass flow:\n" + out[-2000:])
    return abs(float(m[-1]))


shutil.rmtree(dst, ignore_errors=True)
shutil.copytree(src, dst)
m0 = mass_flow()
f = target / m0
path = dst / "0" / "U"
txt = path.read_text()
txt = re.sub(r"\(([-\d.eE+]+) ([-\d.eE+]+) ([-\d.eE+]+)\)",
             lambda m: "({:.6g} {:.6g} {:.6g})".format(*(float(v) * f for v in m.groups())), txt)
path.write_text(txt)
m1 = mass_flow()
print(f"inlet mass flow: mapped profile {m0:.3f} kg/s -> scaled by {f:.4f} -> {m1:.3f} kg/s (target {target})")
