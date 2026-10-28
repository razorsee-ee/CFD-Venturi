"""Internal (wetted) geometry of the Venturi scrubber, taken from the drawing.

All dimensions in mm. x = axis along gas flow, x = 0 at the inlet end of the
cylindrical section, y = up (towards the nozzle stem penetration).
Only the gas-side geometry is kept: manhole M1, nozzle N1, lugs, flanges,
insulation and the 100 mm site-adjustment allowance are ignored where they do
not change the flow path.
"""

WALL = 15.0                     # plate thickness used on the drawing

# --- shell (inner diameters = O.D. - 2*WALL) ---
D_IN = 914.0 - 2 * WALL         # 884: cylindrical inlet section (pos. 8)
D_THROAT = 560.0                # throat I.D. (pos. 9, section D-D)
D_OUT = 1010.0 - 2 * WALL       # ~980: diffuser outlet; O.D. measured off the drawing (~1010), not dimensioned

L_CYL = 240.0 + 1104.0          # 1344: inlet end -> start of convergent (incl. 100 mm site allowance)
L_CONV = 2285.0 - 1104.0        # 1181: convergent (pos. 6)
L_THROAT = 248.0                # throat
L_DIFF = 2134.0                 # diffuser (pos. 5 + 4)

X_CONV = L_CYL
X_THROAT = X_CONV + L_CONV
X_DIFF = X_THROAT + L_THROAT
X_END = X_DIFF + L_DIFF         # 4907

# --- spray nozzle (pos. 1, section D-D) ---
X_STEM = X_DIFF - (90.0 + 438.0)  # stem axis: 528 mm upstream of throat end -> 280 mm before throat
D_STEM = 48.3                   # vertical lance pipe, enters through the top of the convergent
STEM_BELOW_AXIS = 40.0          # lance end below the axis
D_BODY = 28.0                   # nozzle body on the axis, pointing downstream
L_BODY = 69.0                   # stem axis -> nozzle outlet face
GAP = 100.0 - 8.0               # nozzle outlet face -> back face of the deflector disk
D_DISK = 65.0                   # deflector disk
T_DISK = 8.0
N_STRUTS = 3                    # struts nozzle -> disk, 120 deg apart
D_STRUT = 6.0                   # assumed (not on the drawing)
R_STRUT_BODY = 10.0             # strut attachment radius on the nozzle face (assumed)
R_STRUT_DISK = 26.0             # strut attachment radius on the disk (assumed)
D_ORIFICE = 12.0                # water orifice, used for the injector model only (assumed)

X_NOZZLE = X_STEM + L_BODY      # 2314
X_DISK_BACK = X_NOZZLE + GAP    # 2406
X_DISK_FRONT = X_DISK_BACK + T_DISK  # 2414 -> 111 mm before the throat


def wall_radius(x):
    """Inner wall radius at axial position x (mm)."""
    pts = [(0.0, D_IN / 2), (X_CONV, D_IN / 2), (X_THROAT, D_THROAT / 2),
           (X_DIFF, D_THROAT / 2), (X_END, D_OUT / 2)]
    for (x0, r0), (x1, r1) in zip(pts, pts[1:]):
        if x0 <= x <= x1:
            return r0 + (r1 - r0) * (x - x0) / (x1 - x0)
    raise ValueError(x)
