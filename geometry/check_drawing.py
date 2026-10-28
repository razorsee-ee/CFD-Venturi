"""Render a check drawing of the model: axial section with dimensions, nozzle detail, 3D view of the nozzle."""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.patches import Rectangle, Polygon

import params as p

OUT = Path(__file__).parent / "out"
STEEL, WALL = "#868e96", "#212529"


def read_stl(path):
    v = [list(map(float, l.split()[1:])) for l in open(path) if l.strip().startswith("vertex")]
    return np.array(v).reshape(-1, 3, 3)


def nozzle_section(a):
    """Nozzle in the vertical symmetry plane (z = 0)."""
    top = p.wall_radius(p.X_STEM)
    a.add_patch(Rectangle((p.X_STEM - p.D_STEM / 2, -p.STEM_BELOW_AXIS), p.D_STEM, top + p.STEM_BELOW_AXIS, fc=STEEL, ec=WALL))
    a.add_patch(Rectangle((p.X_STEM, -p.D_BODY / 2), p.L_BODY, p.D_BODY, fc=STEEL, ec=WALL))
    a.add_patch(Rectangle((p.X_DISK_BACK, -p.D_DISK / 2), p.T_DISK, p.D_DISK, fc=STEEL, ec=WALL))
    for k in (1, -0.5):
        a.plot([p.X_NOZZLE, p.X_DISK_BACK], [k * p.R_STRUT_BODY, k * p.R_STRUT_DISK], color=WALL, lw=2)


def hdim(a, x0, x1, y, txt, above=False):
    a.annotate("", (x0, y), (x1, y), arrowprops=dict(arrowstyle="<->", lw=0.8))
    a.text((x0 + x1) / 2, y + (8 if above else -8), txt, ha="center", va="bottom" if above else "top", fontsize=9)


def vdim(a, x, r, txt):
    a.annotate("", (x, -r), (x, r), arrowprops=dict(arrowstyle="<->", lw=0.8))
    a.text(x + 15, 0, txt, rotation=90, va="center", fontsize=9, bbox=dict(fc="white", ec="none", pad=1))


fig = plt.figure(figsize=(17, 10.5))
ax = fig.add_axes([0.02, 0.45, 0.96, 0.5])
xs = np.array([0, p.X_CONV, p.X_THROAT, p.X_DIFF, p.X_END])
rs = np.array([p.wall_radius(x) for x in xs])
ax.fill_between(xs, -rs, rs, color="#e7f5ff")
for s in (1, -1):
    ax.plot(xs, s * rs, color=WALL, lw=2)
ax.plot([0, 0], [-rs[0], rs[0]], color="#2f9e44", lw=3)
ax.plot([p.X_END, p.X_END], [-rs[-1], rs[-1]], color="#e8590c", lw=3)
ax.text(-30, 0, "inlet", rotation=90, va="center", ha="right", color="#2f9e44")
ax.text(p.X_END + 30, 0, "outlet", rotation=90, va="center", color="#e8590c")
ax.axhline(0, color="#adb5bd", ls="-.", lw=0.7)
for x in xs[1:-1]:
    ax.axvline(x, color="#ced4da", lw=0.6, ls=":")
nozzle_section(ax)
y = -rs.max() - 60
hdim(ax, 0, p.X_CONV, y, f"цилиндр {p.L_CYL:.0f}")
hdim(ax, p.X_CONV, p.X_THROAT, y, f"конфузор {p.L_CONV:.0f}")
hdim(ax, p.X_THROAT, p.X_DIFF, y - 70, f"горловина {p.L_THROAT:.0f}")
hdim(ax, p.X_DIFF, p.X_END, y, f"диффузор {p.L_DIFF:.0f}")
hdim(ax, 0, p.X_END, y - 150, f"всего {p.X_END:.0f}")
hdim(ax, p.X_STEM, p.X_THROAT, rs.max() + 40, f"{p.X_THROAT - p.X_STEM:.0f}", above=True)
vdim(ax, 300, p.D_IN / 2, f"⌀{p.D_IN:.0f}")
vdim(ax, p.X_THROAT + 124, p.D_THROAT / 2, f"⌀{p.D_THROAT:.0f}")
vdim(ax, p.X_END - 200, p.wall_radius(p.X_END - 200), f"⌀{p.D_OUT:.0f} на выходе (≈, снято с чертежа)")
ang_c = np.degrees(np.arctan((p.D_IN - p.D_THROAT) / 2 / p.L_CONV))
ang_d = np.degrees(np.arctan((p.D_OUT - p.D_THROAT) / 2 / p.L_DIFF))
ax.text(p.X_CONV + 400, -220, f"полуугол {ang_c:.1f}°", fontsize=9)
ax.text(p.X_DIFF + 900, -260, f"полуугол {ang_d:.1f}°", fontsize=9)
ax.set_aspect("equal"); ax.axis("off")
ax.set_xlim(-120, p.X_END + 120); ax.set_ylim(y - 230, rs.max() + 110)
ax.set_title("Газовая полость скруббера (модель), продольный разрез, мм", fontsize=13)

# nozzle detail
ad = fig.add_axes([0.03, 0.03, 0.5, 0.4])
x0, x1 = p.X_STEM - 120, p.X_DIFF + 40
xx = np.linspace(x0, x1, 200)
rr = np.array([p.wall_radius(x) for x in xx])
ad.plot(xx, rr, color=WALL, lw=2); ad.plot(xx, -rr, color=WALL, lw=2)
ad.axhline(0, color="#adb5bd", ls="-.", lw=0.7)
nozzle_section(ad)
hdim(ad, p.X_STEM, p.X_NOZZLE, -70, f"{p.L_BODY:.0f}")
hdim(ad, p.X_NOZZLE, p.X_DISK_FRONT, -70, f"{p.X_DISK_FRONT - p.X_NOZZLE:.0f}")
hdim(ad, p.X_DISK_FRONT, p.X_THROAT, -70, f"{p.X_THROAT - p.X_DISK_FRONT:.0f}")
hdim(ad, p.X_THROAT, p.X_DIFF, -70, f"{p.L_THROAT:.0f}")
ad.text(p.X_STEM + 30, 120, f"лансет ⌀{p.D_STEM}", fontsize=9)
ad.text(p.X_STEM + 5, 22, f"корпус ⌀{p.D_BODY:.0f}", fontsize=9)
ad.text(p.X_DISK_FRONT + 8, 30, f"диск ⌀{p.D_DISK:.0f}×{p.T_DISK:.0f}", fontsize=9)
ad.text(p.X_NOZZLE + 5, -45, "3 стойки ⌀6 (допущение)", fontsize=9)
ad.set_aspect("equal"); ad.axis("off")
ad.set_xlim(x0, x1); ad.set_ylim(-p.D_THROAT / 2 - 50, p.wall_radius(x0) + 20)
ad.set_title("Узел форсунки (в разрезе по оси лансета)", fontsize=12)

# 3D nozzle from the exported STL
a3 = fig.add_axes([0.56, 0.02, 0.42, 0.42], projection="3d")
tri = read_stl(OUT / "nozzle.stl") * 1e3
a3.add_collection3d(Poly3DCollection(tri[:, :, [0, 2, 1]], facecolors="#adb5bd", shade=True,
                                     ))
c = np.array([p.X_NOZZLE, 0, 0])
a3.set_xlim(c[0] - 110, c[0] + 110); a3.set_ylim(-110, 110); a3.set_zlim(-80, 140)
a3.set_box_aspect((1, 1, 1)); a3.view_init(elev=20, azim=-60)
a3.set_xlabel("x (по потоку)"); a3.set_ylabel("z"); a3.set_zlabel("y")
a3.set_title("Форсунка из STL (nozzle.stl)", fontsize=12)

fig.savefig(OUT / "check_drawing.png", dpi=120)
