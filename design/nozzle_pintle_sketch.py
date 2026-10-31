"""Simple sketch: project nozzle with the disk on a central stem held by a cross in the tube."""
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, Polygon, FancyArrowPatch

D_CROSS, D_EXIT, D_STEM, D_DISK = 26.0, 22.4, 10.0, 65.0   # mm
GAP = 92.0                                                 # nozzle exit -> disk face
A0 = math.pi / 4 * 20 ** 2                                 # original D20 orifice
A_exit = math.pi / 4 * (D_EXIT ** 2 - D_STEM ** 2)
A_cross = math.pi / 4 * (D_CROSS ** 2 - D_STEM ** 2) - 4 * 4 * (D_CROSS / 2 - D_STEM / 2)

fig = plt.figure(figsize=(13, 5.6))
ax = fig.add_axes([0.03, 0.08, 0.66, 0.84])
bx = fig.add_axes([0.72, 0.18, 0.26, 0.62])
steel, oil = "0.55", "#3b7dd8"

# side section (x right = gas flow, oil flows the same way)
W = 6.0                                     # tube wall
x0, x_cr, x_tp, x_ex = -120, -60, -25, 0    # tube start, cross, taper start, exit
for s in (1, -1):
    ax.add_patch(Polygon([(x0, s * D_CROSS / 2), (x_tp, s * D_CROSS / 2), (x_ex, s * D_EXIT / 2),
                          (x_ex, s * (D_EXIT / 2 + W)), (x_tp, s * (D_CROSS / 2 + W)), (x0, s * (D_CROSS / 2 + W))],
                         fc=steel, ec="k"))
    ax.add_patch(Rectangle((x_cr - 2, s * D_STEM / 2), 4, s * (D_CROSS - D_STEM) / 2, fc="0.3", ec="k"))
ax.add_patch(Rectangle((x_cr - 2, -D_STEM / 2), GAP - x_cr + 2, D_STEM, fc="0.7", ec="k"))     # stem
ax.add_patch(Rectangle((GAP, -D_DISK / 2), 8, D_DISK, fc=steel, ec="k"))                     # disk
# oil: annular jet and sheet
for s in (1, -1):
    ax.add_patch(Polygon([(x_ex, s * D_STEM / 2), (x_ex, s * D_EXIT / 2), (GAP - 6, s * D_EXIT / 2 * 1.05),
                          (GAP, s * D_EXIT / 2 * 1.6), (GAP, s * D_STEM / 2)], fc=oil, alpha=0.35, ec="none"))
    ax.add_patch(Polygon([(GAP, s * D_EXIT / 2 * 1.6), (GAP - 1, s * D_DISK / 2), (GAP + 6, s * 70),
                          (GAP + 9, s * 70), (GAP + 2, s * D_DISK / 2), (GAP + 0.5, s * D_EXIT / 2 * 1.6)],
                         fc=oil, alpha=0.35, ec="none"))
ax.add_patch(FancyArrowPatch((-115, 0), (-75, 0), arrowstyle="-|>", mutation_scale=18, color=oil, lw=2))
ax.text(-117, 4, "масло", color=oil, fontsize=11)
ax.add_patch(FancyArrowPatch((-110, 60), (-50, 60), arrowstyle="-|>", mutation_scale=18, color="0.4", lw=2))
ax.text(-110, 64, "газ", color="0.3", fontsize=11)
for s in (1, -1):
    ax.add_patch(FancyArrowPatch((GAP + 4, s * 45), (GAP + 9, s * 72), arrowstyle="-|>", mutation_scale=14,
                                 color=oil, lw=1.5))

def dim(x1, x2, y, t):
    ax.annotate("", (x1, y), (x2, y), arrowprops=dict(arrowstyle="<->", lw=0.8))
    ax.text((x1 + x2) / 2, y + 2, t, ha="center", fontsize=9)

dim(x_ex, GAP, -48, f"{GAP:.0f} мм (как в проекте)")
ax.annotate(f"крест: 2 полосы 4 мм\n(лучше 3–5 диаметров\nот выхода, а не у самого носика)",
            (x_cr, D_CROSS / 2), (-118, -75), fontsize=9, arrowprops=dict(arrowstyle="->", lw=0.8))
ax.annotate(f"центральная ножка ⌀{D_STEM:.0f}", (40, D_STEM / 2), (25, 32), fontsize=9,
            arrowprops=dict(arrowstyle="->", lw=0.8))
ax.annotate(f"канал ⌀{D_CROSS:.0f}", (-90, D_CROSS / 2 + W), (-100, 40), fontsize=9,
            arrowprops=dict(arrowstyle="->", lw=0.8))
ax.annotate(f"сужение до ⌀{D_EXIT:.1f} на выходе", (x_ex, D_EXIT / 2 + W), (-30, 48), fontsize=9,
            arrowprops=dict(arrowstyle="->", lw=0.8))
ax.annotate(f"диск ⌀{D_DISK:.0f}, без стоек", (GAP + 8, -D_DISK / 2), (GAP + 15, -62), fontsize=9,
            arrowprops=dict(arrowstyle="->", lw=0.8))
ax.text(x_cr - 4, -D_CROSS / 2 - W - 8, "A", fontsize=12, weight="bold")
ax.plot([x_cr, x_cr], [-D_CROSS / 2 - W - 4, D_CROSS / 2 + W + 4], "k-.", lw=0.8)
ax.set_xlim(-125, 135); ax.set_ylim(-85, 85); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("Форсунка с диском на центральной ножке (разрез вдоль оси, мм)")

# section A-A
bx.add_patch(Circle((0, 0), D_CROSS / 2 + W, fc=steel, ec="k"))
bx.add_patch(Circle((0, 0), D_CROSS / 2, fc=oil, alpha=0.35, ec="k"))
for ang in (0, 90):
    bx.add_patch(Rectangle((-D_CROSS / 2, -2), D_CROSS, 4, angle=ang, rotation_point=(0, 0), fc="0.3", ec="k"))
bx.add_patch(Circle((0, 0), D_STEM / 2, fc="0.7", ec="k"))
bx.set_xlim(-22, 22); bx.set_ylim(-22, 22); bx.set_aspect("equal"); bx.axis("off")
bx.set_title("сечение A: крест держит ножку")
bx.text(0, -27, f"проход по маслу, мм²:\nсопло ⌀20 сейчас: {A0:.0f}\nвыход ⌀{D_EXIT:.1f} − ножка ⌀{D_STEM:.0f}: {A_exit:.0f}\n"
        f"в кресте (⌀{D_CROSS:.0f} − ножка − полосы): {A_cross:.0f}",
        ha="center", va="top", fontsize=9)
fig.savefig("nozzle_pintle_sketch.png", dpi=110)
print(round(A0), round(A_exit), round(A_cross))
