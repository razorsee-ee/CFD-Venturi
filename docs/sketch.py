"""Concept sketch of the test Venturi scrubber (hot air + water, deflector-disk nozzle)."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, Polygon, FancyArrowPatch, Arc

# ---- parameters, mm ----
R_in, R_th = 100.0, 50.0          # inlet / throat radius (constriction 2x by diameter)
L_up, L_conv, L_th, L_div, L_out = 500, 225, 100, 820, 150
d_tube, d_tube_i = 12.0, 8.0      # nozzle tube OD / ID
D_disk, t_disk = 3 * d_tube, 3.0  # deflector disk: 3x tube OD
x_tube_end, x_disk = -45.0, -25.0 # tube outlet / disk front face
x_tube_in, R_bend = -260.0, 25.0  # where tube passes the wall, bend radius

x0 = -L_up; x1 = 0; x2 = L_conv; x3 = x2 + L_th; x4 = x3 + L_div; x5 = x4 + L_out
wall_x = [x0, x1, x2, x3, x4, x5]
wall_r = [R_in, R_in, R_th, R_th, R_in, R_in]

GAS, WATER, STEEL, WALL = "#d9480f", "#1c7ed6", "#495057", "#212529"

fig = plt.figure(figsize=(16, 10))
ax = fig.add_axes([0.03, 0.42, 0.94, 0.53])
axd = fig.add_axes([0.03, 0.03, 0.45, 0.36])
axe = fig.add_axes([0.55, 0.03, 0.40, 0.36])

def draw_tube_and_disk(a, lw=1.5):
    ro, ri = d_tube / 2, d_tube_i / 2
    xb = x_tube_in + R_bend                              # end of bend, horizontal part starts
    # vertical leg (top wall -> bend)
    a.add_patch(Rectangle((x_tube_in - ro, R_bend), 2 * ro, R_in + 15 - R_bend, fc="#ced4da", ec=STEEL, lw=lw))
    # bend (quarter torus in section)
    th = np.linspace(np.pi, 1.5 * np.pi, 40)
    outer = np.c_[xb + (R_bend + ro) * np.cos(th), R_bend + (R_bend + ro) * np.sin(th)]
    inner = np.c_[xb + (R_bend - ro) * np.cos(th[::-1]), R_bend + (R_bend - ro) * np.sin(th[::-1])]
    a.add_patch(Polygon(np.r_[outer, inner], fc="#ced4da", ec=STEEL, lw=lw))
    # horizontal leg along axis
    a.add_patch(Rectangle((xb, -ro), x_tube_end - xb, 2 * ro, fc="#ced4da", ec=STEEL, lw=lw))
    a.add_patch(Rectangle((xb, -ri), x_tube_end - xb, 2 * ri, fc="#a5d8ff", ec="none"))
    # disk
    a.add_patch(Rectangle((x_disk, -D_disk / 2), t_disk, D_disk, fc="#868e96", ec=STEEL, lw=lw))
    # 3 rods at 120 deg: in section we see the one at 90 deg (top) and the projection of two others
    for k in (1.0, -0.5):  # struts at 90 deg (in plane) and 210/330 deg (projected)
        a.plot([x_tube_end - 12, x_disk], [k * ro, k * (D_disk / 2 - 3)], color=STEEL, lw=lw * 1.3)

# ---------- main section ----------
xs = np.array(wall_x); rs = np.array(wall_r)
for s in (1, -1):
    ax.plot(xs, s * rs, color=WALL, lw=2.5)
    ax.fill_between(xs, s * rs, s * (rs + 12), color="#adb5bd")
ax.axhline(0, color="#868e96", ls=(0, (12, 4, 2, 4)), lw=0.8)
draw_tube_and_disk(ax)
# pass-through in wall
ax.add_patch(Rectangle((x_tube_in - 6, R_in), 12, 12, fc="#ced4da", ec=STEEL))
ax.annotate("вода, 20–30 °C", (x_tube_in, R_in + 14), (x_tube_in - 130, R_in + 60),
            color=WATER, fontsize=11, arrowprops=dict(arrowstyle="->", color=WATER, lw=1.5))
# gas arrows
for y in (-70, -35, 35, 70):
    ax.add_patch(FancyArrowPatch((x0 + 10, y), (x0 + 120, y), arrowstyle="-|>", mutation_scale=18, color=GAS, lw=2))
ax.text(x0 + 15, 0, "горячий воздух ~470 °C", color=GAS, fontsize=11, va="center")
ax.add_patch(FancyArrowPatch((x4 + 20, 0), (x5 + 60, 0), arrowstyle="-|>", mutation_scale=18, color=GAS, lw=2))
ax.text(x4 + 30, 12, "газ + капли", color=GAS, fontsize=11)
# water umbrella sheet from disk edge + droplets carried to throat
rng = np.random.default_rng(3)
for s in (1, -1):
    t = np.linspace(0, 1, 30)
    xsheet = x_disk + 2 + 140 * t
    rsheet = s * (D_disk / 2 + 55 * np.sqrt(t))
    ax.plot(xsheet, rsheet, color=WATER, lw=1.8)
    n = 140
    xd = rng.uniform(x_disk + 60, x3 + 250, n)
    # droplets confined by local wall radius
    rwall = np.interp(xd, xs, rs)
    rd = s * rng.uniform(0.15, 0.85, n) * rwall
    ax.scatter(xd, rd, s=4, color=WATER, alpha=0.6)
# dimensions
def dim_v(a, x, r, txt):
    a.annotate("", (x, -r), (x, r), arrowprops=dict(arrowstyle="<->", lw=1))
    a.text(x + 6, -r * 0.5, txt, fontsize=10, rotation=90, va="center", bbox=dict(fc="white", ec="none", pad=1))
dim_v(ax, x0 + 250, R_in, "D = 200")
dim_v(ax, x2 + 50, R_th, "d = 100")
def dim_h(a, xa, xb, y, txt):
    a.annotate("", (xa, y), (xb, y), arrowprops=dict(arrowstyle="<->", lw=1))
    a.text((xa + xb) / 2, y - 5, txt, fontsize=10, ha="center", va="top")
dim_h(ax, x1, x2, -R_in - 35, "конфузор 225\nполуугол 12.5°")
dim_h(ax, x2, x3, -R_in - 35, "горло-\nвина 100")
dim_h(ax, x3, x4, -R_in - 35, "диффузор 820 (полуугол ~3.5°)")
dim_h(ax, x0, x1, -R_in - 35, "подводящая труба 500")
ax.annotate("диск-отражатель ⌀36 = 3×⌀12", (x_disk, D_disk / 2), (x_disk - 120, R_in + 20), fontsize=10, arrowprops=dict(arrowstyle="->"))
ax.set_xlim(x0 - 20, x5 + 90); ax.set_ylim(-R_in - 80, R_in + 70)
ax.set_aspect("equal"); ax.axis("off")
ax.set_title("Тестовый скруббер Вентури — продольный разрез (гипотетические размеры, мм)", fontsize=14)

# ---------- detail of nozzle ----------
draw_tube_and_disk(axd, lw=2)
axd.axhline(0, color="#868e96", ls=(0, (12, 4, 2, 4)), lw=0.8)
axd.plot([x0, x1, x2], [R_in, R_in, R_th], color=WALL, lw=2.5)
axd.add_patch(FancyArrowPatch((x_tube_end + 1, 0), (x_disk - 1, 0), arrowstyle="-|>", mutation_scale=14, color=WATER, lw=2.5))
for s in (1, -1):
    t = np.linspace(0, 1, 30)
    axd.plot(x_disk + 2 + 60 * t, s * (D_disk / 2 + 35 * np.sqrt(t)), color=WATER, lw=2)
    axd.plot(x_disk + 1 + 0 * t, s * np.linspace(1, D_disk / 2, 30), color=WATER, lw=2)
axd.annotate("загнутая трубка ⌀12/⌀8\nструя бьёт в диск", (x_tube_end - 40, d_tube / 2), (-250, -60),
             fontsize=10, arrowprops=dict(arrowstyle="->"))
axd.annotate("3 стойки под 120°", (x_tube_end - 8, d_tube / 2 + 1.5), (-150, 60),
             fontsize=10, arrowprops=dict(arrowstyle="->"))
axd.annotate("зазор 20", (x_tube_end + 10, -2), (-85, -40), fontsize=10, arrowprops=dict(arrowstyle="->"))
axd.annotate("водяная плёнка-«зонтик»\n→ распад на капли", (x_disk + 50, -D_disk / 2 - 30), (20, -80),
             fontsize=10, color=WATER, arrowprops=dict(arrowstyle="->", color=WATER))
axd.annotate("начало конфузора", (0, R_in), (-60, R_in + 15), fontsize=10, arrowprops=dict(arrowstyle="->"))
axd.set_xlim(-300, 120); axd.set_ylim(-R_in - 5, R_in + 30)
axd.set_aspect("equal"); axd.axis("off")
axd.set_title("Узел форсунки (увеличено)", fontsize=12)

# ---------- end view (looking upstream from the throat) ----------
axe.add_patch(Circle((0, 0), R_in, fc="none", ec=WALL, lw=2.5))
axe.add_patch(Circle((0, 0), R_th, fc="none", ec=WALL, lw=1, ls="--"))
axe.add_patch(Circle((0, 0), D_disk / 2, fc="#868e96", ec=STEEL, lw=1.5))
axe.add_patch(Circle((0, 0), d_tube / 2, fc="#a5d8ff", ec=STEEL, lw=1.5))
for ang in (90, 210, 330):
    a = np.radians(ang); r = 11
    axe.add_patch(Circle((r * np.cos(a), r * np.sin(a)), 1.5, fc=WALL))
for ang in np.linspace(0, 2 * np.pi, 16, endpoint=False):
    axe.add_patch(FancyArrowPatch((22 * np.cos(ang), 22 * np.sin(ang)), (60 * np.cos(ang), 60 * np.sin(ang)),
                                  arrowstyle="-|>", mutation_scale=9, color=WATER, lw=1))
axe.text(0, -R_in - 18, "струя из трубки, 3 стойки, за ними диск ⌀36, радиальный выброс воды\n"
         "(разрез между торцом трубки и диском); пунктир — горловина ⌀100",
         ha="center", va="top", fontsize=10)
axe.set_xlim(-R_in - 10, R_in + 10); axe.set_ylim(-R_in - 45, R_in + 10)
axe.set_aspect("equal"); axe.axis("off")
axe.set_title("Сечение A–A, вид по потоку", fontsize=12)

fig.savefig(__file__.replace("sketch.py", "concept_sketch.png"), dpi=130)
