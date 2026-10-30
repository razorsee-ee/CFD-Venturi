"""Pictures and sheet parameters from the axisymmetric VOF run (no ParaView needed).

Needs VTK files from:  foamToVTK -ascii -fields '(alpha.oil U)' [-latestTime]
Writes:
  sheet_<step>.png   oil volume fraction and velocity magnitude, mirrored about
                     the axis so the jet, the disk and the sheet look like a section
  sheet_summary.txt  for the latest time: sheet angle, speed and thickness at
                     several radii (what case-rim needs for its injection)
  sheet.gif          animation of all frames (with --gif)
Usage (in nozzle-vof/):  python plot_sheet.py [--all | --gif]
  (no option: only the latest time; frames already drawn are not redrawn)
"""
import math
import re
import sys
from pathlib import Path

import numpy as np

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

CASE = Path(__file__).resolve().parent
DISK_X, DISK_R, NOZZLE_R = (0.092, 0.100), 0.0325, 0.014


def read_vtk(path):
    """Cell centres (x, r) and cell data of a legacy ascii unstructured-grid VTK."""
    tok = path.read_text(encoding="utf-8", errors="ignore").split()
    i = tok.index("POINTS")
    n = int(tok[i + 1])
    P = np.array(tok[i + 3:i + 3 + 3 * n], float).reshape(n, 3)
    i = tok.index("CELLS")
    nc, ntok = int(tok[i + 1]), int(tok[i + 2])
    conn = np.array(tok[i + 3:i + 3 + ntok], int)
    C = np.zeros((nc, 3))
    k = 0
    for c in range(nc):
        m = conn[k]
        C[c] = P[conn[k + 1:k + 1 + m]].mean(0)
        k += 1 + m
    data = {}
    for name, ncomp in (("alpha.oil", 1), ("U", 3)):
        for j in range(len(tok) - 3):
            if tok[j] == name and tok[j + 1] == str(ncomp) and tok[j + 2] == str(nc):
                data[name] = np.array(tok[j + 4:j + 4 + ncomp * nc], float).reshape(nc, ncomp).squeeze()
                break
    x, r = C[:, 0], np.hypot(C[:, 1], C[:, 2])
    return x, r, data


def sheet_parameters(x, r, a, U):
    """Sheet position, angle, speed, thickness on radial stations beyond the disk rim."""
    rows = []
    stations = [0.035, 0.045, 0.06, 0.08, 0.10, 0.12, 0.14]
    xs = []
    for rs in stations:
        sel = (np.abs(r - rs) < 0.0015) & (a > 0.5)
        if sel.sum() < 2:
            rows.append((rs, None, None, None))
            xs.append(None)
            continue
        xm = np.average(x[sel], weights=a[sel])
        # thickness: extent along x of the liquid at this radius (sheet ~ normal to r)
        thick = (x[sel].max() - x[sel].min()) + 0.4e-3
        speed = np.average(np.linalg.norm(U[sel], axis=1), weights=a[sel])
        rows.append((rs, xm, thick, speed))
        xs.append(xm)
    out = []
    for k, (rs, xm, th, sp) in enumerate(rows):
        if xm is None:
            out.append(f"  r = {rs * 1e3:5.0f} мм: масла нет (плёнка не дошла или распалась)")
            continue
        ang = ""
        if k > 0 and xs[k - 1] is not None:
            dx, dr = xm - xs[k - 1], rs - stations[k - 1]
            ang = f", угол к оси {math.degrees(math.atan2(dr, dx)):5.1f}°"
        out.append(f"  r = {rs * 1e3:5.0f} мм: x = {xm * 1e3:6.1f} мм, толщина ~{th * 1e3:4.1f} мм, "
                   f"скорость {sp:5.1f} м/с{ang}")
    return out


def plot(x, r, data, t, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    X = np.r_[x, x]
    Y = np.r_[r, -r]
    tri = mtri.Triangulation(X * 1e3, Y * 1e3)
    # mask triangles inside the solids / crossing them (long edges)
    xm = X[tri.triangles].mean(1)
    ym = np.abs(Y[tri.triangles]).mean(1)
    solid = ((xm >= DISK_X[0]) & (xm <= DISK_X[1]) & (ym < DISK_R)) | ((xm < 0) & (ym < NOZZLE_R))
    edge = np.max([np.hypot(*(np.c_[X, Y][tri.triangles[:, a]] - np.c_[X, Y][tri.triangles[:, b]]).T)
                   for a, b in ((0, 1), (1, 2), (2, 0))], axis=0)
    tri.set_mask(solid | (edge > 0.006))

    fig, ax = plt.subplots(2, 1, figsize=(8.5, 10), sharex=True, layout="constrained")
    a = np.r_[data["alpha.oil"], data["alpha.oil"]]
    pc = ax[0].tripcolor(tri, a, cmap="Blues", vmin=0, vmax=1, shading="gouraud")
    fig.colorbar(pc, ax=ax[0], label="доля масла (1 = масло)")
    Umag = np.linalg.norm(data["U"], axis=1)
    pc = ax[1].tripcolor(tri, np.r_[Umag, Umag], cmap="viridis", shading="gouraud", vmin=0, vmax=100)  # fixed: frames comparable
    fig.colorbar(pc, ax=ax[1], label="|U|, м/с (выше 100 — тем же цветом)", extend="max")
    ax[1].tricontour(tri, a, levels=[0.5], colors="w", linewidths=0.7)
    for k in range(2):
        a_ = ax[k]
        a_.add_patch(plt.Rectangle((DISK_X[0] * 1e3, -DISK_R * 1e3), (DISK_X[1] - DISK_X[0]) * 1e3,
                                   2 * DISK_R * 1e3, fc="0.4"))
        a_.add_patch(plt.Rectangle((-150, -NOZZLE_R * 1e3), 150, 2 * NOZZLE_R * 1e3, fc="0.4"))
        a_.set_aspect("equal")
        a_.set_xlim(-40, 300)
        a_.set_ylim(-150, 150)
        a_.set_ylabel("r, мм (зеркально)")
    ax[1].set_xlabel("x, мм (сопло — 0, диск — 92…100; газ слева направо)")
    ax[0].set_title(f"Струя масла → диск → плёнка, t = {t * 1e3:.1f} мс")
    fig.savefig(path, dpi=110)
    plt.close(fig)


def main():
    files = sorted((CASE / "VTK").glob("*.vtk"), key=lambda f: int(re.findall(r"_(\d+)\.vtk$", f.name)[0]))
    files = [f for f in files if "_" in f.stem and f.parent.name == "VTK"]
    if not files:
        sys.exit("no VTK files: run  foamToVTK -ascii -fields '(alpha.oil U)'  first")
    # foamToVTK numbers files by time-step index; <time>/uniform/time holds that index
    times = {}
    for d in CASE.iterdir():
        u = d / "uniform" / "time"
        if re.fullmatch(r"[0-9.eE+-]+", d.name) and u.is_file():
            m = re.search(r"^index\s+(\d+);", u.read_text(errors="ignore"), re.M)
            if m:
                times[int(m.group(1))] = float(d.name)
    sel = files if ("--all" in sys.argv or "--gif" in sys.argv) else files[-1:]
    frames = []
    for f in sel:
        idx = int(re.findall(r"_(\d+)\.vtk$", f.name)[0])
        t = times.get(idx, float("nan"))
        png = CASE / f"sheet_{idx:05d}.png"
        frames.append(png)
        if png.exists() and f is not files[-1]:
            continue                    # already drawn by an earlier call
        x, r, data = read_vtk(f)
        plot(x, r, data, t, png)
        print(f"written {png.name} (t = {t * 1e3:.2f} ms)")
    if "--gif" in sys.argv and len(frames) > 1:
        from PIL import Image
        im = [Image.open(f).convert("RGB") for f in frames]
        w = 800
        im = [i.resize((w, round(i.height * w / i.width)), Image.LANCZOS).quantize(64) for i in im]
        # 0.15 s per frame, the last one held for 2 s, endless loop
        im[0].save(CASE / "sheet.gif", save_all=True, append_images=im[1:],
                   duration=[150] * (len(im) - 1) + [2000], loop=0, optimize=True)
        print(f"written sheet.gif ({len(im)} frames)")
    x, r, data = read_vtk(files[-1])
    lines = [f"Плёнка масла за кромкой диска (последний момент, {files[-1].name}):"]
    lines += sheet_parameters(x, r, data["alpha.oil"], data["U"])
    txt = "\n".join(lines)
    print(txt)
    (CASE / "sheet_summary.txt").write_text(txt + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
