"""Frames, animation and sheet break-up numbers for the 3D sector VOF run (no ParaView needed).

Reads postProcessing/surfaces/<time>/{slice,iso}.vtk written every 0.1 ms by system/surfaces:
  slice - plane through the axis along the diagonal of the quarter (x-r section)
  iso   - oil surface alpha.oil = 0.5 (sheet, ligaments, drops)
Writes:
  frames/f_<time>.png  three views: section; "shadow photo" from the side and from
                       downstream (the 90 deg quarter copied around the axis to 360 deg)
  sheet3d.gif          all frames (with --gif)
  sheet3d_summary.txt  per time: radius up to which the sheet is still one piece,
                       number of pieces torn off, oil sheet angle/speed at the rim
Usage (in nozzle-vof-3d/):  python plot3d.py [--gif]   (frames already drawn are kept)
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
SURF = CASE / "postProcessing" / "surfaces"
OUT = CASE / "frames"
DISK_X, DISK_R = (0.092, 0.100), 0.0325
LANCE_X, LANCE_R = 0.082, 0.014


SHADOW_MAX = 4.0                        # fixed colour scale of the shadow pictures (surfaces crossed)
SECTOR = 90.0                           # the mesh is a quarter with rotational cyclic sides


def read_vtk(path):
    """Legacy VTK POLYDATA (binary or ascii) from the OpenFOAM surface writer."""
    raw = path.read_bytes()
    pos = 0
    binary = False

    def line():
        nonlocal pos
        end = raw.index(b"\n", pos)
        s = raw[pos:end].decode("latin-1").strip()
        pos = end + 1
        return s

    def block(n, dt):
        nonlocal pos
        if binary:
            a = np.frombuffer(raw, dtype=">" + dt, count=n, offset=pos)
            pos += n * a.itemsize
            if pos < len(raw) and raw[pos:pos + 1] == b"\n":
                pos += 1
            return a.astype(float if dt[0] == "f" else np.int64)
        vals = []
        while len(vals) < n:
            vals += line().split()
        return np.array(vals, float if dt[0] == "f" else np.int64)

    pts, polys, data = None, None, {}
    while pos < len(raw):
        s = line()
        if not s:
            continue
        w = s.split()
        if w[0] == "BINARY":
            binary = True
        elif w[0] == "POINTS":
            pts = block(3 * int(w[1]), "f4" if w[2] == "float" else "f8").reshape(-1, 3)
        elif w[0] == "POLYGONS":
            polys = block(int(w[2]), "i4")
        elif w[0] == "FIELD":
            for _ in range(int(w[2])):
                h = line().split()
                while not h:
                    h = line().split()
                nc, nt = int(h[1]), int(h[2])
                data[h[0]] = block(nc * nt, "f4" if h[3] == "float" else "f8").reshape(nt, nc).squeeze()
    # polygons -> list of point index arrays
    faces, k = [], 0
    while polys is not None and k < len(polys):
        n = polys[k]
        faces.append(polys[k + 1:k + 1 + n])
        k += 1 + n
    return pts, faces, data


def triangles(faces):
    tri = []
    for f in faces:
        for i in range(1, len(f) - 1):
            tri.append((f[0], f[i], f[i + 1]))
    return np.array(tri, np.int64).reshape(-1, 3)


def components(tri, npts):
    """Connected pieces of the oil surface (label per point), by label propagation."""
    lab = np.arange(npts)
    while True:
        m = lab[tri].min(1)
        new = lab.copy()
        for c in range(3):
            np.minimum.at(new, tri[:, c], m)
        new = new[new]                      # pointer jumping
        if np.array_equal(new, lab):
            return lab
        lab = new


def breakup(pts, tri):
    """Intact sheet radius and torn-off pieces from the iso surface."""
    if len(tri) == 0:
        return None
    lab = components(tri, len(pts))
    r = np.hypot(pts[:, 1], pts[:, 2])
    x = pts[:, 0]
    # the piece fed by the jet: has points in the jet gap in front of the disk
    jet = (x > LANCE_X) & (x < DISK_X[0]) & (r < 0.012)
    if not jet.any():
        return None
    main = np.bincount(lab[jet]).argmax()
    used = np.unique(tri)
    pieces = np.unique(lab[used])
    others = pieces[pieces != main]
    r_main = r[lab == main].max()
    far = [r[lab == p].mean() for p in others]
    return r_main, len(others), np.array(far)


def sheet_rim(pts, faces, data):
    """Oil on the slice just outside the disk rim: angle, speed, thickness."""
    if "alpha.oil" not in data or not faces:
        return ""
    c = np.array([pts[f].mean(0) for f in faces])
    c[:, 1] = np.hypot(c[:, 1], c[:, 2])       # radius
    a = data["alpha.oil"]
    U = data.get("U")
    rows = []
    for r0, r1 in ((0.035, 0.045), (0.045, 0.060)):
        sel = (c[:, 1] > r0) & (c[:, 1] < r1) & (a > 0.5)
        if sel.sum() < 3:
            return ""
        rows.append((np.average(c[sel, 0], weights=a[sel]), np.average(c[sel, 1], weights=a[sel]),
                     np.average(np.linalg.norm(U[sel], axis=1), weights=a[sel]) if U is not None else float("nan")))
    (x0, y0, u0), (x1, y1, u1) = rows
    ang = math.degrees(math.atan2(y1 - y0, x1 - x0))
    return f"угол к оси {ang:4.0f}°, скорость {0.5 * (u0 + u1):4.1f} м/с (r 35–60 мм)"


def shadow(ax, u, v, w, extent, bins):
    """Number of oil surfaces a line of sight crosses (projected area / pixel area).
    Fixed scale for all frames: 0 (gas) ... 4 and more (several sheets/drops in a row)."""
    H, _, _ = np.histogram2d(u, v, bins=bins, range=[extent[:2], extent[2:]], weights=w)
    pix = (extent[1] - extent[0]) / bins[0] * (extent[3] - extent[2]) / bins[1]
    return ax.imshow(H.T / pix, origin="lower", extent=extent, cmap="Greys", vmin=0, vmax=SHADOW_MAX,
                     interpolation="bilinear")


def frame(t, path, sector):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    fig, ax = plt.subplots(1, 3, figsize=(16, 5.6), layout="constrained",
                           gridspec_kw={"width_ratios": [1.15, 1.15, 1]})
    # 1) section
    sl = path / "slice.vtk"
    if sl.exists():
        pts, faces, data = read_vtk(sl)
        a = data.get("alpha.oil", np.zeros(len(faces)))
        xr = np.c_[pts[:, 0], np.hypot(pts[:, 1], pts[:, 2])] * 1e3
        polys = [xr[f] for f in faces]
        polys += [p * np.array([1, -1]) for p in polys]
        pc = PolyCollection(polys, array=np.r_[a, a], cmap="Blues", clim=(0, 1), edgecolors="none")
        ax[0].add_collection(pc)
        fig.colorbar(pc, ax=ax[0], label="доля масла", shrink=0.8)
    ax[0].set_title("разрез через ось (доля масла)")
    # 2), 3) shadow pictures of the oil surface, sector copied around the axis
    iso = path / "iso.vtk"
    info = "масла ещё нет"
    if iso.exists():
        pts, faces, _ = read_vtk(iso)
        tri = triangles(faces)
        if len(tri):
            P = pts[tri]
            cen = P.mean(1)
            nA = 0.5 * np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])   # area vector
            n = max(1, round(360 / sector))
            X, Y, Z, WS, WF = [], [], [], [], []
            for k in range(n):
                c, s_ = math.cos(math.radians(k * sector)), math.sin(math.radians(k * sector))
                X.append(cen[:, 0])
                Y.append(cen[:, 1] * c - cen[:, 2] * s_)
                Z.append(cen[:, 1] * s_ + cen[:, 2] * c)
                WS.append(np.abs(nA[:, 1] * s_ + nA[:, 2] * c))     # seen from the side (along z)
                WF.append(np.abs(nA[:, 0]))                          # seen from downstream (along x)
            X, Y, Z = (np.concatenate(q) * 1e3 for q in (X, Y, Z))
            WS, WF = (np.concatenate(q) * 1e6 for q in (WS, WF))
            shadow(ax[1], X, Y, WS, [40, 250, -150, 150], (210, 300))       # 1 mm pixels
            shadow(ax[2], Z, Y, WF, [-150, 150, -150, 150], (300, 300))
            b = breakup(pts, tri)
            if b:
                info = f"плёнка цела до r ≈ {b[0] * 1e3:.0f} мм, оторвалось кусков (в четверти): {b[1]}"
    # colour bar always drawn (same layout in every frame, even before the oil arrives)
    sm = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(0, SHADOW_MAX), cmap="Greys")
    fig.colorbar(sm, ax=ax[2], label="поверхностей масла на луче", shrink=0.8)
    ax[1].set_title("масло сбоку (тень, 360°)")
    ax[2].set_title("масло со стороны выхода (тень, 360°)")
    for k, a_ in enumerate(ax):
        if k < 2:
            a_.add_patch(plt.Rectangle((DISK_X[0] * 1e3, -DISK_R * 1e3), (DISK_X[1] - DISK_X[0]) * 1e3,
                                       2 * DISK_R * 1e3, fc="0.45", zorder=3))
            a_.add_patch(plt.Rectangle((-50, -LANCE_R * 1e3), LANCE_X * 1e3 + 50, 2 * LANCE_R * 1e3,
                                       fc="0.45", zorder=3))
            a_.set_xlim(40, 250)
            a_.set_xlabel("x, мм (газ слева направо)")
        else:
            a_.add_patch(plt.Circle((0, 0), DISK_R * 1e3, fc="0.45", zorder=3))     # the disk hides the oil on its far side
            a_.set_xlim(-150, 150)
            a_.set_xlabel("z, мм")
        a_.set_ylim(-150, 150)
        a_.set_aspect("equal")
    ax[0].set_ylabel("r, мм")
    fig.suptitle(f"t = {t * 1e3:.1f} мс   ({info})")
    fig.savefig(OUT / f"f_{t * 1e3:07.3f}.png", dpi=90)
    plt.close(fig)
    return info


def main():
    if not SURF.is_dir():
        sys.exit("no postProcessing/surfaces yet: the run writes it every 0.1 ms")
    sector = SECTOR
    OUT.mkdir(exist_ok=True)
    times = sorted((float(d.name), d) for d in SURF.iterdir() if re.fullmatch(r"[0-9.eE+-]+", d.name))
    summary = []
    for t, d in times:
        png = OUT / f"f_{t * 1e3:07.3f}.png"
        if png.exists() and d is not times[-1][1]:
            continue
        info = frame(t, d, sector)
        rim = ""
        if (d / "slice.vtk").exists():
            rim = sheet_rim(*read_vtk(d / "slice.vtk"))
        line = f"t = {t * 1e3:5.2f} мс: {info}" + (f"; {rim}" if rim else "")
        print(line)
        summary.append(line)
    if summary:
        with open(CASE / "sheet3d_summary.txt", "a", encoding="utf-8") as f:
            f.write("\n".join(summary) + "\n")
    if "--gif" in sys.argv:
        try:
            from PIL import Image
        except ImportError:
            sys.exit("frames are drawn; for sheet3d.gif install Pillow: "
                     "pacman -S mingw-w64-x86_64-python-pillow")
        files = sorted(OUT.glob("f_*.png"))
        if len(files) > 1:
            w = 1100
            im = [Image.open(f).convert("RGB") for f in files]
            im = [i.resize((w, round(i.height * w / i.width)), Image.LANCZOS).quantize(64) for i in im]
            im[0].save(CASE / "sheet3d.gif", save_all=True, append_images=im[1:],
                       duration=[120] * (len(im) - 1) + [2000], loop=0, optimize=True)
            print(f"written sheet3d.gif ({len(im)} frames)")


if __name__ == "__main__":
    main()
