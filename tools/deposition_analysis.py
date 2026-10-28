"""Where do the oil droplets hit the wall, and how hot does the oil get there?

Reads, from a finished spray run:
  postProcessing/lagrangian/cloud/outletParcels/<time>/{walls,nozzle,outlet}.post
        every droplet impact of one cloud evaluation (= 1 s of injection in the
        steady cloud): cell, diameter, temperature, number of real droplets
  processor*/<time>/C  (or <time>/C)   cell centres -> impact position
  VTK/walls/walls_<time>.vtk           wall gas temperature and shear stress

and writes deposition.png, deposition_bins.csv and a text summary.

Film temperature estimate (per wall bin, local energy balance, no film flow,
no evaporation, no conduction in the steel):
    m'' cp_oil (T_film - T_impact) = h (T_gas - T_film)
    h from the Chilton-Colburn analogy: h = cp_gas tau_w / (U_bulk Pr^(2/3))
Heavily wetted bins stay near the droplet temperature; sparsely or occasionally
wetted bins approach the gas temperature: thin hot oil films are where coke forms.

Usage (from the case directory, after ./deposition prepared the inputs):
    python ../tools/deposition_analysis.py [case] [--samples 4]
"""
import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "geometry"))
import params as geo  # noqa: E402

CP_OIL = 2100.0          # J/(kg K), as constantProperties in cloudProperties
CP_GAS, PR_GAS = 2835.0, 1.04
MDOT_GAS, RHO_GAS = 21.11, 1.106
L_EXT = 2000.0           # outlet extension in mesh/make_blockmesh.py, mm
DX, NTHETA = 100.0, 24   # bin size: 100 mm along the axis (> axial cell size), 15 deg around

# film temperature bands for heavy oil (rough guide; real limits depend on the
# oil: Conradson carbon, asphaltenes, residence time of the film)
BANDS = [(300.0, "низкий"), (380.0, "умеренный"), (1e9, "высокий")]


def wall_radius(x_mm):
    return geo.wall_radius(min(x_mm, geo.X_END))


# ---------------------------------------------------------------- readers
def read_cell_centres(path):
    """internalField of a volVectorField C, ascii or binary."""
    raw = path.read_bytes()
    i = raw.index(b"internalField")
    m = re.compile(rb"nonuniform\s+List<vector>\s*(\d+)\s*\(").search(raw, i)
    n, start = int(m.group(1)), m.end()
    if b"format      binary" in raw[:i] or b"format binary" in raw[:i]:
        return np.frombuffer(raw[start:start + 24 * n], dtype="<f8").reshape(n, 3)
    vecs = re.findall(r"\(([^()]+)\)", raw[start:].decode(errors="ignore"))[:n]
    return np.array([list(map(float, v.split())) for v in vecs])


def load_centres(case):
    procs = sorted(case.glob("processor[0-9]*"), key=lambda p: int(p.name[9:]))
    dirs = procs if procs else [case]
    centres = []
    for d in dirs:
        files = sorted(d.glob("[0-9]*/C"), key=lambda f: float(f.parent.name))
        if not files:
            sys.exit(f"no cell centres in {d}: run ./deposition first")
        centres.append(read_cell_centres(files[-1]))
    return centres


def read_post(path):
    """Impact records of one .post file -> dict of arrays."""
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        return None
    names = lines[0].lstrip("#").replace("(", " ").replace(")", " ").split()
    rows = [ln.replace("(", " ").replace(")", " ").split() for ln in lines[1:] if ln.strip()]
    rows = [r for r in rows if len(r) == len(names)]
    if not rows:
        return None
    col = {n: k for k, n in enumerate(names)}

    def get(n, typ=float):
        return np.array([typ(r[col[n]]) for r in rows])

    d = get("d")
    return {
        "proc": get("currentProc", int), "cell": get("celli", int), "d": d,
        "T": get("T"), "mass": get("nParticle") * get("rho") * math.pi / 6 * d**3,
    }


def read_vtk_walls(path):
    """Face centres, areas and cell data (T, wallShearStress) of a legacy ascii VTK patch."""
    tok = path.read_text().split()
    i = tok.index("POINTS")
    n = int(tok[i + 1])
    P = np.array(tok[i + 3:i + 3 + 3 * n], float).reshape(n, 3)
    i = tok.index("POLYGONS")
    nf, nt = int(tok[i + 1]), int(tok[i + 2])
    conn = list(map(int, tok[i + 3:i + 3 + nt]))
    C, A, k = np.zeros((nf, 3)), np.zeros(nf), 0
    for f in range(nf):
        q = conn[k + 1:k + 1 + conn[k]]
        k += 1 + conn[k]
        p = P[q]
        C[f] = p.mean(0)
        a = np.zeros(3)
        for j in range(1, len(q) - 1):
            a += np.cross(p[j] - p[0], p[j + 1] - p[0])
        A[f] = 0.5 * np.linalg.norm(a)
    data = {}
    for name, nc in (("T", 1), ("wallShearStress", 3)):
        for j in range(len(tok) - 3):
            if tok[j] == name and tok[j + 1] == str(nc) and tok[j + 2] == str(nf):
                data[name] = np.array(tok[j + 4:j + 4 + nc * nf], float).reshape(nf, nc).squeeze()
                break
    return C, A, data


# ---------------------------------------------------------------- analysis
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", nargs="?", default=".")
    ap.add_argument("--samples", type=int, default=4,
                    help="use the last N written cloud evaluations (skip the start-up)")
    args = ap.parse_args()
    case = Path(args.case)

    post_root = case / "postProcessing" / "lagrangian" / "cloud" / "outletParcels"
    times = sorted((d for d in post_root.glob("*") if d.is_dir()), key=lambda d: float(d.name))
    times = times[-args.samples:]
    if not times:
        sys.exit("no droplet impact data: run ./Allrun.spray first")
    centres = load_centres(case)

    hits = {}
    for patch in ("walls", "nozzle", "outlet"):
        recs = [r for r in (read_post(t / f"{patch}.post") for t in times if (t / f"{patch}.post").exists()) if r]
        if recs:
            h = {k: np.concatenate([r[k] for r in recs]) for k in recs[0]}
            h["xyz"] = np.array([centres[p][c] for p, c in zip(h["proc"], h["cell"])]) * 1e3  # mm
            hits[patch] = h
    nS = len(times)
    total = sum(h["mass"].sum() for h in hits.values()) / nS

    vtk = sorted((case / "VTK" / "walls").glob("walls_*.vtk"))
    if not vtk:
        sys.exit("no VTK/walls/*.vtk: run ./deposition first")
    Cw, Aw, wdat = read_vtk_walls(vtk[-1])
    Cw = Cw * 1e3

    # ---- bins on the unrolled wall: x along the axis, theta from the top (lance side)
    x_end = geo.X_END + L_EXT
    # 100 mm bins in the apparatus, 250 mm on the coarse outlet extension
    xb = np.r_[np.arange(0, geo.X_END, DX), np.arange(geo.X_END, x_end + 250.0, 250.0)]
    tb = np.linspace(-180, 180, NTHETA + 1)
    nx, nt = len(xb) - 1, NTHETA

    def theta(xyz):
        return np.degrees(np.arctan2(xyz[:, 2], xyz[:, 1]))

    def bin_of(xyz):
        ix = np.clip(np.digitize(xyz[:, 0], xb) - 1, 0, nx - 1)
        it = np.clip(np.digitize(theta(xyz), tb) - 1, 0, nt - 1)
        return ix, it

    area = np.zeros((nx, nt))
    Tg = np.zeros((nx, nt))
    tau = np.zeros((nx, nt))
    ix, it = bin_of(Cw)
    np.add.at(area, (ix, it), Aw)
    np.add.at(Tg, (ix, it), wdat["T"] * Aw)
    np.add.at(tau, (ix, it), np.linalg.norm(wdat["wallShearStress"], axis=1) * Aw)
    ok = area > 0
    Tg[ok] /= area[ok]
    tau[ok] /= area[ok]

    w = hits.get("walls")
    mflux = np.zeros((nx, nt))
    Ti = np.zeros((nx, nt))
    if w is not None:
        ix, it = bin_of(w["xyz"])
        np.add.at(mflux, (ix, it), w["mass"])
        np.add.at(Ti, (ix, it), w["mass"] * w["T"])
        wet = mflux > 0
        Ti[wet] /= mflux[wet]
        mflux[ok] /= area[ok] * nS                     # kg/(m2 s)

    xc = 0.5 * (xb[1:] + xb[:-1])
    Ubulk = MDOT_GAS / RHO_GAS / (math.pi * (np.array([wall_radius(x) for x in xc]) * 1e-3) ** 2)
    hcoef = CP_GAS * tau / (Ubulk[:, None] * PR_GAS ** (2 / 3))
    wet = (mflux > 0) & ok
    Tf = np.full((nx, nt), np.nan)
    Tf[wet] = (mflux[wet] * CP_OIL * Ti[wet] + hcoef[wet] * Tg[wet]) / (mflux[wet] * CP_OIL + hcoef[wet])

    # ---- summary
    out = []
    p = lambda s="": out.append(s)  # noqa: E731
    p(f"Разобрано {nS} вычислений облака (итерации {', '.join(t.name for t in times)}), "
      f"масса на одно вычисление {total:.2f} кг/с")
    for patch, h in hits.items():
        m = h["mass"].sum() / nS
        Tm = (h["mass"] * h["T"]).sum() / h["mass"].sum() - 273.15
        dm = (h["mass"] * h["d"]).sum() / h["mass"].sum() * 1e6
        p(f"  {patch:7s}: {m:6.2f} кг/с ({100 * m / total:5.1f} %), капли в момент удара: "
          f"T = {Tm:5.0f} °C, d = {dm:4.0f} мкм (средние по массе)")
    if "nozzle" in hits:
        h = hits["nozzle"]
        disk = h["xyz"][:, 0] >= geo.X_DISK_BACK - 5
        p(f"  из них на диск {100 * h['mass'][disk].sum() / h['mass'].sum():.0f} %, "
          f"на лансет и корпус {100 * h['mass'][~disk].sum() / h['mass'].sum():.0f} %")
    p()
    if w is not None:
        a_wet = area[wet].sum()
        p(f"Смоченная площадь стенки: {a_wet:.2f} м² из {area[ok].sum():.2f} м²")
        lo = 0.0
        for hi, name in BANDS:
            sel = wet & (Tf - 273.15 >= lo) & (Tf - 273.15 < hi)
            p(f"  T плёнки {lo:3.0f}–{hi if hi < 1e8 else 999:3.0f} °C (риск {name:9s}): "
              f"{area[sel].sum():5.2f} м², {100 * (mflux[sel] * area[sel]).sum() / max((mflux[wet] * area[wet]).sum(), 1e-12):5.1f} % масла")
            lo = hi
        real = xc <= geo.X_END
        m_all = (mflux * area).sum()
        p(f"  из осевшего на обечайку: в аппарате (до конца диффузора) "
          f"{100 * (mflux[real] * area[real]).sum() / max(m_all, 1e-12):.0f} %, "
          f"на расчётном удлинении за диффузором {100 * (mflux[~real] * area[~real]).sum() / max(m_all, 1e-12):.0f} % "
          f"(в реальности это следующий аппарат/газоход)")
        xs = [x for x in xc[wet.any(1)]]
        p(f"Масло попадает на стенку на участке x = {min(xs):.0f}…{max(xs):.0f} мм "
          f"(конфузор {geo.X_CONV:.0f}…{geo.X_THROAT:.0f}, горловина до {geo.X_DIFF:.0f})")
        p(f"Теплоотдача газ→стенка (по аналогии Кольборна): {np.nanmin(hcoef[ok]):.0f}…{np.nanmax(hcoef[ok]):.0f} Вт/(м²·К)")
    summary = "\n".join(out)
    print(summary)
    (case / "deposition_summary.txt").write_text(summary + "\n")

    with open(case / "deposition_bins.csv", "w") as f:
        f.write("x_mm,theta_deg,area_m2,oil_flux_kg_m2s,T_impact_C,T_gas_C,h_W_m2K,T_film_C\n")
        for i in range(nx):
            for j in range(nt):
                if ok[i, j]:
                    f.write(f"{xc[i]:.0f},{0.5 * (tb[j] + tb[j + 1]):.1f},{area[i, j]:.5f},{mflux[i, j]:.5g},"
                            f"{Ti[i, j] - 273.15 if wet[i, j] else float('nan'):.1f},{Tg[i, j] - 273.15:.1f},"
                            f"{hcoef[i, j]:.1f},{Tf[i, j] - 273.15 if wet[i, j] else float('nan'):.1f}\n")

    plot(case, xb, tb, mflux, Tf, wet, ok, hits, nS, xc, area, summary)


def plot(case, xb, tb, mflux, Tf, wet, ok, hits, nS, xc, area, summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, ListedColormap, BoundaryNorm

    # show only where something happens: from 300 mm before the nozzle to the end of the throat region
    x0 = max(0.0, geo.X_STEM - 300)
    xs_wet = xc[wet.any(1)]
    x1 = max(geo.X_DIFF + 300, (xs_wet.max() + 100) if len(xs_wet) else 0)
    sel = (xb[:-1] >= x0) & (xb[1:] <= x1)
    X = xb[np.r_[sel, False] | np.r_[False, sel]]

    fig = plt.figure(figsize=(15, 10.5))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.9])
    marks = [(geo.X_DISK_FRONT, "диск"), (geo.X_THROAT, "горловина"), (geo.X_DIFF, "диффузор"),
             (geo.X_END, "конец диффузора")]

    def deco(a):
        for x, name in marks:
            if x0 <= x <= x1:
                a.axvline(x, color="k", lw=0.7, ls="--")
                a.text(x + 15, 170, name, fontsize=8, ha="left", va="top", rotation=90,
                       bbox=dict(fc="white", ec="none", pad=1, alpha=0.8))
        a.set_ylabel("угол от верха, ° (0 = сторона лансета)")
        a.set_xlabel("x вдоль оси, мм")
        a.set_ylim(-180, 180)
        a.set_yticks(range(-180, 181, 60))

    a = fig.add_subplot(gs[0, :])
    Z = np.ma.masked_where(~wet[sel], mflux[sel])
    pc = a.pcolormesh(X, tb, Z.T, cmap="viridis", norm=LogNorm(vmin=max(Z.min(), 1e-3), vmax=Z.max()), shading="flat")
    a.pcolormesh(X, tb, np.ma.masked_where(wet[sel] | ~ok[sel], np.ones_like(mflux[sel])).T,
                 cmap=ListedColormap(["#e9ecef"]), shading="flat")
    fig.colorbar(pc, ax=a, label="поток масла на стенку, кг/(м²·с)")
    a.set_title("Где масло попадает на стенку (развёртка обечайки; серое: капли не долетают)", pad=4)
    deco(a)

    a = fig.add_subplot(gs[1, :])
    bounds = [0, 200, 250, 300, 340, 380, 420, 500]
    cmap = plt.get_cmap("RdYlGn_r", len(bounds) - 1)
    Z = np.ma.masked_where(~wet[sel], Tf[sel] - 273.15)
    pc = a.pcolormesh(X, tb, Z.T, cmap=cmap, norm=BoundaryNorm(bounds, cmap.N), shading="flat")
    a.pcolormesh(X, tb, np.ma.masked_where(wet[sel] | ~ok[sel], np.ones_like(mflux[sel])).T,
                 cmap=ListedColormap(["#e9ecef"]), shading="flat")
    fig.colorbar(pc, ax=a, label="оценка T масляной плёнки, °C")
    a.set_title("Оценка T масла на стенке: зелёное = плёнку охлаждает свежее масло, красное = редкие капли на горячей стенке (кокс)", pad=4)
    deco(a)

    a = fig.add_subplot(gs[2, 0])
    w = hits.get("walls")
    if w is not None:
        share = (mflux * area).sum(1) * 100 / max((mflux * area).sum(), 1e-12)
        width = (xb[1:] - xb[:-1])[sel]
        a.bar(xc[sel], share[sel] * DX / width, width=width, color="#1c7ed6")
        a2 = a.twinx()
        r = [geo.wall_radius(min(x, geo.X_END)) for x in xc[sel]]
        a2.plot(xc[sel], r, color="gray", lw=1)
        a2.set_ylabel("радиус стенки, мм", color="gray")
        a2.set_ylim(0, max(r) * 1.1)
    a.set_xlabel("x вдоль оси, мм")
    a.set_ylabel("% масла на стенке на 100 мм длины")
    a.set_title("Доля осевшего масла вдоль оси")

    a = fig.add_subplot(gs[2, 1])
    for patch, color in (("walls", "#1c7ed6"), ("nozzle", "#e8590c")):
        if patch in hits:
            h = hits[patch]
            a.hist(h["T"] - 273.15, bins=30, weights=h["mass"] / nS, alpha=0.6, color=color,
                   label={"walls": "на обечайку", "nozzle": "на диск/лансет"}[patch])
    a.set_xlabel("температура капли в момент удара, °C")
    a.set_ylabel("кг/с")
    a.legend()
    a.set_title("Насколько капли успевают нагреться до удара")

    fig.suptitle("Осаждение масла и риск коксования (оценка; плёнка не моделируется)", fontsize=13)
    fig.tight_layout()
    fig.savefig(case / "deposition.png", dpi=110)
    print(f"written {case / 'deposition.png'}, deposition_bins.csv, deposition_summary.txt")


if __name__ == "__main__":
    main()
