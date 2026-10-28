"""Convergence plots for the Venturi case: residuals, pressure drop, outlet temperature, mass flow.

Usage (from the case directory or with the path to it):
    python ../tools/plot_convergence.py            # write convergence.png once
    python ../tools/plot_convergence.py --watch 30 # live window, refresh every 30 s
Needs: pip install matplotlib
"""
import argparse
import math
import time
from pathlib import Path

import matplotlib


def read_dat(name, case):
    """Read postProcessing/<name>/<time>/*.dat from all restart folders -> (columns, rows)."""
    root = case / "postProcessing" / name
    if not root.is_dir():
        # #includeFunc without name= writes to e.g. "residuals(p,U,h,k,omega)"
        alt = sorted((case / "postProcessing").glob(name + "(*")) if (case / "postProcessing").is_dir() else []
        if not alt:
            return [], []
        root = alt[0]
    dirs = sorted((d for d in root.iterdir() if d.is_dir()), key=lambda d: float(d.name))
    cols, rows = [], []
    for d in dirs:
        for f in sorted(d.glob("*.dat")):
            header = None
            for line in f.read_text(errors="ignore").splitlines():
                if line.startswith("#"):
                    header = line.lstrip("#").split()
                    continue
                vals = []
                for tok in line.split():
                    try:
                        vals.append(float(tok.strip("()")))
                    except ValueError:
                        vals.append(math.nan)
                if vals:
                    rows.append(vals)
            if header:
                cols = header
    # a newer restart overwrites the tail of the previous one
    rows.sort(key=lambda r: r[0])
    dedup = {}
    for r in rows:
        dedup[r[0]] = r
    return cols, list(dedup.values())


def column(cols, rows, key):
    i = next((k for k, c in enumerate(cols) if key in c), None)
    if i is None:
        return [], []
    xs = [r[0] for r in rows if len(r) > i]
    ys = [r[i] for r in rows if len(r) > i]
    return xs, ys


def smooth(a, xs, ys, color, n=100):
    """Thin raw curve plus a running mean over n iterations (spray results are noisy)."""
    a.plot(xs, ys, color=color, lw=0.6, alpha=0.45)
    if len(ys) > n:
        mean = [sum(ys[max(0, i - n + 1):i + 1]) / len(ys[max(0, i - n + 1):i + 1]) for i in range(len(ys))]
        a.plot(xs, mean, color=color, lw=2, label=f"среднее за {n} итераций")
        return mean[-1]
    return ys[-1]


def zoom(a, ys):
    """Set y limits from the last 90 % of the data, so start-up spikes do not flatten the curve."""
    tail = [y for y in ys[len(ys) // 10:] if not math.isnan(y)]
    if len(tail) > 5:
        lo, hi = min(tail), max(tail)
        pad = max(0.1 * (hi - lo), 0.01 * abs(hi), 1e-3)
        a.set_ylim(lo - pad, hi + pad)


def draw(fig, case):
    fig.clf()
    ax = fig.subplots(2, 2)
    fig.suptitle(f"Сходимость: {case.resolve().name}")

    a = ax[0][0]
    cols, rows = read_dat("residuals", case)
    for k, c in enumerate(cols[1:], start=1):
        ys = [r[k] if len(r) > k else math.nan for r in rows]
        if any(y > 0 for y in ys if not math.isnan(y)):
            a.semilogy([r[0] for r in rows], ys, label=c, lw=1)
    a.set_title("Невязки (начальные)")
    a.set_xlabel("итерация")
    a.grid(True, which="both", alpha=0.3)
    if cols:
        a.legend(fontsize=8)

    ci, ri = read_dat("inletAverage", case)
    co, ro = read_dat("outletAverage", case)
    xi, pi = column(ci, ri, "(p)")
    xo, po = column(co, ro, "(p)")
    a = ax[0][1]
    po_map = dict(zip(xo, po))
    dp = [(x, p - po_map[x]) for x, p in zip(xi, pi) if x in po_map]
    if dp:
        last = smooth(a, [d[0] for d in dp], [d[1] / 1000 for d in dp], "#1c7ed6")
        zoom(a, [d[1] / 1000 for d in dp])
        a.set_title(f"Перепад давления вход-выход: {last:.2f} кПа")
    else:
        a.set_title("Перепад давления (нет данных)")
    a.set_xlabel("итерация"); a.set_ylabel("кПа"); a.grid(True, alpha=0.3)

    a = ax[1][0]
    xt, tout = column(co, ro, "(T)")
    if tout:
        last = smooth(a, xt, [t - 273.15 for t in tout], "#e8590c")
        zoom(a, [t - 273.15 for t in tout])
        a.set_title(f"Средняя T на выходе: {last:.1f} °C")
    else:
        a.set_title("T на выходе (нет данных)")
    a.set_xlabel("итерация"); a.set_ylabel("°C"); a.grid(True, alpha=0.3)

    a = ax[1][1]
    cm, rm = read_dat("outletMassFlow", case)
    xm, mf = column(cm, rm, "(phi)")
    if mf:
        smooth(a, xm, [abs(m) for m in mf], "#2f9e44")
        zoom(a, [abs(m) for m in mf] + [21.11])
        a.axhline(21.11, color="gray", ls="--", lw=1, label="вход 21.11 кг/с")
        a.set_title(f"Расход газа на выходе: {abs(mf[-1]):.2f} кг/с")
        a.legend(fontsize=8)
    else:
        a.set_title("Расход на выходе (нет данных)")
    a.set_xlabel("итерация"); a.set_ylabel("кг/с"); a.grid(True, alpha=0.3)
    fig.tight_layout()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", nargs="?", default=".")
    ap.add_argument("--watch", type=float, default=0, help="refresh period, s (opens a window)")
    args = ap.parse_args()
    case = Path(args.case)
    if not args.watch:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(13, 8))
    if args.watch:
        try:
            fig.canvas.manager.set_window_title(f"Сходимость: {case.resolve().name}")
        except Exception:
            pass
        plt.ion()
        plt.show(block=False)
        while plt.fignum_exists(fig.number):
            try:
                draw(fig, case)
                fig.savefig(case / "convergence.png", dpi=100)
                fig.canvas.draw_idle()
                fig.canvas.flush_events()
            except Exception as e:  # e.g. a .dat file being written right now: try again next time
                print(f"{time.strftime('%H:%M:%S')} skipped one refresh: {e!r}", flush=True)
            plt.pause(args.watch)
    else:
        draw(fig, case)
        fig.savefig(case / "convergence.png", dpi=100)
        print(f"written {case / 'convergence.png'}")


if __name__ == "__main__":
    main()
