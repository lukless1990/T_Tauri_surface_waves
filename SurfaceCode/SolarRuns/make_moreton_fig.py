#!/usr/bin/env python3
"""make_moreton_fig.py — Appendix B figure: simulated Moreton front vs Francile+2013 Eq. (1).

A&A single-column vector PDF, rendered through the paper's shared style.
Usage: python3 make_moreton_fig.py [run_dir] [outfile.pdf]
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Paper", "figscripts"))
import aastyle
aastyle.setup()
from aastyle import COL, C
import matplotlib.pyplot as plt

C1, DELTA, C2 = 15.287, 0.578627, -108.609      # Francile et al. 2013, Eq. (1)
T_OBS = (86.0, 545.0)                            # 18:43:26-18:51:05 UT, s since t_i
francile = lambda t: C1 * np.asarray(t, float)**DELTA + C2

def track_snaps(run, thr=0.30):
    """Front position from the SNAPSHOTS — usable while a run is still in flight
    (the 1 s trace stacks are only written when the time loop completes)."""
    import glob
    x = np.loadtxt(f"{run}/grid_x_Mm.txt"); z = np.loadtxt(f"{run}/grid_z_Mm.txt")
    meta = dict(l.strip().split("=") for l in open(f"{run}/meta.txt"))
    xfl = float(meta["x_fl_Mm"]); jh = np.argmin(abs(z - float(meta["z_ha_Mm"])))
    ts = np.atleast_1d(np.loadtxt(f"{run}/snap_times.txt"))
    tt, dd = [], []
    for k, f in enumerate(sorted(glob.glob(f"{run}/snap_vz_*.txt")), start=1):
        v = np.loadtxt(f)[:, jh]
        i = np.where((np.abs(v) > thr) & (x > xfl + 30))[0]
        if len(i): tt.append(ts[k-1]); dd.append(x[i[-1]] - xfl)
    return np.array(tt), np.array(dd)

def track(run, thr=0.30):
    if not os.path.exists(f"{run}/trace_times.txt"):
        return track_snaps(run, thr)
    trt = np.loadtxt(f"{run}/trace_times.txt"); trx = np.loadtxt(f"{run}/trace_x_Mm.txt")
    meta = dict(l.strip().split("=") for l in open(f"{run}/meta.txt")); xfl = float(meta["x_fl_Mm"])
    st = np.loadtxt(f"{run}/trace_chromo_vz.txt")
    pos = np.full(len(trt), np.nan)
    for k in range(len(trt)):
        i = np.where((np.abs(st[k]) > thr) & (trx > xfl + 30))[0]
        if len(i): pos[k] = trx[i[-1]] - xfl
    kmax = int(np.nanargmax(pos)); keep = np.zeros(len(pos), bool); best = -np.inf
    for k in range(kmax + 1):
        if np.isfinite(pos[k]) and pos[k] > best: keep[k], best = True, pos[k]
    return trt[keep], pos[keep]

run = sys.argv[1] if len(sys.argv) > 1 else "run_krfull"
TMAX = float(sys.argv[3]) if len(sys.argv) > 3 else float(os.environ.get("MORETON_TMAX", "inf"))
out = sys.argv[2] if len(sys.argv) > 2 else "../Paper/Surface_Wave_paper/figs/figB_moreton.pdf"

fig, ax = plt.subplots(figsize=(COL, COL * 0.78))
tt = np.linspace(30, 620, 600)
sel = (tt >= T_OBS[0]) & (tt <= T_OBS[1])
ax.plot(tt, francile(tt), color="k", lw=0.7, ls=(0, (4, 2.5)), alpha=0.55, zorder=4)
ax.plot(tt[sel], francile(tt[sel]), color="k", lw=1.5, zorder=5, label="Francile+2013 fit")
t, d = track(run)
m = t <= TMAX; t, d = t[m], d[m]        # freeze the comparison at TMAX (obs coverage ends 545 s)
partial = not os.path.exists(f"{run}/trace_times.txt")
ax.plot(t, d, color=C["red"], lw=1.3, marker="o" if partial else None,
        ms=3, zorder=3, label="Model")
ax.set_xlabel("$t$ [s]"); ax.set_ylabel("$d$ [Mm]")
ax.set_xlim(0, 620); ax.set_ylim(0, 650)
ax.legend(loc="upper left", frameon=False)
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out); print("wrote", out, "from", run)
