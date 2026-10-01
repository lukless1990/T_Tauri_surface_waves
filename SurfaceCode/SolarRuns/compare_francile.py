#!/usr/bin/env python3
"""compare_francile.py — simulated Moreton wavefronts vs Francile et al. 2013 Eq. (1).

Francile+2013 (A&A 552, A3) Eq. 1, the power-law fit to the 2D-averaged chromospheric
distance d from the radiant point Q0:
     d(t) = c1 (t - t_i)^delta + c2,   delta=0.578627, c1=15.287, c2=-108.609
with t_i = 18:42:00 UT, d in Mm and (t - t_i) in seconds. Observed coverage is
18:43:26-18:51:05 UT, i.e. t - t_i = 86-545 s (d = 93-477 Mm).

Time alignment: simulation t=0 is flare ignition, mapped directly onto (t - t_i).
Krause+2015 set h=35 Mm precisely so this delay reproduces the observed ~100 s between
flare ignition and the emergence of the chromospheric disturbance.

Usage: python3 compare_francile.py [out.png]
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C1, DELTA, C2 = 15.287, 0.578627, -108.609
T_OBS = (86.0, 545.0)                       # observed coverage, s since t_i
def francile(t):
    t = np.asarray(t, float)
    return np.where(t > 0, C1*np.maximum(t, 1e-9)**DELTA + C2, np.nan)

def track(run, which="chromo", thr=None):
    """Front distance vs time from the 1 s trace stacks; monotone advancing part only."""
    d = run
    trt = np.loadtxt(f"{d}/trace_times.txt")
    trx = np.loadtxt(f"{d}/trace_x_Mm.txt")
    meta = dict(l.strip().split("=") for l in open(f"{d}/meta.txt"))
    xfl = float(meta["x_fl_Mm"])
    if which == "chromo":
        stack, thr = np.loadtxt(f"{d}/trace_chromo_vz.txt"), (0.30 if thr is None else thr)
    else:
        stack, thr = np.loadtxt(f"{d}/trace_corona_dpp.txt"), (0.13 if thr is None else thr)
    pos = np.full(len(trt), np.nan)
    for k in range(len(trt)):
        i = np.where((np.abs(stack[k]) > thr) & (trx > xfl + 30))[0]
        if len(i):
            pos[k] = trx[i[-1]] - xfl
    kmax = int(np.nanargmax(pos))
    keep, best = np.zeros(len(pos), bool), -np.inf
    for k in range(kmax + 1):
        if np.isfinite(pos[k]) and pos[k] > best:
            keep[k], best = True, pos[k]
    return trt[keep], pos[keep]

fig, ax = plt.subplots(figsize=(6.0, 4.2))

# Francile+2013 Eq. (1): solid over the observed coverage, thin dashed where extrapolated
tt = np.linspace(30, 620, 600)
sel = (tt >= T_OBS[0]) & (tt <= T_OBS[1])
ax.plot(tt, francile(tt), color="k", lw=0.9, ls=(0, (5, 3)), alpha=0.5, zorder=4)
ax.plot(tt[sel], francile(tt[sel]), color="k", lw=2.2, zorder=5, label="Francile+2013 fit")

t, d = track("run_krfull", "chromo")
ax.plot(t, d, color="#c1272d", lw=1.8, zorder=3, label="Model")

ax.set_xlabel("$t$  [s]")
ax.set_ylabel("$d$  [Mm]")
ax.set_xlim(0, 620); ax.set_ylim(0, 650)
ax.grid(alpha=0.25, lw=0.5)
ax.legend(loc="upper left", fontsize=9, frameon=False)
fig.tight_layout()
out = sys.argv[1] if len(sys.argv) > 1 else "plots_krfull/fig_francile_eq1.png"
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=170)
print("wrote", out)
for tq in (100, 200, 300, 400, 500, 545):
    print(f"  Eq.(1) at t-ti={tq:4d} s : d={francile(tq):6.1f} Mm")
