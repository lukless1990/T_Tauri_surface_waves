#!/usr/bin/env python3
# make_solar_fig.py — appendix figure: SurfaceCode solar Moreton-wave benchmark vs the
# observed 2006 Dec 6 event (Balasubramaniam et al. 2010). Distance–time of the coronal
# front for the impulsive (production) and sustained-piston (extreme) drivers, with the
# constant-deceleration fit and the observed kinematic fit.
import os, json
import numpy as np
import aastyle
aastyle.setup()
from aastyle import COL, C
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
sr   = os.path.join(here, "..", "..", "SolarRuns")
out  = os.path.join(here, "..", "Surface_Wave_paper", "figs")
nums = {}

def load(run):
    d = os.path.join(sr, run)
    meta = {}
    for ln in open(os.path.join(d, "meta.txt")):
        k, v = ln.strip().split("=")
        try: meta[k] = float(v)
        except ValueError: meta[k] = v
    trt = np.loadtxt(os.path.join(d, "trace_times.txt"))
    trx = np.loadtxt(os.path.join(d, "trace_x_Mm.txt"))
    trdp = np.loadtxt(os.path.join(d, "trace_corona_dpp.txt"))
    return meta, trt, trx, trdp

def front(trx, stack, thresh, x_fl):
    pos = np.full(stack.shape[0], np.nan)
    for k in range(stack.shape[0]):
        idx = np.where((np.abs(stack[k]) > thresh) & (trx > x_fl + 10))[0]
        if len(idx): pos[k] = trx[idx[-1]]
    return pos

def mono_track(trt, pos, x_fl, xmax):
    k_end = int(np.nanargmax(pos))
    mono = np.zeros_like(pos, bool); best = -np.inf
    for k in range(k_end + 1):
        if np.isfinite(pos[k]) and pos[k] > best:
            mono[k] = True; best = pos[k]
    msk = mono & (pos > x_fl + 30) & (pos < xmax - 20)
    return trt[msk], pos[msk] - x_fl

fig, ax = plt.subplots(figsize=(COL, 2.7))
runs = [("run_production", C["blue"], "impulsive blast (production)"),
        ("run_extreme",    C["orange"], "sustained CME-piston driver")]
for run, cc, lab in runs:
    meta, trt, trx, trdp = load(run)
    x_fl = meta["x_fl_Mm"]
    pos = front(trx, trdp, 0.03, x_fl)
    tf, df = mono_track(trt, pos, x_fl, trx[-1])
    ax.plot(tf, df, ".", color=cc, ms=2.2, label=lab)
    if run == "run_production":
        c2, c1, c0 = np.polyfit(tf, df, 2)
        tt = np.linspace(tf[0], tf[-1], 100)
        ax.plot(tt, c0 + c1*tt + c2*tt**2, "-", color="k", lw=0.9,
                label="constant-deceleration fit")
        v_at = lambda t: (c1 + 2*c2*t)*1e3
        nums.update(v0=v_at(tf[0]), v_end=v_at(tf[-1]), a=2*c2*1e3,
                    v_mean=(df[-1]-df[0])/(tf[-1]-tf[0])*1e3, d_max_prod=df[-1])
    else:
        nums.update(d_max_ext=np.nanmax(df))
# observed kinematic fit (Balasubramaniam et al. 2010: v0=1125 km/s, a=-0.87 km/s^2),
# aligned to launch at the production fit's origin
t_obs = np.linspace(0, 700, 200)
d_obs = 1.125*t_obs + 0.5*(-0.87e-3)*t_obs**2
d_obs[np.gradient(d_obs, t_obs) < 0] = np.nan
ax.plot(t_obs, d_obs, "--", color=C["gray"], lw=1.1,
        label="observed 2006 Dec 6 fit")
ax.axhline(550, color=C["gray"], lw=0.6, ls=":")
ax.text(15, 562, "observed visibility limit (~550 Mm)", fontsize=6.2, color="0.35")
ax.set_xlabel("$t$ [s]"); ax.set_ylabel("distance from source [Mm]")
ax.set_xlim(0, 720); ax.set_ylim(0, 620)
ax.legend(loc="lower right", fontsize=6.2)
fig.savefig(os.path.join(out, "figB_moreton.pdf"))
print(json.dumps({k: round(float(v), 2) for k, v in nums.items()}, indent=0))
