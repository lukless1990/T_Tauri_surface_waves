#!/usr/bin/env python3
"""compare_drivers.py — the referee test, visualized.

Overlays the front kinematics of the impulsive run (run_production) and the
sustained CME-piston / large-domain run (run_extreme) against the observed
2006 Dec 6 curve (Balasubramaniam et al. 2010). Same solver, same dissipation
physics, same resolution — only the driver and domain change. If the visibility
range stretches toward the observed ~550 Mm under sustained driving, the short
range of the impulsive run is DRIVER-limited, not dissipation-limited.

Usage: python3 compare_drivers.py [run_production run_extreme]
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

runs = sys.argv[1:] or ["run_production", "run_extreme"]

def load(rundir):
    rundir = rundir.rstrip("/")
    meta = {}
    with open(os.path.join(rundir, "meta.txt")) as f:
        for ln in f:
            k, v = ln.strip().split("=")
            try: meta[k] = float(v)
            except ValueError: meta[k] = v
    d = dict(
        label=os.path.basename(rundir).replace("run_", ""), meta=meta,
        x_fl=meta["x_fl_Mm"],
        t=np.loadtxt(os.path.join(rundir, "trace_times.txt")),
        x=np.loadtxt(os.path.join(rundir, "trace_x_Mm.txt")),
        vz=np.loadtxt(os.path.join(rundir, "trace_chromo_vz.txt")),
        dp=np.loadtxt(os.path.join(rundir, "trace_corona_dpp.txt")),
    )
    return d

def track(d, stack_key, thresh):
    stack = d[stack_key]; x = d["x"]; x_fl = d["x_fl"]
    pos = np.full(stack.shape[0], np.nan)
    for k in range(stack.shape[0]):
        idx = np.where((np.abs(stack[k]) > thresh) & (x > x_fl + 30))[0]
        if len(idx): pos[k] = x[idx[-1]]
    # monotone advancing segment only
    if np.all(np.isnan(pos)): return d["t"], pos, np.nan
    k_end = int(np.nanargmax(pos)); dmax = pos[k_end] - x_fl
    mono = np.full_like(pos, np.nan); best = -np.inf
    for k in range(k_end + 1):
        if np.isfinite(pos[k]) and pos[k] > best:
            mono[k] = pos[k]; best = pos[k]
    return d["t"], mono - x_fl, dmax

OBS = dict(v0=1125., a=-0.87, v_mean=850., d_last=550.)

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
colors = {"production": "C0", "extreme": "C1"}
names  = {"production": "impulsive blast (A=100, 210 Mm box)",
          "extreme":    "sustained CME-piston (A=150, 600 Mm box)"}
summary = []
for r in runs:
    d = load(r); lab = d["label"]; c = colors.get(lab, "C2")
    t, pos, dmax = track(d, "dp", 0.03)
    ax[0].plot(t, pos, ".", color=c, ms=3, label=names.get(lab, lab) + " — coronal front")
    # fit the clean advancing part beyond the guard zone
    m = np.isfinite(pos) & (pos > 30) & (pos < (d["x"][-1] - d["x_fl"]) - 20)
    vmean = vfit = a = np.nan
    if m.sum() > 5:
        c2, c1, c0 = np.polyfit(t[m], pos[m], 2)
        vmean = (pos[m][-1] - pos[m][0]) / (t[m][-1] - t[m][0]) * 1e3
        vfit = c1 * 1e3; a = 2 * c2 * 1e3
        tt = np.linspace(t[m][0], t[m][-1], 80)
        ax[0].plot(tt, c0 + c1 * tt + c2 * tt**2, "-", color=c, lw=1.0)
    summary.append((lab, vmean, a, dmax))

t_obs = np.linspace(0, 700, 120)
d_obs = OBS["v0"] * 1e-3 * t_obs + 0.5 * OBS["a"] * 1e-3 * t_obs**2
ax[0].plot(t_obs, d_obs, "--", color="gray", lw=2,
           label=f"observed 2006 Dec 6 (v$_0$={OBS['v0']:.0f}, a={OBS['a']} km/s$^2$)")
ax[0].axhline(OBS["d_last"], color="gray", ls=":", lw=1)
ax[0].text(10, OBS["d_last"] + 8, "observed last-seen ~550 Mm", color="gray", fontsize=8)
ax[0].set_xlabel("t [s]"); ax[0].set_ylabel("front distance from source [Mm]")
ax[0].set_title("same solver, same dissipation physics — only the driver changes")
ax[0].legend(fontsize=8, loc="upper left"); ax[0].grid(alpha=0.3)
ax[0].set_ylim(0, 620)

ax[1].axis("off")
cells = [[f"{v:.0f}" if np.isfinite(v) else "—",
          f"{a:+.2f}" if np.isfinite(a) else "—",
          f"{dm:.0f}"] for (_, v, a, dm) in summary]
rowlab = [names.get(l, l).split(" (")[0] for (l, *_ ) in summary]
cells.append([f"{OBS['v_mean']:.0f}", f"{OBS['a']:+.2f}", f"~{OBS['d_last']:.0f}"])
rowlab.append("observed 2006 Dec 6")
tbl = ax[1].table(cellText=cells, rowLabels=rowlab,
                  colLabels=["mean speed\n[km/s]", "deceleration\n[km/s$^2$]",
                             "front reached\n[Mm]"],
                  loc="center", cellLoc="center")
tbl.auto_set_font_size(False); tbl.set_fontsize(9); tbl.scale(0.85, 2.0)
ax[1].set_title("driver-limited, not dissipation-limited", fontsize=11)
fig.tight_layout()
out = os.path.join(os.path.dirname(runs[0]) or ".", "plots_extreme")
os.makedirs(out, exist_ok=True)
fig.savefig(os.path.join(out, "fig6_driver_test.png"), dpi=150); plt.close(fig)
for (l, v, a, dm) in summary:
    print(f"  {l:12s}  mean {v:6.0f} km/s   decel {a:+.2f} km/s^2   front reached {dm:.0f} Mm")
print(f"  observed      mean {OBS['v_mean']:6.0f} km/s   decel {OBS['a']:+.2f} km/s^2   last seen ~{OBS['d_last']:.0f} Mm")
print(f"figure -> {out}/fig6_driver_test.png")
