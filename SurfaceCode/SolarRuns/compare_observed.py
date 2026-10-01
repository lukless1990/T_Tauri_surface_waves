#!/usr/bin/env python3
"""compare_observed.py — compare the simulated Moreton wave kinematics with the
observed 2006 December 6 event (Balasubramaniam et al. 2010, ApJ 723, 587 —
PDF in RefMoreton/), the same event simulated by Krause et al. 2015.

Usage: python3 compare_observed.py run_<label>
Writes plots_<label>/fig5_obs_comparison.png and prints the comparison table.

Method (mirrors the paper): track the front's leading edge, fit the constant-
deceleration model d(t) = d0 + v0*t + a*t^2/2 (their Eq. 1), and read off the
mean speed over the tracked interval. Front tracking uses the coronal dp/p
trace (z=15 Mm; the coronal shock IS the wave; the chromospheric signature is
its skirt) with a fixed detection threshold, plus the chromospheric vz trace
for the visibility range.
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rundir = (sys.argv[1] if len(sys.argv) > 1 else "run_production").rstrip("/")
label  = os.path.basename(rundir).replace("run_", "")
outdir = os.path.join(os.path.dirname(rundir) or ".", f"plots_{label}")
os.makedirs(outdir, exist_ok=True)

meta = {}
with open(os.path.join(rundir, "meta.txt")) as f:
    for ln in f:
        k, v = ln.strip().split("=")
        try: meta[k] = float(v)
        except ValueError: meta[k] = v
x_fl = meta["x_fl_Mm"]

trt  = np.loadtxt(os.path.join(rundir, "trace_times.txt"))
trx  = np.loadtxt(os.path.join(rundir, "trace_x_Mm.txt"))
trvz = np.loadtxt(os.path.join(rundir, "trace_chromo_vz.txt"))
trdp = np.loadtxt(os.path.join(rundir, "trace_corona_dpp.txt"))

# ---- front tracking: rightmost crossing of a fixed threshold ----------------
def front(stack, thresh):
    pos = np.full(stack.shape[0], np.nan)
    for k in range(stack.shape[0]):
        idx = np.where((np.abs(stack[k]) > thresh) & (trx > x_fl + 10))[0]
        if len(idx): pos[k] = trx[idx[-1]]
    return pos

pos_cor = front(trdp, 0.03)          # coronal front: dp/p > 3%
pos_ha  = front(trvz, 0.30)          # chromospheric swing: |vz| > 0.3 km/s

# keep only the monotonically advancing part of the coronal track: once the
# front's dp/p decays below threshold the tracker falls back to trailing
# structure (that decay distance is itself a result, reported below)
k_end = int(np.nanargmax(pos_cor))
d_cor_max = pos_cor[k_end] - x_fl
mono = np.zeros_like(pos_cor, bool); best = -np.inf
for k in range(k_end+1):
    if np.isfinite(pos_cor[k]) and pos_cor[k] > best:
        mono[k] = True; best = pos_cor[k]
msk = mono & (pos_cor > x_fl + 30) & (pos_cor < trx[-1] - 20)
tfit, dfit = trt[msk], (pos_cor[msk] - x_fl)          # distance from source [Mm]

# constant-deceleration fit d(t) = d0 + v0 t + 0.5 a t^2   (paper Eq. 1)
c2, c1, c0 = np.polyfit(tfit, dfit, 2)
a_fit  = 2*c2*1e3                    # Mm/s^2 -> km/s^2
v0_fit = c1*1e3                      # km/s  (at t=0 of the fit's time axis)
v_at   = lambda t: (c1 + 2*c2*t)*1e3
v_start, v_end = v_at(tfit[0]), v_at(tfit[-1])
v_mean = (dfit[-1]-dfit[0])/(tfit[-1]-tfit[0])*1e3

# chromospheric visibility range (distance from source where the swing exceeds
# the threshold, analog of "first seen / lost" in the observations)
vis = np.isfinite(pos_ha) & (pos_ha > x_fl + 30)      # outside the flare guard zone
d_first = np.nanmin(pos_ha[vis]) - x_fl if vis.any() else np.nan
d_last  = np.nanmax(pos_ha[vis]) - x_fl if vis.any() else np.nan
at_edge = vis.any() and np.nanmax(pos_ha) > trx[-1] - 30
# far-field chromospheric swing amplitude (x-x_fl > 50 Mm, clear of the source)
far = trx > x_fl + 50
swing = np.nanmax(np.abs(trvz[:, far]))

# ---- observed values: Balasubramaniam et al. 2010 (2006 Dec 6, X6.5/3B) ----
OBS = dict(v_mean=850., v0=1125., a=-0.87, d_first=120., d_last=550.)

# ------------------------------------------------------------------ figure
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
ax[0].plot(trt, pos_cor - x_fl, "C0.", ms=3, label="simulated coronal front (dp/p > 3%)")
ax[0].plot(trt, pos_ha - x_fl, "C3.", ms=3, label="simulated chromospheric swing (|v$_z$| > 0.3 km/s)")
tt = np.linspace(tfit[0], tfit[-1], 100)
ax[0].plot(tt, c0 + c1*tt + c2*tt**2, "k-", lw=1.2,
           label=f"decel. fit: v$_0$={v_start:.0f} km/s, a={a_fit:.2f} km/s$^2$")
# observed reference curve, shifted to launch at our fit's origin
t_obs = np.linspace(0, 480, 100)
d_obs = OBS["v0"]*1e-3*t_obs + 0.5*OBS["a"]*1e-3*t_obs**2
t0 = tfit[0] - (dfit[0] - 0)/ (OBS["v0"]*1e-3)   # rough alignment at first detection
ax[0].plot(t_obs + max(t0, 0), d_obs, "--", color="gray", lw=1.5,
           label=f"observed 2006 Dec 6 fit (v$_0$={OBS['v0']:.0f}, a={OBS['a']} km/s$^2$)")
ax[0].set_xlabel("t [s]"); ax[0].set_ylabel("distance from source [Mm]")
ax[0].set_ylim(0, max(np.nanmax(pos_cor - x_fl)*1.15, 260)); ax[0].set_xlim(0, trt[-1])
ax[0].legend(fontsize=7.5, loc="upper left"); ax[0].grid(alpha=0.3)
ax[0].set_title("front kinematics: simulation vs 2006 Dec 6")

rows = [
    ("mean front speed [km/s]",        f"{v_mean:.0f}",               f"{OBS['v_mean']:.0f}"),
    ("initial speed v0 [km/s]",        f"{v_start:.0f}",              f"{OBS['v0']:.0f}"),
    ("deceleration [km/s$^2$]",        f"{a_fit:.2f}",                f"{OBS['a']:.2f}"),
    ("final tracked speed [km/s]",     f"{v_end:.0f}",                "~700-800"),
    ("chromo. signal first seen [Mm]", f"{d_first:.0f}",              f"~{OBS['d_first']:.0f}"),
    ("chromo. signal tracked to [Mm]", (f">{d_last:.0f} (domain edge)" if at_edge else f"{d_last:.0f}"),
                                       f"~{OBS['d_last']:.0f}"),
    ("far-field chromo. swing [km/s]", f"{swing:.1f} peak",           "1-4 (typical events)"),
    ("coronal signal decayed at [Mm]", f"{d_cor_max:.0f} (dp/p < 3%)",  "EUV front fades ~R_sun scale"),
]
ax[1].axis("off")
tbl = ax[1].table(cellText=[[r[1], r[2]] for r in rows],
                  rowLabels=[r[0] for r in rows],
                  colLabels=[f"simulation ({label})", "observed"],
                  loc="center", cellLoc="center")
tbl.auto_set_font_size(False); tbl.set_fontsize(8.5); tbl.scale(0.72, 1.55)
ax[1].set_title("Balasubramaniam et al. 2010 (ApJ 723, 587)", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(outdir, "fig5_obs_comparison.png"), dpi=150); plt.close(fig)

print(f"fit interval: t = {tfit[0]:.0f}-{tfit[-1]:.0f} s, d = {dfit[0]:.0f}-{dfit[-1]:.0f} Mm from source")
for r in rows: print(f"  {r[0]:34s} sim: {r[1]:22s} obs: {r[2]}")
print(f"figure -> {outdir}/fig5_obs_comparison.png")
