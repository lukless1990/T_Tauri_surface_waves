#!/usr/bin/env python3
"""tradeoff.py — range vs chromospheric-swing trade-off across the driver scan.

For each run: driver amplitude A, mean front speed, range the coronal front stays
above dp/p=3% (proxy for observable extent), and the far-field (x-x_fl>50 Mm)
chromospheric vertical swing amplitude. The observed 2006 Dec 6 target is
~550 Mm range AND ~2-4 km/s swing; the scan shows whether a single sustained
driver can hit both.

Usage: python3 tradeoff.py run_mid20 run_mid50 run_extreme [...]
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

runs = sys.argv[1:] or ["run_mid20", "run_mid50", "run_extreme"]

def summarize(rundir):
    rundir = rundir.rstrip("/")
    meta = {}
    with open(os.path.join(rundir, "meta.txt")) as f:
        for ln in f:
            k, v = ln.strip().split("=")
            try: meta[k] = float(v)
            except ValueError: meta[k] = v
    x_fl = meta["x_fl_Mm"]; A = meta["A_pulse"]
    lab = os.path.basename(rundir).replace("run_", "")
    # meta only recorded 'drive' from this run on; the mid*/extreme scan runs are piston
    drive = meta.get("drive", "piston" if (lab.startswith("mid") or lab == "extreme") else "pulse")
    t  = np.loadtxt(os.path.join(rundir, "trace_times.txt"))
    x  = np.loadtxt(os.path.join(rundir, "trace_x_Mm.txt"))
    vz = np.loadtxt(os.path.join(rundir, "trace_chromo_vz.txt"))
    dp = np.loadtxt(os.path.join(rundir, "trace_corona_dpp.txt"))
    # coronal front: rightmost dp/p>3%, then max advancing distance (the range)
    pos = np.full(dp.shape[0], np.nan)
    for k in range(dp.shape[0]):
        idx = np.where((np.abs(dp[k]) > 0.03) & (x > x_fl + 30))[0]
        if len(idx): pos[k] = x[idx[-1]]
    rng = (np.nanmax(pos) - x_fl) if np.any(np.isfinite(pos)) else np.nan
    # mean speed over the monotone advancing segment
    k_end = int(np.nanargmax(pos)) if np.any(np.isfinite(pos)) else 0
    mono = []
    best = -np.inf
    for k in range(k_end + 1):
        if np.isfinite(pos[k]) and pos[k] > best:
            mono.append((t[k], pos[k] - x_fl)); best = pos[k]
    mono = np.array(mono)
    m = mono[(mono[:,1] > 30)] if len(mono) else np.empty((0,2))
    vmean = (m[-1,1]-m[0,1])/(m[-1,0]-m[0,0])*1e3 if len(m) > 3 else np.nan
    # far-field chromospheric swing
    far = x > x_fl + 50
    swing = np.nanmax(np.abs(vz[:, far])) if far.any() else np.nan
    return dict(label=os.path.basename(rundir).replace("run_",""), A=A, drive=drive,
                vmean=vmean, rng=rng, swing=swing)

rows = []
for r in runs:
    try: rows.append(summarize(r))
    except FileNotFoundError: pass
# include the impulsive production point if present
if os.path.isdir("run_production"):
    try: rows.append(summarize("run_production"))
    except Exception: pass
rows.sort(key=lambda d: d["A"])

print(f"{'run':12s} {'drive':7s} {'A':>5s} {'mean[km/s]':>11s} {'range[Mm]':>10s} {'swing[km/s]':>12s}")
for d in rows:
    print(f"{d['label']:12s} {d['drive']:7s} {d['A']:5.0f} {d['vmean']:11.0f} {d['rng']:10.0f} {d['swing']:12.1f}")

# ---- trade-off figure -------------------------------------------------------
pis = [d for d in rows if d["drive"] == "piston"]
fig, ax = plt.subplots(figsize=(7.2, 5.2))
if pis:
    A   = np.array([d["A"] for d in pis])
    rng = np.array([d["rng"] for d in pis])
    sw  = np.array([d["swing"] for d in pis])
    ax.plot(rng, sw, "o-", color="C1", ms=8, label="sustained CME-piston scan")
    for d in pis:
        ax.annotate(f"A={d['A']:.0f}", (d["rng"], d["swing"]),
                    textcoords="offset points", xytext=(8, 4), fontsize=9)
# observed target box
ax.axvspan(500, 600, color="gray", alpha=0.12)
ax.axhspan(1, 4, color="C2", alpha=0.15)
ax.axvline(550, color="gray", ls=":", lw=1);
ax.text(552, ax.get_ylim()[1]*0.5 if pis else 50, "observed range ~550 Mm",
        rotation=90, va="center", fontsize=8, color="gray")
ax.axhspan(1, 4, color="C2", alpha=0.0)
ax.set_yscale("log")
ax.set_xlabel("observable range [Mm]  (coronal front dp/p > 3%)")
ax.set_ylabel("far-field chromospheric swing [km/s]")
ax.set_title("range vs chromospheric swing: the sustained-driver trade-off")
ax.axhline(4, color="C2", lw=1); ax.axhline(1, color="C2", lw=1)
ax.text(60, 2, "observed swing\n1-4 km/s", color="C2", fontsize=8, va="center")
ax.grid(alpha=0.3, which="both"); ax.legend(loc="upper left")
fig.tight_layout()
os.makedirs("plots_tradeoff", exist_ok=True)
fig.savefig("plots_tradeoff/fig7_tradeoff.png", dpi=150); plt.close(fig)
print("figure -> plots_tradeoff/fig7_tradeoff.png")
