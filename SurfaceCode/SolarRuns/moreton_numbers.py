#!/usr/bin/env python3
"""moreton_numbers.py — Appendix B numbers: model front vs Francile+2013 Eq. (1).

Fills the placeholders: % agreement over the observed coverage, the fitted power-law
exponent delta vs their 0.578627, the mean front speed vs the observed 838 km/s, and
the trackable range. Uses the 1 s trace stacks when the run has finished, else the
50 s snapshots (front positions are identical, the cadence is not).

Usage: python3 moreton_numbers.py run_kr270full [TMAX]
"""
import os, sys, glob, numpy as np
from scipy.optimize import curve_fit
import importlib.util
spec = importlib.util.spec_from_file_location("mmf", "make_moreton_fig.py")

C1, DELTA, C2 = 15.287, 0.578627, -108.609
T_OBS = (86.0, 545.0)
francile = lambda t: C1*np.asarray(t, float)**DELTA + C2

run  = sys.argv[1] if len(sys.argv) > 1 else "run_kr270full"
TMAX = float(sys.argv[2]) if len(sys.argv) > 2 else np.inf

# reuse the figure script's tracker (identical thresholds => identical numbers)
src = open("make_moreton_fig.py").read().split("run = sys.argv")[0]
ns = {"__file__": os.path.abspath("make_moreton_fig.py")}; exec(src, ns)
t, d = ns["track"](run)
m = t <= TMAX; t, d = t[m], d[m]
src_kind = "1 s traces" if os.path.exists(f"{run}/trace_times.txt") else "50 s snapshots"

print(f"=== {run}  [{src_kind}, {len(t)} points, t = {t[0]:.0f}-{t[-1]:.0f} s] ===\n")
print(f"{'t [s]':>7} {'model':>8} {'Eq.(1)':>8} {'diff':>8}")
show = t if len(t) < 30 else t[np.minimum(np.searchsorted(t, np.arange(50, t[-1], 50)), len(t)-1)]
for tk in show:
    k = int(np.argmin(abs(t-tk))); f = francile(t[k])
    print(f"{t[k]:7.0f} {d[k]:8.1f} {f:8.1f} {100*(d[k]-f)/f:+7.1f}%")

sel = (t >= T_OBS[0]) & (t <= min(T_OBS[1], TMAX))
res = 100*(d[sel]-francile(t[sel]))/francile(t[sel])
print(f"\nover observed coverage {T_OBS[0]:.0f}-{min(T_OBS[1],TMAX):.0f} s:")
print(f"  agreement      {res.mean():+.1f}% mean, {res.min():+.1f} to {res.max():+.1f}% range, "
      f"{np.abs(res).max():.1f}% max |diff|")

# power-law fit, same functional form as their Eq. (1)
try:
    p, _ = curve_fit(lambda x,c1,de,c2: c1*x**de+c2, t[sel], d[sel], p0=[C1, DELTA, C2], maxfev=20000)
    print(f"  fitted         d = {p[0]:.3f} t^{p[1]:.4f} {p[2]:+.2f}   (delta {p[1]:.3f} vs their {DELTA:.3f})")
except Exception as e:
    print("  fit failed:", e)

ta, tb = T_OBS[0], min(T_OBS[1], t[-1])
va = (np.interp(tb, t, d)-np.interp(ta, t, d))/(tb-ta)*1e3
vo = (francile(tb)-francile(ta))/(tb-ta)*1e3
print(f"  mean speed     {va:.0f} km/s   (observed {vo:.0f} km/s, {100*(va-vo)/vo:+.1f}%)")
print(f"  trackable range {d[-1]:.0f} Mm at t = {t[-1]:.0f} s   (observed {francile(545):.0f} Mm at 545 s)")
