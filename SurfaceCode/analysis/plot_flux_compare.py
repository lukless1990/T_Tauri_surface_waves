#!/usr/bin/env python3
# plot_flux_compare.py — overlay the energy-flux decay of the two single-pulse runs (1 kG vs 0.1 kG
# stellar field), both launched at the same 10 km/s. A weaker field decays the fast-mode wave more
# slowly (less magnetic steepening -> more acoustic), so it carries further before being buried in the
# ambient convective-wave background -- but still nowhere near the pole.
#
# Metric: the single-exponential fit is fragile for the weak field (its decay is two-slope: a shallow
# shoulder while the smooth pulse steepens, then a steep drop once the shock forms). We therefore quote
# the definition-free MEASURED ambient-crossing x_amb (where Phi(x) drops below the convective floor),
# and a robust effective decay length L_eff = x_amb / ln(Phi0/Phi_amb) averaged over that drop.
import numpy as np, matplotlib, os
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(__file__)
outd = os.path.join(here, "plots"); os.makedirs(outd, exist_ok=True)
pole = 0.61
runs = [("paper_run_1kG",   "tab:blue", "1 kG"),
        ("paper_run_0p1kG", "tab:red",  "0.1 kG")]
v_amb = 0.28

def load(lbl):
    dd = os.path.join(here, "..", "..", "output", lbl)
    x, phi = np.loadtxt(os.path.join(dd, "flux_final.txt")).T
    pk = np.atleast_2d(np.loadtxt(os.path.join(dd, "snap_times.txt")))[:, 1]
    meta = dict(l.split("=") for l in open(os.path.join(dd, "meta.txt")).read().split() if "=" in l)
    return x, np.abs(phi)/np.abs(phi).max(), pk[0], meta

def ambient_crossing(x, phin, floor):
    # first x>0.02 where the (monotone-ish) curve drops below the floor; linear-in-log interpolation
    below = np.where((x > 0.02) & (phin < floor))[0]
    if not len(below): return np.nan
    i = below[0]; x0, x1 = x[i-1], x[i]; y0, y1 = np.log(phin[i-1]), np.log(phin[i])
    return x0 + (np.log(floor) - y0) * (x1 - x0) / (y1 - y0)

fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
floor = (v_amb / 10.0)**2
summ = {}
for lbl, col, name in runs:
    x, phin, vlaunch, meta = load(lbl)
    bkg = float(meta.get("Bstar_kG", "nan"))
    xamb = ambient_crossing(x, phin, floor)
    Leff = xamb / np.log(1.0/floor)
    summ[name] = (xamb, Leff)
    ax.semilogy(x, phin, '-', color=col, lw=1.6, alpha=0.9,
                label=fr"$B_\star = {bkg:.1f}$ kG")
    ax.plot([xamb], [floor], 'o', color=col, ms=8, zorder=5, mec="k", mew=0.6)

ax.axhline(floor, color="tab:green", lw=1.5)
ax.axhspan(1e-8, floor, color="tab:green", alpha=0.10)
ax.text(0.004, floor*1.8, r"$v_\perp\!\sim\!0.28$ km/s", fontsize=12, color="darkgreen")
# clip just below the ambient floor: below it the wave is buried (undetectable vs the convective
# background), so the sub-floor pulse wake is not physically meaningful to show.
ax.set_xlim(0, 0.1); ax.set_ylim(2e-4, 3)
ax.set_xlabel(r"lateral distance from footpoint  $x/R_\star$", fontsize=13)
ax.set_ylabel(r"energy-flux $\Phi(x)/\Phi_0$", fontsize=13)
ax.tick_params(which='both', direction='in', top=True, right=True, labelsize=12)
ax.legend(loc="upper right", fontsize=12)
plt.savefig(os.path.join(outd, "flux_decay_field_compare.png"), dpi=140)
for name,(xa,le) in summ.items():
    print(f"{name}: buried at x={xa:.4f} R*, L_eff={le:.4f} R*, pole {pole/xa:.1f}x further")
print("wrote flux_decay_field_compare.png")
