#!/usr/bin/env python3
# plot_tau_compare.py — driver-duration test: does the decay length scale with the pulse wavelength?
#
# The paper's production runs drive the footpoint for tau = 150 s. The clump passage time implied by
# Cranmer (2008)'s clump scale (r_c ~ 1e-2 R*, v_ff ~ 300 km/s) is tau ~ 50 s, i.e. a 3x SHORTER pulse
# and hence a 3x shorter wavelength.
#
# RESULT (2026-07-25): L is UNCHANGED (0.0062 vs 0.0060 R*, within the fit-window scatter). The naive
# expectation L ~ lambda from weak-shock theory does NOT hold here, because both pulses steepen within
# x_sh ~ 1e-4 R* -- far inside the first fitted point -- and the resulting N-wave broadens as it decays,
# erasing the launch duration. L is a property of the medium + amplitude, not of the driver. Same
# conclusion the continuously driven run reaches from the opposite extreme of driving history.
#
# Both runs are launched at the same peak amplitude (v_inj ~ 10 km/s, M0 recalibrated per tau) on the
# same domain, so the ONLY difference is the driver duration. Fit window identical to plot_paper_run.py.
import numpy as np, matplotlib, os, sys
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
outd = os.path.join(here, "plots"); os.makedirs(outd, exist_ok=True)
pole = 0.61
v_amb = 0.28
runs = [("paper_run_1kG",   "tab:blue", r"$\tau = 150$ s (production)"),
        ("paper_run_tau50", "tab:orange", r"$\tau = 50$ s (clump-consistent)")]

def load(lbl):
    d = os.path.join(here, "..", "output", lbl)
    xf, phi = np.loadtxt(os.path.join(d, "flux_final.txt")).T
    pk = np.atleast_2d(np.loadtxt(os.path.join(d, "snap_times.txt")))[:, 1]
    meta = dict(l.split("=") for l in open(os.path.join(d, "meta.txt")).read().split() if "=" in l)
    return xf, np.abs(phi) / np.abs(phi).max(), pk[0], meta

def fit_L(xf, phin):
    # identical window to plot_paper_run.py: skip the near field (x<0.004) and the arriving front (x>0.065)
    m = (xf > 0.004) & (xf < 0.065) & (phin > 1e-6)
    sl, b = np.polyfit(xf[m], np.log(phin[m]), 1)
    return -1.0 / sl, b, sl

def ambient_crossing(xf, phin, floor):
    below = np.where((xf > 0.02) & (phin < floor))[0]
    if not len(below): return np.nan
    i = below[0]; x0, x1 = xf[i-1], xf[i]; y0, y1 = np.log(phin[i-1]), np.log(phin[i])
    return x0 + (np.log(floor) - y0) * (x1 - x0) / (y1 - y0)

fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
summ = {}
for lbl, col, name in runs:
    try:
        xf, phin, v_inj, meta = load(lbl)
    except OSError as e:
        print(f"!! missing run {lbl}: {e}"); sys.exit(1)
    tau = float(meta.get("tau_s", 150.0))
    L, b, sl = fit_L(xf, phin)
    floor = (v_amb / v_inj) ** 2
    xamb = ambient_crossing(xf, phin, floor)
    summ[name] = dict(tau=tau, L=L, xamb=xamb, v_inj=v_inj, M0=float(meta["M0"]), floor=floor,
                      lam=2 * tau * float(meta["cf_kms"]))
    ax.semilogy(xf, phin, '-', color=col, lw=1.6, alpha=0.9, label=fr"{name}: $L={L:.4f}\,R_\star$")
    xe = np.linspace(0.004, 0.09, 100)
    ax.semilogy(xe, np.exp(b + sl * xe), '--', color=col, lw=1.1, alpha=0.7)
    if np.isfinite(xamb):
        ax.plot([xamb], [floor], 'o', color=col, ms=8, zorder=5, mec="k", mew=0.6)

fl = np.mean([s["floor"] for s in summ.values()])
ax.axhline(fl, color="tab:green", lw=1.5)
ax.axhspan(1e-9, fl, color="tab:green", alpha=0.10)
ax.text(0.004, fl * 1.8, r"$v_\perp\!\sim\!0.28$ km/s (ambient)", fontsize=11, color="darkgreen")
ax.set_xlim(0, 0.1); ax.set_ylim(1e-6, 3)
ax.set_xlabel(r"lateral distance from footpoint  $x/R_\star$", fontsize=13)
ax.set_ylabel(r"energy-flux $\Phi(x)/\Phi_0$", fontsize=13)
ax.tick_params(which='both', direction='in', top=True, right=True, labelsize=12)
ax.legend(loc="upper right", fontsize=11)
ax.set_title(r"Decay length is independent of driver duration ($\lambda$ differs by 3$\times$)", fontsize=11)
plt.savefig(os.path.join(outd, "flux_decay_tau_compare.png"), dpi=140)

print(f"{'run':32s} {'tau[s]':>7s} {'lam[km]':>9s} {'M0':>6s} {'v_inj':>7s} {'L/R*':>8s} {'x_amb/R*':>9s}")
for name, s in summ.items():
    print(f"{name:32s} {s['tau']:7.0f} {s['lam']:9.0f} {s['M0']:6.2f} {s['v_inj']:7.2f} "
          f"{s['L']:8.4f} {s['xamb']:9.4f}")
ks = list(summ)
if len(ks) == 2:
    a, c = summ[ks[0]], summ[ks[1]]
    print(f"\nratio tau: {a['tau']/c['tau']:.2f}   ratio L: {a['L']/c['L']:.2f}   "
          f"ratio x_amb: {a['xamb']/c['xamb']:.2f}")
    print(f"pole / L(tau=50) = {pole/c['L']:.0f} e-foldings;  pole / x_amb(tau=50) = {pole/c['xamb']:.0f}x")
print("wrote plots/flux_decay_tau_compare.png")
