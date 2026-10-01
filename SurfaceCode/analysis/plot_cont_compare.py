#!/usr/bin/env python3
# plot_cont_compare.py — continuous-driving vs single-pulse decay comparison.
#   Overlays the steady-state time-averaged energy flux <F(x)> (continuous driving, output/paper_run_cont/)
#   on the single-pulse Phi(x) (output/paper_run/), both normalised. The point: continuous driving reaches
#   a steady state whose flux decays with the SAME length L — so sustained accretion driving does not
#   deliver the wave any further toward the pole.
import numpy as np, matplotlib, os, glob
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm

here = os.path.dirname(__file__)
dc   = os.path.join(here, "..", "..", "output", "paper_run_cont")
dp   = os.path.join(here, "..", "..", "output", "paper_run_1kG")   # single-pulse run (1 kG, same field as the cont run)
outd = os.path.join(here, "plots"); os.makedirs(outd, exist_ok=True)
pole = 0.61

def fit_L(x, f, xhi=0.065):
    fn = np.abs(f)/np.abs(f).max()
    m = (x > 0.004) & (x < xhi) & (fn > 1e-6)
    sl, b = np.polyfit(x[m], np.log(fn[m]), 1)
    return fn, -1/sl, sl, b

xc, Fc = np.loadtxt(os.path.join(dc, "flux_steady.txt")).T   # continuous: <F(x)>
xp, Fp = np.loadtxt(os.path.join(dp, "flux_final.txt")).T    # pulse: Phi(x)
Fcn, Lc, slc, bc = fit_L(xc, Fc)
Fpn, Lp, slp, bp = fit_L(xp, Fp)

# ambient floor from the launched amplitude (same 10 km/s calibration for both)
pk = np.atleast_2d(np.loadtxt(os.path.join(dc, "snap_times.txt")))[:, 1]
v_amb = 0.28; v_inj = pk.max()
floor = (v_amb / v_inj)**2
xcross_c = -np.log(floor) * Lc

fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
ax.semilogy(xp, Fpn, '-', color="tab:blue", lw=1.4, alpha=0.9, label=r"single pulse $\Phi(x)/\Phi_0$")
ax.semilogy(xc, Fcn, '-', color="tab:orange", lw=1.4, alpha=0.9,
            label=r"continuous drive $\langle F(x)\rangle/\langle F\rangle_0$ (steady state)")
xe = np.linspace(0, 0.09, 200)
ax.semilogy(xe, np.exp(bp + slp*xe), '--', color="tab:blue",   lw=1.2)
ax.semilogy(xe, np.exp(bc + slc*xe), '--', color="tab:orange", lw=1.2)
ax.text(0.05, 0.5, fr"$L_{{\rm pulse}}={Lp:.4f}\,R_\star$""\n"fr"$L_{{\rm cont}}={Lc:.4f}\,R_\star$",
        fontsize=10, va="top")
ax.axhline(floor, color="tab:green", lw=1.5)
ax.axhspan(1e-8, floor, color="tab:green", alpha=0.10)
ax.text(pole*0.4, floor*1.7, r"ambient convective-wave background ($v_\perp\!\sim\!0.28$ km/s)",
        fontsize=9, color="darkgreen")
ax.axvline(pole, color="k", ls=":", lw=1.3)
ax.text(pole-0.008, 1e-6, f"pole\n(0.61 $R_\\star$,\n{pole/xcross_c:.0f}× further)", ha="right", fontsize=9)
ax.axvspan(0, xc[-1], color="0.8", alpha=0.20)
ax.text(xc[-1]*0.5, 1.6, "simulated patch (0.1 $R_\\star$)", ha="center", fontsize=8, color="0.4")
ax.set_xlim(0, pole*1.02); ax.set_ylim(1e-7, 3)
ax.set_xlabel(r"lateral distance from footpoint  $x/R_\star$")
ax.set_ylabel(r"normalised energy flux")
ax.set_title("Continuous accretion driving decays with the same length as a single pulse\n"
             f"(res×2, 0.1 $R_\\star$; $L_{{\\rm cont}}$≈{Lc:.4f} vs $L_{{\\rm pulse}}$≈{Lp:.4f} $R_\\star$ — "
             "sustained driving does not reach the pole)")
ax.legend(loc="upper right", fontsize=9)
plt.savefig(os.path.join(outd, "cont_vs_pulse_decay.png"), dpi=140)
print(f"wrote cont_vs_pulse_decay.png (L_cont={Lc:.4f}, L_pulse={Lp:.4f} R*, xcross_cont={xcross_c:.3f})")

# ---- 3-panel v_x(x,z) for the continuous run: the steady wave train filling the domain,
#      amplitude decaying laterally. z/Hp in [0,2]; uniform symlog colour across panels. ----
x  = np.loadtxt(os.path.join(dc, "grid_x.txt"))
z  = np.loadtxt(os.path.join(dc, "grid_z.txt"))
st = np.atleast_2d(np.loadtxt(os.path.join(dc, "snap_times.txt")))
tt = st[:, 0]
snaps = sorted(glob.glob(os.path.join(dc, "snap_vx_*.txt")))
def fmt_t(s): return f"{s:.0f} s" if s<600 else (f"{s/60:.0f} min" if s<7200 else f"{s/3600:.2f} h")
targets = [225.0, 5400.0, 14400.0]                      # initial (1.5τ), 90 min, 4 h — match the pulse panels
sel = [int(np.argmin(np.abs(tt - tg))) for tg in targets]
zmask = (z >= 0.0) & (z <= 2.0); zc = z[zmask]
fields = [np.loadtxt(snaps[k]).T[zmask, :] for k in sel]
vmax = max(np.abs(f).max() for f in fields)
lin  = 0.03
norm = SymLogNorm(linthresh=lin, linscale=0.8, vmin=-vmax, vmax=vmax, base=10)
fig, axes = plt.subplots(3, 1, figsize=(9, 5.4), constrained_layout=True, sharex=True)
ext = [x[0], x[-1], zc[0], zc[-1]]
dxx = x[1] - x[0]; nw = max(1, int(round(0.002/dxx)))   # ~1 wavelength window behind the front
thr = 0.02                                              # front-detection threshold [km/s] (>> numerical noise)
for ax, k, f in zip(axes, sel, fields):
    im = ax.imshow(f, origin="lower", aspect="auto", extent=ext, cmap="RdBu_r", norm=norm)
    env = np.abs(f).max(axis=0)                         # z-max |vx| envelope vs x
    idx = np.where(env > thr)[0]
    if len(idx):
        i_f = idx[-1]; x_f = x[i_f]
        front_v = env[max(0, i_f-nw):i_f+1].max()       # leading-crest amplitude at the front
        ann = f"wave front: |$v_x$| = {front_v:.2f} km/s\nat x = {x_f:.3f} $R_\\star$"
        ax.axvline(x_f, color="0.25", ls="--", lw=0.9)
    else:
        ann = "wave front below threshold"
    ax.text(0.012, 0.80, fmt_t(tt[k]), transform=ax.transAxes, fontsize=9,
            bbox=dict(fc="white", alpha=0.8, ec="none"))
    ax.text(0.988, 0.80, ann, transform=ax.transAxes, fontsize=8.5, ha="right",
            color="firebrick", bbox=dict(fc="white", alpha=0.8, ec="none"))
    ax.set_ylabel(r"$z/H_p$"); ax.set_ylim(0, 2)
axes[-1].set_xlabel(r"lateral distance from footpoint  $x/R_\star$")
fig.colorbar(im, ax=axes, shrink=0.9, pad=0.015, label=r"$v_x$ [km/s]")
plt.savefig(os.path.join(outd, "cont_pulse_panels.png"), dpi=140)
print(f"wrote cont_pulse_panels.png (3 panels @ t={[round(tt[k]) for k in sel]} s; vmax={vmax:.2f})")
