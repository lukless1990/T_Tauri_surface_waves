#!/usr/bin/env python3
# plot_paper_run.py — paper figures from output/paper_run/ (the res×2, 0.1 R*, 10 km/s run).
#   Fig A (pulse_decay_panels.png): v_x(x,z) at a selection of 30-min snapshots — the pulse skimming
#     the surface and collapsing. Per-panel auto-scaled (each shows its own structure); the decay is
#     read from the annotated peak |v_x| and absolute time on each panel.
#   Fig B (flux_decay_0p1.png): height/time-integrated energy flux Phi(x)/Phi0 out to 0.1 R*, with an
#     exponential fit (decay length L), the ambient convective-wave floor, and the pole marker.
#
# Usage: python SurfaceCode/analysis/plot_paper_run.py [n_panels]   (default 6; use "all" for every snapshot)
import numpy as np, matplotlib, os, glob, sys
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

label = os.environ.get("RUN_LABEL", "paper_run")            # data subdir (RUN_LABEL=paper_run_1kG etc.)
d   = os.path.join(os.path.dirname(__file__),"..","..","output", label)
outd= os.path.join(os.path.dirname(__file__), "plots"); os.makedirs(outd, exist_ok=True)
Rstar_pole = 0.61                                  # ring->pole transit in R*

meta = dict(l.split("=") for l in open(os.path.join(d,"meta.txt")).read().split() if "=" in l)
bkg  = float(meta.get("Bstar_kG", "1.0")); tag = f"{bkg:g}kG".replace(".", "p")   # 1.0->1kG, 0.1->0p1kG
B0G  = float(meta.get("B0_G", "nan"))
x  = np.loadtxt(os.path.join(d,"grid_x.txt"))      # x/R*  (Nx)
z  = np.loadtxt(os.path.join(d,"grid_z.txt"))      # z/Hp  (Nz)
st = np.atleast_2d(np.loadtxt(os.path.join(d,"snap_times.txt")))   # cols: t[s], domain peak |vx| [km/s]
tt, pk = st[:,0], st[:,1]
snaps = sorted(glob.glob(os.path.join(d,"snap_vx_*.txt")))
nsnap = len(snaps)

def fmt_t(s):  # absolute time label
    return f"{s:.0f} s" if s<600 else (f"{s/60:.0f} min" if s<7200 else f"{s/3600:.2f} h")

# ---- Fig A: three v_x(x,z) panels (launch, 90 min, 4 h); z/Hp in [0,2];
#      UNIFORM symlog colour scale across all panels (so the decay is visible directly) ----
from matplotlib.colors import SymLogNorm
targets = [225.0, 5400.0, 14400.0]                      # launch, 90 min, 4 h  [s]
sel = [int(np.argmin(np.abs(tt - tg))) for tg in targets]
zmask = (z >= 0.0) & (z <= 2.0); zc = z[zmask]
fields = [np.loadtxt(snaps[k]).T[zmask, :] for k in sel]  # each (nz_crop, Nx), km/s
vmax = max(np.abs(f).max() for f in fields)             # shared across panels
lin  = 0.03                                             # symlog linear threshold [km/s]
norm = SymLogNorm(linthresh=lin, linscale=0.8, vmin=-vmax, vmax=vmax, base=10)
fig, axes = plt.subplots(3, 1, figsize=(9, 5.4), constrained_layout=True, sharex=True)
ext = [x[0], x[-1], zc[0], zc[-1]]
dxx = x[1] - x[0]; nw = max(1, int(round(0.002/dxx)))    # ~1 wavelength window behind the front
thr = 0.006                                              # front-detection threshold [km/s]: low enough to
                                                         # catch the faint oscillatory leading crests (which
                                                         # sit below 0.02 but are still visible) so the marker
                                                         # reaches the true visible front
for ax, k, f in zip(axes, sel, fields):
    im = ax.imshow(f, origin="lower", aspect="auto", extent=ext, cmap="RdBu_r", norm=norm)
    env = np.abs(f).max(axis=0)                          # z-max |vx| envelope vs x
    idx = np.where(env > thr)[0]
    if len(idx):
        i_f = idx[-1]; x_f = x[i_f]
        front_v = env[max(0, i_f-nw):i_f+1].max()        # leading-crest |vx| at the front
        ann = f"wave front: |$v_x$| = {front_v:.2f} km/s\nat x = {x_f:.3f} $R_\\star$"
        ax.axvline(x_f, color="0.25", ls="--", lw=0.9)
    else:
        ann = "wave front below threshold"
    lbl = "initial pulse" if tt[k] < 400 else fmt_t(tt[k])
    ax.text(0.012, 0.80, lbl, transform=ax.transAxes, fontsize=10,
            bbox=dict(fc="white", alpha=0.8, ec="none"))
    ax.text(0.988, 0.80, ann, transform=ax.transAxes, fontsize=9, ha="right",
            color="firebrick", bbox=dict(fc="white", alpha=0.8, ec="none"))
    ax.set_ylabel(r"$z/H_p$", fontsize=12); ax.set_ylim(0, 2)
    ax.tick_params(which='both', direction='in', top=True, right=True, labelsize=11)
axes[-1].set_xlabel(r"lateral distance from footpoint  $x/R_\star$", fontsize=12)
cb = fig.colorbar(im, ax=axes, shrink=0.9, pad=0.015, label=r"$v_x$ [km/s]")
plt.savefig(os.path.join(outd,f"pulse_decay_panels_{tag}.png"), dpi=140)
print(f"wrote pulse_decay_panels_{tag}.png (3 panels @ t={[round(tt[k]) for k in sel]} s; vmax={vmax:.2f})")

# ---- Fig B: Phi(x) decay to the ambient floor, out to 0.1 R* ----
xf, phi = np.loadtxt(os.path.join(d,"flux_final.txt")).T
phin = np.abs(phi)/np.abs(phi).max()
v_amb = 0.28; v_inj = pk[0]                         # ambient conv-wave amplitude vs launched peak
floor = (v_amb/v_inj)**2
# Fit ONLY the clean, fully-transited near-field decay. The run stops when the pulse reaches 0.1 R*, so
# its leading edge sits at x~0.08-0.1 with full instantaneous flux (Phi upturns there) — that is the
# pulse front just arriving, NOT decayed residual, and must be excluded from the decay fit.
x_edge = 0.075                                       # start of the leading-edge/arrival region
m = (xf > 0.004) & (xf < 0.065) & (phin > 1e-6)
sl, b = np.polyfit(xf[m], np.log(phin[m]), 1); L = -1/sl
xcross = -np.log(floor)*L
fig, ax = plt.subplots(figsize=(9,5), constrained_layout=True)
ax.semilogy(xf, phin, '-', color="tab:blue", lw=1.4, label=r"measured $\Phi(x)/\Phi_0$ (MHD, res×2)")
xe = np.linspace(0, Rstar_pole*1.02, 300)
ax.semilogy(xe, np.exp(b+sl*xe), '--', color="tab:red", lw=1.5, label=fr"exp. fit, $L={L:.4f}\,R_\star$")
ax.axhline(floor, color="tab:green", lw=1.6)
ax.axhspan(1e-8, floor, color="tab:green", alpha=0.10)
ax.text(Rstar_pole*0.45, floor*1.7, r"ambient convective-wave background ($v_\perp\!\sim\!0.28$ km/s)",
        fontsize=9, color="darkgreen")
ax.plot([xcross],[floor],'o',color="darkgreen",ms=7,zorder=5)
ax.annotate(f"signal buried in ambient\nat x≈{xcross:.3f} $R_\\star$", xy=(xcross,floor),
            xytext=(xcross+0.05,3e-2), fontsize=9, color="darkgreen",
            arrowprops=dict(arrowstyle="->",color="darkgreen"))
ax.axvspan(0, x_edge, color="tab:blue", alpha=0.05)
ax.axvspan(x_edge, x[-1], color="0.6", alpha=0.18)
ax.annotate("pulse front\n(just reached edge;\nnot residual)", xy=(0.096, 3e-4),
            xytext=(0.135, 6e-3), fontsize=8, color="0.30", ha="left",
            arrowprops=dict(arrowstyle="->", color="0.45"))
ax.axvline(Rstar_pole, color="k", ls=":", lw=1.4)
ax.text(Rstar_pole-0.008, 1e-6, f"pole\n(0.61 $R_\\star$,\n{Rstar_pole/xcross:.0f}× further)", ha="right", fontsize=9)
ax.set_xlim(0, Rstar_pole*1.02); ax.set_ylim(1e-7, 3)
ax.set_xlabel(r"lateral distance from footpoint  $x/R_\star$")
ax.set_ylabel(r"energy-flux $\Phi(x)/\Phi_0$")
ax.set_title(f"Fast-mode surface wave decays into the ambient background long before the pole\n"
             f"(res×2, 0.1 $R_\\star$ domain; $L$≈{L:.4f} $R_\\star$; buried by x≈{xcross:.2f} $R_\\star$)")
ax.legend(loc="lower left", fontsize=9)
plt.savefig(os.path.join(outd,f"flux_decay_{tag}.png"), dpi=140)
print(f"wrote flux_decay_{tag}.png (L={L:.4f} R*, xcross={xcross:.3f} R*, floor={floor:.1e})")
