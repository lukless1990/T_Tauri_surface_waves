#!/usr/bin/env python3
# make_paper_figs.py — regenerate all main-text + appendix figures of the surface-wave
# paper as A&A-sized vector PDFs into Paper/Surface_Wave_paper/figs/.
# Data: ../../../output/{paper_run_tau50,paper_run_0p1kG_tau50,paper_run_cont_tau50,shocktests,wavetests}
#       (run labels overridable via RUN_1KG / RUN_0P1KG / RUN_CONT — see below)
# (same data the exploratory plots_/*.png were made from; this is a re-render, not a re-run).
import os, glob, json
import numpy as np
import aastyle
aastyle.setup()
from aastyle import COL, TWOCOL, C
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm

here = os.path.dirname(os.path.abspath(__file__))
sc   = os.path.join(here, "..", "..")                      # SurfaceCode/
out  = os.path.join(here, "..", "Surface_Wave_paper", "figs")
odat = os.path.join(sc, "..", "output")
os.makedirs(out, exist_ok=True)
POLE  = 0.61        # ring->pole arc for theta_ring=35 deg [R*]
V_AMB = 0.28        # ambient convective-wave amplitude [km/s]

# Production runs. Default = the clump-consistent driver duration tau = 50 s (r_c ~ 1e-2 R*,
# paper Sect. 2.2); the tau = 150 s runs are kept as the driver-duration comparison of Sect. 4.3
# and can be re-selected with e.g.  RUN_1KG=paper_run_1kG RUN_0P1KG=paper_run_0p1kG RUN_CONT=paper_run_cont
R_1KG   = os.environ.get("RUN_1KG",   "paper_run_tau50")
R_0P1KG = os.environ.get("RUN_0P1KG", "paper_run_0p1kG_tau50")
R_CONT  = os.environ.get("RUN_CONT",  "paper_run_cont_tau50")
nums  = {}          # numbers quoted in the text -> figs/fig_numbers.json

def meta_of(run):
    d = os.path.join(odat, run)
    return dict(l.split("=") for l in open(os.path.join(d, "meta.txt")).read().split() if "=" in l)

def savefig(fig, name):
    fig.savefig(os.path.join(out, name))
    plt.close(fig)
    print("wrote", name)

# ================================================================ Fig: background
bg = np.loadtxt(os.path.join(sc, "background_35deg.txt"))
zb, rho, T, P, VA1, cs, Hp = bg[:,0], bg[:,2], bg[:,3], bg[:,4], bg[:,11], bg[:,12], bg[:,13]
Hp0 = 8.396e5                                              # Hp at z=0 [m] (paper-run meta)
zH  = zb / Hp0
msk = (zH >= -4.2) & (zH <= 3.2)                           # the simulated patch is [-4,3] Hp
fig, axes = plt.subplots(2, 1, figsize=(COL, 3.4), sharex=True)
ax = axes[0]
ax.semilogy(zH[msk], rho[msk], color=C["blue"])
ax.set_ylabel(r"$\rho\ \mathrm{[kg\,m^{-3}]}$", color=C["blue"])
ax.tick_params(axis="y", colors=C["blue"], which="both")
ax2 = ax.twinx()
ax2.plot(zH[msk], T[msk], color=C["red"])
ax2.set_ylabel(r"$T\ \mathrm{[K]}$", color=C["red"])
ax2.tick_params(axis="y", colors=C["red"], which="both")
ax = axes[1]
ax.semilogy(zH[msk], cs[msk]/1e3, color=C["gray"], label=r"$c_\mathrm{s}$")
ax.semilogy(zH[msk], VA1[msk]/1e3, color=C["blue"], label=r"$V_\mathrm{A}$ ($B_\star=1$ kG)")
ax.semilogy(zH[msk], 0.1*VA1[msk]/1e3, color=C["blue"], ls="--",
            label=r"$V_\mathrm{A}$ ($B_\star=0.1$ kG)")
ax.axvline(0, color="0.85", lw=0.7)
ax.set_xlabel(r"$z/H_p$")
ax.set_ylabel(r"signal speed [km s$^{-1}$]")
ax.legend(loc="upper left")
ax.set_xlim(-4.2, 3.2)
savefig(fig, "fig_background.pdf")

# ================================================================ Fig: pulse panels (2-col)
def load_run_panels(run, targets):
    d = os.path.join(odat, run)
    x = np.loadtxt(os.path.join(d, "grid_x.txt"))
    z = np.loadtxt(os.path.join(d, "grid_z.txt"))
    st = np.atleast_2d(np.loadtxt(os.path.join(d, "snap_times.txt")))
    tt = st[:,0]
    snaps = sorted(glob.glob(os.path.join(d, "snap_vx_*.txt")))
    sel = [int(np.argmin(np.abs(tt - tg))) for tg in targets]
    zmask = (z >= 0.0) & (z <= 2.0)
    fields = [np.loadtxt(snaps[k]).T[zmask, :] for k in sel]
    return x, z[zmask], tt, sel, fields

targets = [225.0, 5400.0, 14400.0]
runs = [(R_1KG, r"$B_\star = 1$ kG"), (R_0P1KG, r"$B_\star = 0.1$ kG")]
data = {r: load_run_panels(r, targets) for r, _ in runs}
vmax = max(np.abs(f).max() for r, _ in runs for f in data[r][4])
norm = SymLogNorm(linthresh=0.03, linscale=0.8, vmin=-vmax, vmax=vmax, base=10)

fig, axes = plt.subplots(3, 2, figsize=(TWOCOL, 3.6), sharex=True, sharey=True)
def fmt_t(s): return f"{s:.0f} s" if s < 600 else (f"{s/60:.0f} min" if s < 7200 else f"{s/3600:.1f} h")
for col, (run, blab) in enumerate(runs):
    x, zc, tt, sel, fields = data[run]
    ext = [x[0], x[-1], zc[0], zc[-1]]
    dxx = x[1]-x[0]; nw = max(1, int(round(0.002/dxx)))
    for row, (k, f) in enumerate(zip(sel, fields)):
        ax = axes[row, col]
        im = ax.imshow(f, origin="lower", aspect="auto", extent=ext, cmap="RdBu_r",
                       norm=norm, interpolation="nearest", rasterized=True)
        env = np.abs(f).max(axis=0)
        idx = np.where(env > 0.006)[0]
        if len(idx):
            i_f = idx[-1]; x_f = x[i_f]
            fv = env[max(0, i_f-nw):i_f+1].max()
            ax.axvline(x_f, color="0.25", ls="--", lw=0.7)
            ax.text(0.985, 0.76, rf"front: $|v_x|$ = {fv:.2f} km s$^{{-1}}$",
                    transform=ax.transAxes, fontsize=6.2, ha="right", color="0.15",
                    bbox=dict(fc="white", alpha=0.75, ec="none", pad=1.2))
        lbl = "launch" if tt[k] < 400 else fmt_t(tt[k])
        ax.text(0.015, 0.76, lbl, transform=ax.transAxes, fontsize=6.8,
                bbox=dict(fc="white", alpha=0.75, ec="none", pad=1.2))
        if row == 0: ax.set_title(blab, fontsize=8.5)
        if col == 0: ax.set_ylabel(r"$z/H_p$")
        ax.set_ylim(0, 2)
for ax in axes[-1]: ax.set_xlabel(r"$x/R_\star$")
cb = fig.colorbar(im, ax=axes, shrink=0.92, pad=0.012)
cb.set_label(r"$v_x$ [km s$^{-1}$]")
savefig(fig, "fig_pulse_panels.pdf")

# ================================================================ Fig: wake zoom (1-col, 2 stacked)
dW  = os.path.join(odat, R_1KG)
xw  = np.loadtxt(os.path.join(dW, "grid_x.txt"))
zw2 = np.loadtxt(os.path.join(dW, "grid_z.txt"))
stw = np.atleast_2d(np.loadtxt(os.path.join(dW, "snap_times.txt")))
snw = sorted(glob.glob(os.path.join(dW, "snap_vx_*.txt")))
HpR = 8.396e5 / 1.391e9                                    # Hp in R* units
CF  = 8.3                                                  # photospheric c_f [km/s], 1 kG
wtimes = [5400.0, 7200.0, 9000.0]
wks = [int(np.argmin(np.abs(stw[:,0] - t))) for t in wtimes]
zmk = (zw2 >= 0.0) & (zw2 <= 2.0)
iz1 = int(np.argmin(np.abs(zw2 - 1.0)))
fronts, profs = [], []
for k in wks:                                              # front = same env>0.006 criterion as fig_pulse
    F = np.loadtxt(snw[k])
    env = np.abs(F[:, zmk]).max(axis=1)
    idx = np.where(env > 0.006)[0]
    fronts.append(xw[idx[-1]] if len(idx) else np.nan)
    profs.append(F[:, iz1])
fig, ax = plt.subplots(figsize=(COL, 2.1))
ax.axhline(0, color="0.85", lw=0.6)
for t, fx, pv, cc in zip(wtimes, fronts, profs, (C["blue"], C["orange"], C["green"])):
    ax.plot(xw, pv, color=cc, lw=0.9, label=f"{t/60:.0f} min")
    ax.axvline(fx, color=cc, ls=":", lw=0.7)
ax.set_xlim(0, 0.06); ax.set_ylim(-0.13, 0.13)
ax.set_xlabel(r"$x/R_\star$"); ax.set_ylabel(r"$v_x(z\!=\!H_p)$ [km s$^{-1}$]")
ax.legend(loc="upper left")
savefig(fig, "fig_wake.pdf")

# wake numbers quoted in the text
dxg = xw[1] - xw[0]; nwd = max(1, int(round(0.002/dxg)))
v90 = profs[0]; i_f = int(np.argmin(np.abs(xw - fronts[0])))
sg = np.sign(v90[:i_f+1]); zcx = np.where(np.diff(sg) != 0)[0]
zcx = zcx[xw[zcx] > 0.2*fronts[0]]
lam = 2.0*np.median(np.diff(xw[zcx]))
w = (xw > 0.010) & (xw < 0.030)                            # rear-tail pattern speed 90->150 min by cross-correlation
a, b = profs[0][w], profs[2][w]
sh = np.arange(-160, 161)
ccr = [np.corrcoef(np.roll(b, -s)[80:-80], a[80:-80])[0, 1] for s in sh]
dxsh = sh[int(np.argmax(ccr))]*dxg
dtw = stw[wks[2], 0] - stw[wks[0], 0]
nums["wake_lambda_Hp"]        = round(lam/HpR, 1)
nums["wake_tail_over_cf"]     = round(dxsh/(CF*1e3*dtw/1.391e9), 2)
nums["wake_crest90_kms"]      = round(float(np.abs(v90[max(0, i_f-nwd):i_f+1]).max()), 2)
nums["wake_maxlobe90_kms"]    = round(float(np.abs(v90[:i_f-nwd]).max()), 2)

# ================================================================ Fig: flux decay (1-col, 2 stacked)
def load_flux(run):
    d = os.path.join(odat, run)
    x, phi = np.loadtxt(os.path.join(d, "flux_final.txt")).T
    pk = np.atleast_2d(np.loadtxt(os.path.join(d, "snap_times.txt")))[:,1]
    return x, np.abs(phi)/np.abs(phi).max(), pk[0]

def ambient_crossing(x, phin, floor):
    below = np.where((x > 0.02) & (phin < floor))[0]
    if not len(below): return np.nan
    i = below[0]
    y0, y1 = np.log(phin[i-1]), np.log(phin[i])
    return x[i-1] + (np.log(floor)-y0)*(x[i]-x[i-1])/(y1-y0)

def amp_crossing(run):
    # last x where the max-over-time-and-height |vx| still exceeds the ambient amplitude
    d = os.path.join(odat, run)
    x = np.loadtxt(os.path.join(d, "grid_x.txt"))
    A = np.zeros_like(x)
    for s in sorted(glob.glob(os.path.join(d, "snap_vx_*.txt"))):
        A = np.maximum(A, np.abs(np.loadtxt(s)).max(axis=1))
    idx = np.where(A > V_AMB)[0]
    return x[idx[-1]] if len(idx) else np.nan

x1, p1, v1 = load_flux(R_1KG)
x0, p0, v0 = load_flux(R_0P1KG)
xamp1 = amp_crossing(R_1KG)
xamp0 = amp_crossing(R_0P1KG)
floor = (V_AMB/10.0)**2
# exponential fit on the clean transited part of the 1 kG run
m = (x1 > 0.004) & (x1 < 0.065) & (p1 > 1e-6)
sl, b = np.polyfit(x1[m], np.log(p1[m]), 1); L1 = -1/sl
xa1 = ambient_crossing(x1, p1, floor)
xa0 = ambient_crossing(x0, p0, floor)
L0eff = xa0/np.log(1/floor)
nums.update(L_fit_1kG=L1, x_amb_1kG=xa1, x_amb_0p1kG=xa0, L_eff_0p1kG=L0eff,
            floor=floor, pole_over_xamb_1kG=POLE/xa1, pole_over_xamb_0p1kG=POLE/xa0,
            x_ampcross_1kG=xamp1, x_ampcross_0p1kG=xamp0)

fig, axA = plt.subplots(figsize=(COL, 2.4))
axA.axhline(floor, color=C["green"], lw=1.0)
axA.axhspan(1e-9, floor, color=C["green"], alpha=0.09)
axA.semilogy(x1, p1, color=C["blue"], label=r"$B_\star=1$ kG")
axA.semilogy(x0, p0, color=C["red"], label=r"$B_\star=0.1$ kG")
xe = np.linspace(0, 0.1, 200)
axA.semilogy(xe, np.exp(b+sl*xe), "--", color=C["blue"], lw=0.9)
for xa, cc in ((xa1, C["blue"]), (xa0, C["red"])):
    axA.plot([xa], [floor], "o", color=cc, ms=5, mec="k", mew=0.5, zorder=5)
for xam, cc in ((xamp1, C["blue"]), (xamp0, C["red"])):
    axA.plot([xam], [floor], "v", mfc="none", mec=cc, mew=1.0, ms=6, zorder=5)
axA.text(0.003, floor*0.32, r"$\Phi_\mathrm{amb}$", fontsize=7.5, color=C["green"])
axA.set_xlim(0, 0.1); axA.set_ylim(2e-4, 2.5)
axA.set_xlabel(r"lateral distance from footpoint  $x/R_\star$")
axA.set_ylabel(r"$\Phi(x)/\Phi_0$")
axA.legend(loc="upper right")
savefig(fig, "fig_flux_decay.pdf")

# ================================================================ Fig: continuous driving (1-col)
dc = os.path.join(odat, R_CONT)
xc, Fc = np.loadtxt(os.path.join(dc, "flux_steady.txt")).T
def fit_L(x, f, xhi=0.065):
    fn = np.abs(f)/np.abs(f).max()
    mm = (x > 0.004) & (x < xhi) & (fn > 1e-6)
    s, bb = np.polyfit(x[mm], np.log(fn[mm]), 1)
    return fn, -1/s, s, bb
Fcn, Lc, slc, bc = fit_L(xc, Fc)
Fpn, Lp, slp, bp = fit_L(x1, p1)
nums.update(L_cont=Lc, L_pulse=Lp)

zc_ = np.loadtxt(os.path.join(dc, "grid_z.txt")); xg = np.loadtxt(os.path.join(dc, "grid_x.txt"))
st = np.atleast_2d(np.loadtxt(os.path.join(dc, "snap_times.txt"))); tt = st[:,0]
snaps = sorted(glob.glob(os.path.join(dc, "snap_vx_*.txt")))
k4 = int(np.argmin(np.abs(tt-14400.0)))
zmask = (zc_ >= 0) & (zc_ <= 2.0)
f4 = np.loadtxt(snaps[k4]).T[zmask, :]
fig, axS = plt.subplots(figsize=(COL, 1.75))
vm = np.abs(f4).max()
im = axS.imshow(f4, origin="lower", aspect="auto",
                extent=[xg[0], xg[-1], zc_[zmask][0], zc_[zmask][-1]],
                cmap="RdBu_r", norm=SymLogNorm(linthresh=0.03, linscale=0.8,
                vmin=-vm, vmax=vm, base=10), interpolation="nearest", rasterized=True)
axS.set_ylim(0, 2); axS.set_ylabel(r"$z/H_p$")
axS.set_xlim(0, 0.1)
axS.set_xlabel(r"lateral distance from footpoint  $x/R_\star$")
axS.text(0.015, 0.80, "continuous driving, $t$ = 4 h", transform=axS.transAxes,
         fontsize=6.8, bbox=dict(fc="white", alpha=0.75, ec="none", pad=1.2))
cb = fig.colorbar(im, ax=axS, pad=0.02, fraction=0.08, aspect=13)
cb.set_label(r"$v_x$ [km s$^{-1}$]", fontsize=6.5, labelpad=1)
cb.ax.tick_params(labelsize=5.8)
savefig(fig, "fig_cont.pdf")

# ================================================================ Appendix: shock tubes (2-col)
dS = os.path.join(odat, "shocktests")
import importlib.util
spec = importlib.util.spec_from_file_location("pst", os.path.join(sc, "plot_shocktests.py"))
# reuse the exact Sod Riemann solver from the exploratory plotter without executing its __main__
src = open(os.path.join(sc, "plot_shocktests.py")).read()
sod_src = src[src.index("def sod_exact"):src.index("# ---------- Fig 1")]
ns_ = {"np": np}; exec(sod_src, ns_); sod_exact = ns_["sod_exact"]

xs, rs, us, ps = np.loadtxt(os.path.join(dS, "sod.txt")).T
xe_ = np.linspace(0, 1, 2000); re_, ue_, pe_ = sod_exact(xe_, 0.2)
x4, r4, v4, p4, By4 = np.loadtxt(os.path.join(dS, "briowu_400.txt")).T
xR, rR, vR, pR, ByR = np.loadtxt(os.path.join(dS, "briowu_ref.txt")).T

fig, axes = plt.subplots(2, 4, figsize=(TWOCOL, 3.3))
for ax, (yn, ye, lab) in zip(axes[0, :3], [(rs, re_, r"$\rho$"), (us, ue_, r"$u$"), (ps, pe_, r"$p$")]):
    ax.plot(xe_, ye, "-", color="k", lw=0.9, label="exact", zorder=1)
    ax.plot(xs, yn, "o", color=C["blue"], ms=1.7, mew=0, label=r"$N_x=200$", zorder=2)
    ax.set_xlabel("$x$"); ax.set_ylabel(lab); ax.set_xlim(0, 1)
axes[0, 3].axis("off")
h, l = axes[0, 0].get_legend_handles_labels()
axes[0, 3].legend(h, l, loc="center left", title="Sod (hydro)\n$t=0.2$, $\\gamma=1.4$",
                  title_fontsize=7.5)
for ax, (y4, yR, lab) in zip(axes[1], [(r4, rR, r"$\rho$"), (v4, vR, r"$v_x$"),
                                       (p4, pR, r"$p$"), (By4, ByR, r"$B_y$")]):
    ax.plot(xR, yR, "-", color="k", lw=0.9, label=r"$N_x=1600$", zorder=1)
    ax.plot(x4, y4, "o", color=C["orange"], ms=1.7, mew=0, label=r"$N_x=400$", zorder=2)
    ax.set_xlabel("$x$"); ax.set_ylabel(lab); ax.set_xlim(0, 1)
axes[1, 0].legend(loc="upper right", title="Brio--Wu (MHD)\n$t=0.1$, $\\gamma=2$",
                  title_fontsize=6.5)
savefig(fig, "figA_shocktubes.pdf")

# ================================================================ Appendix: wave tests
dW = os.path.join(odat, "wavetests")
H = 1.0
z, un, ua, Bn, Ba = np.loadtxt(os.path.join(dW, "test1_mhd.txt")).T
V = np.max(np.abs(ua/np.exp(z/(2*H))))
cv = np.loadtxt(os.path.join(dW, "test1_convergence.txt"))
fig, axes = plt.subplots(1, 2, figsize=(TWOCOL, 2.4))
ax = axes[0]
ax.plot(z, ua/V, "-", color="k", lw=0.9, label="exact", zorder=1)
ax.plot(z[::6], un[::6]/V, "o", color=C["blue"], ms=1.8, mew=0,
        label=r"SurfaceCode (512 c/$\lambda$)", zorder=2)
ax.plot(z, np.exp(z/(2*H)), "--", color="0.6", lw=0.7, label=r"$\pm e^{z/2H}$")
ax.plot(z, -np.exp(z/(2*H)), "--", color="0.6", lw=0.7)
ax.set_xlabel("$z/H$"); ax.set_ylabel(r"$u_z/V$")
ax.legend(loc="upper left")
ax = axes[1]
for B00, lab, cc, mk in ((0.0, "hydro", C["blue"], "o"), (1.0, r"MHD ($\beta=2$)", C["orange"], "s")):
    mm = cv[:,0] == B00
    ax.loglog(cv[mm,1], cv[mm,2], mk+"-", color=cc, ms=3, label=lab)
cp = np.array([32., 512.])
ax.loglog(cp, 0.35*(cp/32)**-1, ":", color="0.6", lw=0.8)
ax.loglog(cp, 0.35*(cp/32)**-2, ":", color="0.6", lw=0.8)
ax.axhline(0.02, color="k", lw=0.6, ls="--")
ax.text(35, 0.023, "2% criterion", fontsize=6.2)
ax.set_xlabel(r"cells per wavelength"); ax.set_ylabel(r"$\varepsilon(u_z)$")
ax.legend(loc="lower left")
savefig(fig, "figA_wavetest1.pdf")

rows = [l.split() for l in open(os.path.join(dW, "test2_numdiss.txt"))]
modes = np.array([r[0].strip('"') for r in rows])
dat = np.array([[float(v) for v in r[1:]] for r in rows])
Rstar = 2*6.957e8; Lmeas = 0.0062*Rstar
# Scale for the mapping is the SHOCKED pulse, not the driver: the launch wavelength
# (2*tau*cf = 830 km at tau=50 s) is forgotten within x_sh, and the N-wave lobe width
# measured from the snapshots is ~1.2e3 km for both drivers => lam_eff ~ 2.4e3 km.
lam_sci = 2400e3
fig, ax = plt.subplots(figsize=(COL, 2.5))
for mode, lab, cc, mk, mfc in (("fast", r"fast ($\mathbf{k}\perp\mathbf{B}$)", C["blue"], "o", C["blue"]),
                               ("alfven", r"Alfvén ($\mathbf{k}\parallel\mathbf{B}$)", C["orange"], "s", "none")):
    mm = modes == mode
    ax.loglog(dat[mm,0], dat[mm,2], mk, ls="-", color=cc, ms=3.5, mfc=mfc, label=lab)
ax.axhline(Lmeas/lam_sci, color="k", ls="--", lw=0.8)
ax.text(8.5, 1.2*Lmeas/lam_sci, r"measured $L(\Phi)$ in units of $\lambda_\mathrm{eff}$", fontsize=6.2)
for cpl, tag in ((48, "production"), (192, "finest")):
    ax.axvline(cpl, color="0.65", lw=0.6, ls=":")
    ax.text(cpl*1.05, 0.55, tag, fontsize=6.2, color="0.4")
ax.set_xlabel(r"cells per wavelength"); ax.set_ylabel(r"$L_\mathrm{num}(\Phi)/\lambda$")
ax.legend(loc="upper left")
savefig(fig, "figA_wavetest2.pdf")

with open(os.path.join(out, "fig_numbers.json"), "w") as f:
    json.dump({k: float(v) for k, v in nums.items()}, f, indent=1)
print(json.dumps({k: round(float(v), 5) for k, v in nums.items()}, indent=0))
