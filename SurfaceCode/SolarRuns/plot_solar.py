#!/usr/bin/env python3
"""plot_solar.py — figures for the SolarRuns Moreton-wave simulation.

Usage: python3 plot_solar.py run_<label> [outdir]
Writes PNGs into <outdir> (default: plots_<label>/ next to the run dir).

Figures:
  fig1_background.png   T, rho, c_s, V_A, V_fast vs height (the Moreton waveguide)
  fig2_snapshots.png    v_z maps at several times (coronal front + chromospheric skirt)
  fig3_stackplot.png    distance-time stack of chromospheric v_z (the "Halpha" observable)
                        + coronal dp/p stack (the "EIT" front)
  fig4_kinematics.png   front position & speed vs time, against the observed
                        Moreton range (500-2000 km/s, decelerating; Warmuth 2004/2015)
"""
import sys, os, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rundir = sys.argv[1] if len(sys.argv) > 1 else "run_preview"
rundir = rundir.rstrip("/")
label  = os.path.basename(rundir).replace("run_", "")
outdir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(rundir) or ".", f"plots_{label}")
os.makedirs(outdir, exist_ok=True)

meta = {}
with open(os.path.join(rundir, "meta.txt")) as f:
    for ln in f:
        k, v = ln.strip().split("=")
        try: meta[k] = float(v)
        except ValueError: meta[k] = v

bg   = np.loadtxt(os.path.join(rundir, "background.txt"))
xs   = np.loadtxt(os.path.join(rundir, "grid_x_Mm.txt"))      # subsampled x [Mm]
zs   = np.loadtxt(os.path.join(rundir, "grid_z_Mm.txt"))      # z [Mm]
snapt= np.atleast_1d(np.loadtxt(os.path.join(rundir, "snap_times.txt")))
trt  = np.loadtxt(os.path.join(rundir, "trace_times.txt"))
trx  = np.loadtxt(os.path.join(rundir, "trace_x_Mm.txt"))     # full-res x [Mm]
trvz = np.loadtxt(os.path.join(rundir, "trace_chromo_vz.txt"))# [t, x] km/s
trdp = np.loadtxt(os.path.join(rundir, "trace_corona_dpp.txt"))

x_fl = meta["x_fl_Mm"]

# ---------------------------------------------------------------- fig 1: background
z, rho, p, T, cs, va, vf = bg[:,0]/1e6, bg[:,1], bg[:,2], bg[:,3], bg[:,4]/1e3, bg[:,5]/1e3, bg[:,6]/1e3
fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
ax[0].semilogx(T, z, "k-", lw=2)
ax[0].set_xlabel("T [K]"); ax[0].set_ylabel("z [Mm]"); ax[0].set_title("temperature")
axr = ax[0].twiny(); axr.semilogx(rho, z, "C0--", lw=1.5); axr.set_xlabel(r"$\rho$ [kg m$^{-3}$]", color="C0")
ax[1].semilogx(cs, z, "C1-", lw=2, label="$c_s$")
ax[1].semilogx(va, z, "C2-", lw=2, label="$V_A$")
ax[1].semilogx(vf, z, "k--", lw=2, label=r"$V_{f\perp}=\sqrt{c_s^2+V_A^2}$")
ax[1].set_xlabel("speed [km/s]"); ax[1].set_title("wave speeds (the Moreton waveguide)")
ax[1].legend(loc="center right"); ax[1].grid(alpha=0.3)
for a in ax: a.axhline(0, color="gray", lw=0.5); a.axhline(meta["z_ha_Mm"], color="C3", lw=0.8, ls=":")
fig.suptitle(f"solar background  (B$_0$={meta['B0_G']:.0f} G)")
fig.tight_layout(); fig.savefig(os.path.join(outdir, "fig1_background.png"), dpi=150); plt.close(fig)

# ---------------------------------------------------------------- fig 2: vz snapshots
snaps = sorted(glob.glob(os.path.join(rundir, "snap_vz_*.txt")))
pick  = [i for i in (2, 4, 7, 10, len(snaps)-1) if 0 <= i < len(snaps)]
pick  = sorted(set(pick))
fig, axs = plt.subplots(len(pick), 1, figsize=(11, 2.1*len(pick)+1), sharex=True)
axs = np.atleast_1d(axs)
for a, i in zip(axs, pick):
    v = np.loadtxt(snaps[i]).T                      # [z, x]
    vm = max(0.5, np.percentile(np.abs(v), 99.5))
    im = a.pcolormesh(xs, zs, v, cmap="RdBu_r", vmin=-vm, vmax=vm, shading="auto")
    a.set_ylabel("z [Mm]")
    a.text(0.995, 0.93, f"t = {snapt[i]:.0f} s", ha="right", va="top",
           transform=a.transAxes, fontsize=9, bbox=dict(fc="w", alpha=0.8, ec="none"))
    plt.colorbar(im, ax=a, label="$v_z$ [km/s]", pad=0.01)
    a.axhline(meta["z_ha_Mm"], color="k", lw=0.6, ls=":")
axs[-1].set_xlabel("x [Mm]")
axs[0].set_title(f"flare-driven fast wave: $v_z$  (run {label}; dotted = chromospheric trace height)")
fig.tight_layout(); fig.savefig(os.path.join(outdir, "fig2_snapshots.png"), dpi=150); plt.close(fig)

# ------------------------------------------------------- fig 3: distance-time stacks
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
vm = max(0.2, np.percentile(np.abs(trvz), 99.5))
im0 = ax[0].pcolormesh(trx, trt, trvz, cmap="RdBu_r", vmin=-vm, vmax=vm, shading="auto")
plt.colorbar(im0, ax=ax[0], label="$v_z$ [km/s]")
ax[0].set_title(f"chromospheric $v_z$ at z={meta['z_ha_Mm']:.1f} Mm  (Moreton / H$\\alpha$ proxy)")
dm = max(0.05, np.percentile(np.abs(trdp), 99.5))
im1 = ax[1].pcolormesh(trx, trt, trdp, cmap="PuOr_r", vmin=-dm, vmax=dm, shading="auto")
plt.colorbar(im1, ax=ax[1], label=r"$\Delta p/p_0$")
ax[1].set_title(f"coronal $\\Delta p/p_0$ at z={meta['z_eit_Mm']:.0f} Mm  (EUV-wave proxy)")
for a in ax:
    a.set_xlabel("x [Mm]"); a.axvline(x_fl, color="k", lw=0.6, ls="--")
ax[0].set_ylabel("t [s]")
for vref, lab in ((500, "500 km/s"), (1000, "1000 km/s")):
    tt = np.linspace(0, trt[-1], 50)
    ax[0].plot(x_fl + vref*1e-3*tt, tt, "k:", lw=0.8)
for a in ax:
    a.set_xlim(trx[0], trx[-1]); a.set_ylim(trt[0], trt[-1])
fig.tight_layout(); fig.savefig(os.path.join(outdir, "fig3_stackplot.png"), dpi=150); plt.close(fig)

# ------------------------------------------------------------- fig 4: front kinematics
def track_front(stack, x, thresh_frac=0.10, floor=0.05):
    """rightmost x where |signal| exceeds max(thresh_frac*row peak, floor)"""
    pos = np.full(stack.shape[0], np.nan)
    for k in range(stack.shape[0]):
        row = np.abs(stack[k]); pk = row.max()
        th = max(thresh_frac*pk, floor)
        idx = np.where((row > th) & (x > x_fl))[0]
        if len(idx) and pk > floor: pos[k] = x[idx[-1]]
    return pos

pos_ha  = track_front(trvz, trx, floor=0.05)     # km/s floor for vz
pos_eit = track_front(trdp, trx, floor=0.02)     # dp/p floor
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].plot(trt, pos_ha,  "C3.", ms=3, label="chromospheric front (Moreton)")
ax[0].plot(trt, pos_eit, "C0.", ms=3, label="coronal front (EUV)")
ax[0].set_xlabel("t [s]"); ax[0].set_ylabel("front position x [Mm]"); ax[0].legend(); ax[0].grid(alpha=0.3)
ax[0].set_title("front position")
# speeds: finite difference over a sliding window
def speed(tt, pp, w=8):
    v = np.full(len(tt), np.nan)
    for k in range(w, len(tt)-w):
        if np.isfinite(pp[k-w]) and np.isfinite(pp[k+w]) and pp[k+w] > pp[k-w]:
            v[k] = (pp[k+w]-pp[k-w])/(tt[k+w]-tt[k-w])*1e3   # Mm/s -> km/s
    return v
ax[1].plot(trt, speed(trt, pos_ha),  "C3-", lw=1.5, label="Moreton front speed")
ax[1].plot(trt, speed(trt, pos_eit), "C0-", lw=1.5, label="coronal front speed")
ax[1].axhspan(500, 2000, color="gray", alpha=0.15, label="observed Moreton range\n(Warmuth 2004/2015)")
vfc = meta["cf_cor_kms"]
ax[1].axhline(vfc, color="k", ls="--", lw=1, label=f"ambient $V_f$ = {vfc:.0f} km/s")
ax[1].set_xlabel("t [s]"); ax[1].set_ylabel("speed [km/s]"); ax[1].set_ylim(0, 2500)
ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3); ax[1].set_title("front speed (decelerating shock)")
fig.tight_layout(); fig.savefig(os.path.join(outdir, "fig4_kinematics.png"), dpi=150); plt.close(fig)

# quick numbers for the console / PDF
msk = np.isfinite(pos_ha) & (trt > 60) & (trt < 0.7*trt[-1])
if msk.sum() > 10:
    vfit = np.polyfit(trt[msk], pos_ha[msk], 1)[0]*1e3
    print(f"mean chromospheric front speed (60 s < t < {0.7*trt[-1]:.0f} s): {vfit:.0f} km/s")
print(f"figures -> {outdir}")
