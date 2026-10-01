import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import os
d=os.path.join(os.path.dirname(__file__),"..","..","output","plots","surfacecode")   # viz data (gitignored)
outd=os.path.join(os.path.dirname(__file__),"plots"); os.makedirs(outd,exist_ok=True)  # figures (tracked)
x=np.loadtxt(os.path.join(d,"viz_flux.txt"))[:,0]
phi=np.loadtxt(os.path.join(d,"viz_flux.txt"))[:,1]
transit=0.61   # ring->pole arc in R*

# ---- Plot 1: decay curve down to the ambient wave background ----
# physical floor: the ever-present convection-driven wave field (v_perp~0.28 km/s) relative to the
# MEASURED launched-pulse peak (from the trigger snapshot). Phi ~ v_perp^2 => floor = (0.28/v_peak)^2.
v_amb = 0.28
v_inj = float(np.abs(np.loadtxt(os.path.join(d,"viz_vx_1.txt"))).max())   # measured domain pulse peak [km/s]
floor=(v_amb/v_inj)**2                                   # ambient convective-wave background (self-consistent)
phin=np.abs(phi)/np.abs(phi).max()
# fit ONLY the clean exponential-decay region (above the numerical noise floor, past the injection)
m=(x>0.10*x.max())&(phin<3e-2)&(phin>1e-7)
sl,b=np.polyfit(x[m],np.log(phin[m]),1); L=-1/sl
xcross=-np.log(floor)*L                                 # x where the decay meets the ambient background
fig,ax=plt.subplots(figsize=(9,5),constrained_layout=True)
ax.semilogy(x,phin,'-',color="tab:blue",lw=1.4,label=r"measured $\Phi(x)/\Phi_0$ (MHD patch)")
xe=np.linspace(0,transit*1.03,200)
ax.semilogy(xe,np.exp(b+sl*xe),'--',color="tab:red",lw=1.5,label=fr"exp. fit, $L={L:.4f}\,R_\star$")
# ambient background band + crossing
ax.axhline(floor,color="tab:green",lw=1.6)
ax.axhspan(1e-6,floor,color="tab:green",alpha=0.10)
ax.text(transit*0.5,floor*1.6,r"ambient convective-wave background  ($v_\perp\!\sim\!0.28$ km/s)",
        fontsize=9,color="darkgreen")
ax.plot([xcross],[floor],'o',color="darkgreen",ms=7,zorder=5)
ax.annotate(f"accretion signal buried in\nambient background at x≈{xcross:.3f} $R_\\star$",
            xy=(xcross,floor),xytext=(xcross+0.06,3e-2),fontsize=9,color="darkgreen",
            arrowprops=dict(arrowstyle="->",color="darkgreen"))
ax.axvline(transit,color="k",ls=":",lw=1.4)
ax.text(transit-0.008,2e-5,f"pole\n(0.61 $R_\\star$,\n{transit/xcross:.0f}× further)",ha="right",fontsize=9)
ax.axvspan(0,x.max(),color="tab:blue",alpha=0.05)
ax.text(x.max()*0.5,1.5,"simulated patch",ha="center",fontsize=8,color="tab:blue")
ax.set_xlim(0,transit*1.03); ax.set_ylim(1e-6,3)
ax.set_xlabel(r"lateral distance from footpoint  $x / R_\star$")
ax.set_ylabel(r"energy-flux $\Phi(x)/\Phi_0$")
ax.set_title("Fast-mode surface wave decays into the ambient background long before the pole\n"
             f"(MHD patch, T Tauri; L≈{L:.4f} $R_\\star$; signal < ambient by x≈{xcross:.2f} $R_\\star$, "
             f"pole is {transit/xcross:.0f}× beyond)")
ax.legend(loc="lower left",fontsize=9)
plt.savefig(os.path.join(outd,"decay_to_pole.png"),dpi=140); print(f"wrote decay_to_pole.png (L={L:.4f} R*, xcross={xcross:.3f} R*, floor={floor:.1e})")

# ---- Plot 2: 2D v_x snapshots (pulse skimming and fading) ----
xc=np.loadtxt(os.path.join(d,"viz_x.txt")); zc=np.loadtxt(os.path.join(d,"viz_z.txt"))
snaps=[1,2,3]
try:
    ts=np.atleast_1d(np.loadtxt(os.path.join(d,"viz_snaptimes.txt")))
    def fmt(s): return f"t = {s:.0f} s" if s<600 else (f"t = {s/60:.0f} min" if s<7200 else f"t = {s/3600:.1f} h")
    labels=[fmt(ts[i]) for i in range(len(snaps))]
    labels[0]+="  (initial trigger)"
except Exception:
    labels=["initial trigger","mid","late"]
fig,axes=plt.subplots(3,1,figsize=(9,7),constrained_layout=True,sharex=True)
ext=[xc[0],xc[-1],zc[0],zc[-1]]
for ax,s,lb in zip(axes,snaps,labels):
    f=np.loadtxt(os.path.join(d,f"viz_vx_{s}.txt")).T          # (Nz,Nx)
    vpk=max(np.abs(f).max(),1e-3)                               # per-panel scale (each shows its own pulse)
    im=ax.imshow(f,origin="lower",aspect="auto",extent=ext,cmap="RdBu_r",
                 norm=TwoSlopeNorm(vmin=-vpk,vcenter=0,vmax=vpk))
    ax.set_ylabel(r"height $z/H_p$")
    ax.text(0.02,0.80,lb,transform=ax.transAxes,fontsize=10,bbox=dict(fc="white",alpha=0.7,ec="none"))
    ax.text(0.98,0.80,f"peak $|v_x|$ = {vpk:.2f} km/s",transform=ax.transAxes,fontsize=9,ha="right",
            color="firebrick",bbox=dict(fc="white",alpha=0.7,ec="none"))
    ax.axhline(0,color="k",lw=0.4,ls=":")
    plt.colorbar(im,ax=ax,shrink=0.9,label=r"$v_x$ [km/s]")
axes[-1].set_xlabel(r"lateral distance  $x/R_\star$")
axes[0].set_title("Fast-mode pulse skimming the surface: peak amplitude collapses "
                  "~60× (below the 0.28 km/s ambient) within ~0.04 $R_\\star$\n"
                  "(each panel auto-scaled; absolute times and peak-$|v_x|$ annotated)")
plt.savefig(os.path.join(outd,"pulse_snapshots.png"),dpi=140); print("wrote pulse_snapshots.png")
