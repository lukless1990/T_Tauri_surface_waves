#!/usr/bin/env python3
# plot_wavetests.py — appendix verification figures for the wave tests
# (single-fluid analogs of Popescu Braileanu et al. 2019, A&A 627, A25).
#   wavetest_stratified.png : test 1 — stratified fast wave vs the exact solution
#                             (profile at 512 c/λ + convergence ladder).
#   wavetest_numdiss.png    : test 2 — numerical-dissipation floor L_num(Φ) vs
#                             cells/λ, against the measured L(Φ)=0.006 R⋆.
import numpy as np, matplotlib, os
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(__file__)
d    = os.path.join(here, "..", "..", "output", "wavetests")
outd = os.path.join(here, "plots"); os.makedirs(outd, exist_ok=True)

H=1.0; lam=H/4

# ---------- Fig 1: test 1 (stratified wave) ----------
fig,axes=plt.subplots(1,3,figsize=(13,3.8),constrained_layout=True)

z,un,ua,Bn,Ba = np.loadtxt(os.path.join(d,"test1_mhd.txt")).T
V=np.max(np.abs(ua/np.exp(z/(2*H))))          # driving amplitude, inferred from the exact profile
ax=axes[0]
ax.plot(z,ua/V,'-',color="k",lw=1.2,label="exact",zorder=1)
ax.plot(z[::6],un[::6]/V,'o',color="tab:blue",ms=2.4,mew=0,label="SurfaceCode (512 c/$\\lambda$)",zorder=2)
ax.plot(z, np.exp(z/(2*H)),'--',color="0.5",lw=0.9,label=r"$\pm e^{z/2H}$")
ax.plot(z,-np.exp(z/(2*H)),'--',color="0.5",lw=0.9)
ax.set_xlabel("z / H"); ax.set_ylabel(r"$u_z\,/\,V$")
ax.set_title(r"MHD ($\beta=2$) fast wave, stationary state",fontsize=10)
ax.legend(loc="upper left",fontsize=8,framealpha=0.9)

ax=axes[1]
ax.plot(z,(un-ua)/(V*np.exp(z/(2*H))),'-',color="tab:red",lw=0.8)
ax.set_xlabel("z / H"); ax.set_ylabel(r"residual $(u_z-u_z^{\rm ex})\,/\,V e^{z/2H}$")
ax.set_title("carrier-normalized residual (512 c/$\\lambda$)",fontsize=10)
ax.axhline(0,color="k",lw=0.6)

cv=np.loadtxt(os.path.join(d,"test1_convergence.txt"))
ax=axes[2]
for B00,lab,c,mk in ((0.0,"hydro (acoustic-gravity)","tab:blue","o"),(1.0,r"MHD ($\beta=2$, fast)","tab:orange","s")):
    m=cv[:,0]==B00
    ax.loglog(cv[m,1],cv[m,2],mk+"-",color=c,ms=4,label=lab)
cp=np.array([32.,512.])
ax.loglog(cp,0.35*(cp/32)**-1,':',color="0.5",lw=1,label=r"$\propto N^{-1}$ / $N^{-2}$")
ax.loglog(cp,0.35*(cp/32)**-2,':',color="0.5",lw=1)
ax.axhline(0.02,color="k",lw=0.8,ls="--"); ax.text(34,0.022,"2% criterion",fontsize=8)
ax.set_xlabel(r"cells per wavelength"); ax.set_ylabel(r"$\varepsilon(u_z)$ (Eq. 48 norm)")
ax.set_title("convergence to the exact solution",fontsize=10)
ax.legend(loc="lower left",fontsize=8,framealpha=0.9)
fig.suptitle("Verification test 1: linear fast wave in a stratified atmosphere "
             r"($\rho_0,p_0\propto e^{-z/H}$, $B_{x0}\propto e^{-z/2H}$) vs the exact solution",fontsize=11)
plt.savefig(os.path.join(outd,"wavetest_stratified.png"),dpi=140); print("wrote wavetest_stratified.png")

# ---------- Fig 2: test 2 (numerical dissipation) ----------
rows=[l.split() for l in open(os.path.join(d,"test2_numdiss.txt"))]
modes=np.array([r[0].strip('"') for r in rows]); dat=np.array([[float(v) for v in r[1:]] for r in rows])
Rstar=2*6.957e8; Lmeas=0.006*Rstar; lam_sci=2489e3   # m (2 tau c_f, 1 kG run)
fig,ax=plt.subplots(figsize=(6.2,4.4),constrained_layout=True)
for mode,lab,c,mk,mfc in (("fast",r"fast mode ($k\perp B$, science-pulse carrier)","tab:blue","o","tab:blue"),
                          ("alfven",u"Alfvén mode ($k\\parallel B$)","tab:orange","s","none")):
    m=modes==mode
    ax.loglog(dat[m,0],dat[m,2],mk,ls="-",color=c,ms=5,mfc=mfc,label=lab)
ax.axhline(Lmeas/lam_sci,color="k",ls="--",lw=1)
ax.text(9,1.15*Lmeas/lam_sci,r"measured $L(\Phi)=0.006\,R_\star$ (in units of $\lambda_{\rm pulse}$)",fontsize=8)
for cpl,tag in ((50,"production\n dx=50 km"),(191,"finest ladder\n dx=13 km")):
    ax.axvline(cpl,color="0.6",lw=0.8,ls=":")
    ax.text(cpl*1.04,0.6,tag,fontsize=7.5,color="0.35")
ax.set_xlabel(r"cells per wavelength"); ax.set_ylabel(r"$L_{\rm num}(\Phi)\;/\;\lambda$")
ax.set_title("Verification test 2: numerical-dissipation floor\n(ideal eigenmodes — analytic damping = 0)",fontsize=10)
ax.legend(loc="upper left",fontsize=8,framealpha=0.9)
plt.savefig(os.path.join(outd,"wavetest_numdiss.png"),dpi=140); print("wrote wavetest_numdiss.png")
