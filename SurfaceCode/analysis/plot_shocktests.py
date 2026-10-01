#!/usr/bin/env python3
# plot_shocktests.py — appendix validation figures for the SurfaceCode solver.
#   shocktube_sod.png    : Sod (hydro) numerical vs the EXACT Riemann solution (rho, u, p).
#   shocktube_briowu.png : Brio-Wu (MHD) Nx=400 vs an Nx=1600 reference (rho, vx, p, By).
import numpy as np, matplotlib, os
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(__file__)
d    = os.path.join(here, "..", "..", "output", "shocktests")
outd = os.path.join(here, "plots"); os.makedirs(outd, exist_ok=True)

# ---------- exact Sod Riemann solution (Toro) ----------
def sod_exact(x, t, gamma=1.4, x0=0.5, WL=(1.0,0.0,1.0), WR=(0.125,0.0,0.1)):
    rhoL,uL,pL = WL; rhoR,uR,pR = WR
    aL=np.sqrt(gamma*pL/rhoL); aR=np.sqrt(gamma*pR/rhoR); g=gamma
    def f(p,rhoK,pK,aK):
        if p>pK:
            A=2/((g+1)*rhoK); B=(g-1)/(g+1)*pK; return (p-pK)*np.sqrt(A/(p+B))
        return 2*aK/(g-1)*((p/pK)**((g-1)/(2*g))-1)
    def fp(p,rhoK,pK,aK):
        if p>pK:
            A=2/((g+1)*rhoK); B=(g-1)/(g+1)*pK; return np.sqrt(A/(B+p))*(1-(p-pK)/(2*(B+p)))
        return 1/(rhoK*aK)*(p/pK)**(-(g+1)/(2*g))
    p=0.5*(pL+pR)
    for _ in range(100):
        gg=f(p,rhoL,pL,aL)+f(p,rhoR,pR,aR)+(uR-uL); gpp=fp(p,rhoL,pL,aL)+fp(p,rhoR,pR,aR)
        pn=p-gg/gpp
        if abs(pn-p)<1e-12: p=pn; break
        p=max(pn,1e-9)
    pstar=p; ustar=0.5*(uL+uR)+0.5*(f(pstar,rhoR,pR,aR)-f(pstar,rhoL,pL,aL))
    rho=np.zeros_like(x); u=np.zeros_like(x); pr=np.zeros_like(x)
    for i,xi in enumerate(x):
        S=(xi-x0)/t
        if S<=ustar:                                    # left of contact
            if pstar>pL:                                # left shock
                SL=uL-aL*np.sqrt((g+1)/(2*g)*pstar/pL+(g-1)/(2*g))
                if S<=SL: rho[i],u[i],pr[i]=rhoL,uL,pL
                else:
                    rho[i]=rhoL*((pstar/pL+(g-1)/(g+1))/((g-1)/(g+1)*pstar/pL+1)); u[i],pr[i]=ustar,pstar
            else:                                       # left rarefaction
                rsl=rhoL*(pstar/pL)**(1/g); asl=aL*(pstar/pL)**((g-1)/(2*g))
                SHL=uL-aL; STL=ustar-asl                 # head < tail
                if S<=SHL: rho[i],u[i],pr[i]=rhoL,uL,pL   # unshocked left state
                elif S<STL:                              # inside the rarefaction fan
                    u[i]=2/(g+1)*(aL+(g-1)/2*uL+S); a_=2/(g+1)*(aL+(g-1)/2*(uL-S))
                    rho[i]=rhoL*(a_/aL)**(2/(g-1)); pr[i]=pL*(a_/aL)**(2*g/(g-1))
                else: rho[i],u[i],pr[i]=rsl,ustar,pstar   # star-left (tail..contact)
        else:                                           # right of contact
            if pstar>pR:                                # right shock
                SR=uR+aR*np.sqrt((g+1)/(2*g)*pstar/pR+(g-1)/(2*g))
                if S>=SR: rho[i],u[i],pr[i]=rhoR,uR,pR
                else:
                    rho[i]=rhoR*((pstar/pR+(g-1)/(g+1))/((g-1)/(g+1)*pstar/pR+1)); u[i],pr[i]=ustar,pstar
            else:                                       # right rarefaction
                rsr=rhoR*(pstar/pR)**(1/g); asr=aR*(pstar/pR)**((g-1)/(2*g))
                SHR=uR+aR; STR=ustar+asr                 # head > tail
                if S>=SHR: rho[i],u[i],pr[i]=rhoR,uR,pR   # unshocked right state
                elif S>STR:                              # inside the rarefaction fan
                    u[i]=2/(g+1)*(-aR+(g-1)/2*uR+S); a_=2/(g+1)*(aR-(g-1)/2*(uR-S))
                    rho[i]=rhoR*(a_/aR)**(2/(g-1)); pr[i]=pR*(a_/aR)**(2*g/(g-1))
                else: rho[i],u[i],pr[i]=rsr,ustar,pstar   # star-right (contact..tail)
    return rho,u,pr

# ---------- Fig 1: Sod ----------
x,rho,u,p = np.loadtxt(os.path.join(d,"sod.txt")).T
xe = np.linspace(0,1,2000); re,ue,pe = sod_exact(xe, 0.2)
fig,axes=plt.subplots(1,3,figsize=(12,3.6),constrained_layout=True)
for ax,(yn,ye,lab) in zip(axes,[(rho,re,r"density $\rho$"),(u,ue,r"velocity $u$"),(p,pe,r"pressure $p$")]):
    ax.plot(xe,ye,'-',color="k",lw=1.4,label="exact",zorder=1)
    ax.plot(x,yn,'o',color="tab:blue",ms=2.6,mew=0,label="SurfaceCode (Nx=200)",zorder=2)
    ax.set_xlabel("x"); ax.set_ylabel(lab); ax.set_xlim(0,1)
axes[0].legend(loc="upper right",fontsize=8,framealpha=0.9)
fig.suptitle(r"Sod shock tube (hydro), $t=0.2$, $\gamma=1.4$ — SurfaceCode vs exact Riemann solution",fontsize=11)
plt.savefig(os.path.join(outd,"shocktube_sod.png"),dpi=140); print("wrote shocktube_sod.png")

# ---------- Fig 2: Brio-Wu ----------
x4,r4,v4,p4,By4 = np.loadtxt(os.path.join(d,"briowu_400.txt")).T
xR,rR,vR,pR,ByR = np.loadtxt(os.path.join(d,"briowu_ref.txt")).T
fig,axes=plt.subplots(2,2,figsize=(11,6.4),constrained_layout=True)
panels=[(r4,rR,r"density $\rho$"),(v4,vR,r"velocity $v_x$"),(p4,pR,r"pressure $p$"),(By4,ByR,r"transverse field $B_y$")]
for ax,(y4,yR,lab) in zip(axes.flat,panels):
    ax.plot(xR,yR,'-',color="k",lw=1.2,label="reference (Nx=1600)",zorder=1)
    ax.plot(x4,y4,'o',color="tab:orange",ms=2.6,mew=0,label="SurfaceCode (Nx=400)",zorder=2)
    ax.set_xlabel("x"); ax.set_ylabel(lab); ax.set_xlim(0,1)
axes[0,0].legend(loc="upper right",fontsize=8,framealpha=0.9)
fig.suptitle(r"Brio--Wu shock tube (MHD), $t=0.1$, $\gamma=2$ — SurfaceCode vs high-resolution reference",fontsize=11)
plt.savefig(os.path.join(outd,"shocktube_briowu.png"),dpi=140); print("wrote shocktube_briowu.png")
