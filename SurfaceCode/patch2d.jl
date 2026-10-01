# ============================================================================
#  patch2d.jl — explicit 2D compressible solver for the SurfaceCode pilot.
#
#  Purpose-built (no external code available). Built INCREMENTALLY; each stage is
#  validated before the next is added (README.md §10 staged plan):
#    Stage A [this file] : 2D compressible Euler core, HLL + MUSCL(minmod) + SSP-RK2.
#                          Validated on a Sod shock tube.
#    Stage B [next]      : gravity + well-balanced stratified atmosphere (ingest
#                          background_*.txt); a hydrostatic column must stay put.
#    Stage C [next]      : MHD (fast mode) via GLM divergence cleaning.
#    Stage D [next]      : ambipolar diffusion; then the pulse-transmission measurement.
#
#  Conservative U=(ρ, ρu, ρw, E). x=lateral (ring→pole), z=height. Finite volume,
#  method of lines: dU/dt = -∂F/∂x - ∂G/∂z. Padded arrays: interior = NG+1 : NG+N.
# ============================================================================
using Printf
module Patch2D

const NG   = 2
const IRHO=1; const IMX=2; const IMZ=3; const IE=4; const NVAR=4

struct Grid
    Nx::Int; Nz::Int; dx::Float64; dz::Float64
    xc::Vector{Float64}; zc::Vector{Float64}
end
function Grid(Nx,Nz,Lx,Lz; z0=0.0)
    dx=Lx/Nx; dz=Lz/Nz
    Grid(Nx,Nz,dx,dz,[(i-0.5)*dx for i in 1:Nx],[z0+(j-0.5)*dz for j in 1:Nz])
end
alloc(g::Grid) = zeros(NVAR, g.Nx+2NG, g.Nz+2NG)

const ρFLOOR=1e-14; const pFLOOR=1e-12    # positivity floors for reconstructed face states (SI)
@inline pressure(γ,ρ,mx,mz,E) = (γ-1)*(E - 0.5*(mx^2+mz^2)/max(ρ,ρFLOOR))
@inline soundspeed(γ,ρ,p) = sqrt(γ*max(p,pFLOOR)/max(ρ,ρFLOOR))
@inline minmod(a,b) = (a*b ≤ 0) ? 0.0 : (abs(a)<abs(b) ? a : b)

"conservative → primitive (ρ,u,w,p) at padded cell (i,j)"
@inline function prim(U,γ,i,j)
    ρ=U[IRHO,i,j]; mx=U[IMX,i,j]; mz=U[IMZ,i,j]
    (ρ, mx/ρ, mz/ρ, pressure(γ,ρ,mx,mz,U[IE,i,j]))
end

"HLL flux, direction dir (1=x,2=z), from L/R primitive states"
@inline function hll(γ, ρL,uL,wL,pL, ρR,uR,wR,pR, dir)
    vnL = dir==1 ? uL : wL;  vnR = dir==1 ? uR : wR
    cL=soundspeed(γ,ρL,pL); cR=soundspeed(γ,ρR,pR)
    SL=min(vnL-cL, vnR-cR); SR=max(vnL+cL, vnR+cR)
    sf(ρ,u,w,p) = begin
        mx=ρ*u; mz=ρ*w; E=p/(γ-1)+0.5*ρ*(u^2+w^2); vn = dir==1 ? u : w
        ((ρ,mx,mz,E),
         (ρ*vn, mx*vn+(dir==1 ? p : 0.0), mz*vn+(dir==2 ? p : 0.0), (E+p)*vn))
    end
    UL,FL=sf(ρL,uL,wL,pL); UR,FR=sf(ρR,uR,wR,pR)
    SL≥0 && return FL
    SR≤0 && return FR
    ntuple(m->(SR*FL[m]-SL*FR[m]+SL*SR*(UR[m]-UL[m]))/(SR-SL), NVAR)
end

"Fill ghost cells. bcx,bcz ∈ (:outflow,:periodic,:reflect). reflect flips normal momentum."
function apply_bc!(U,g::Grid,bcx,bcz)
    Nx,Nz=g.Nx,g.Nz
    @inbounds for j in 1:Nz+2NG, m in 1:NVAR, k in 1:NG
        # x-lo ghost (NG-k+1 ..) mirror of interior; x-hi
        il=NG+1; ih=NG+Nx
        if bcx===:periodic
            U[m, il-k, j]=U[m, ih-k+1, j]; U[m, ih+k, j]=U[m, il+k-1, j]
        else # outflow / reflect share zero-gradient except reflect flips IMX
            s = (bcx===:reflect && m==IMX) ? -1.0 : 1.0
            U[m, il-k, j]=s*U[m, il+k-1, j]; U[m, ih+k, j]=s*U[m, ih-k+1, j]
        end
    end
    @inbounds for i in 1:Nx+2NG, m in 1:NVAR, k in 1:NG
        jl=NG+1; jh=NG+Nz
        if bcz===:periodic
            U[m, i, jl-k]=U[m, i, jh-k+1]; U[m, i, jh+k]=U[m, i, jl+k-1]
        else
            s = (bcz===:reflect && m==IMZ) ? -1.0 : 1.0
            U[m, i, jl-k]=s*U[m, i, jl+k-1]; U[m, i, jh+k]=s*U[m, i, jh-k+1]
        end
    end
    U
end

"dU/dt into R (method of lines). W is a scratch primitive array (NVAR × padded)."
function rhs!(R,U,W,g::Grid,γ,bcx,bcz)
    Nx,Nz=g.Nx,g.Nz
    apply_bc!(U,g,bcx,bcz)
    @inbounds for j in 1:Nz+2NG, i in 1:Nx+2NG
        ρ,u,w,p=prim(U,γ,i,j); W[1,i,j]=ρ; W[2,i,j]=u; W[3,i,j]=w; W[4,i,j]=p
    end
    fill!(R,0.0)
    # x-direction interfaces
    @inbounds for j in NG+1:NG+Nz, i in NG:NG+Nx      # interface between i and i+1
        WL=ntuple(m->W[m,i,j]  +0.5*minmod(W[m,i,j]-W[m,i-1,j], W[m,i+1,j]-W[m,i,j]),4)
        WR=ntuple(m->W[m,i+1,j]-0.5*minmod(W[m,i+1,j]-W[m,i,j], W[m,i+2,j]-W[m,i+1,j]),4)
        F=hll(γ, WL...,WR..., 1)
        for m in 1:NVAR
            R[m,i,j]   -= F[m]/g.dx
            R[m,i+1,j] += F[m]/g.dx
        end
    end
    # z-direction interfaces
    @inbounds for j in NG:NG+Nz, i in NG+1:NG+Nx
        WL=ntuple(m->W[m,i,j]  +0.5*minmod(W[m,i,j]-W[m,i,j-1], W[m,i,j+1]-W[m,i,j]),4)
        WR=ntuple(m->W[m,i,j+1]-0.5*minmod(W[m,i,j+1]-W[m,i,j], W[m,i,j+2]-W[m,i,j+1]),4)
        G=hll(γ, WL...,WR..., 2)
        for m in 1:NVAR
            R[m,i,j]   -= G[m]/g.dz
            R[m,i,j+1] += G[m]/g.dz
        end
    end
    R
end

# ==========================================================================
#  Stage B: gravity + WELL-BALANCED stratified atmosphere.
#  Background ρ0(z), p0(z) (u0=w0=0), gravity -g ẑ. In z we reconstruct the
#  DEVIATION from the background so the HLL dissipation term SL·SR·(U_R−U_L)
#  vanishes in hydrostatic balance (no spurious mass/energy flux), and the
#  gravity source is defined to cancel the background pressure-flux exactly:
#     S_mz[j] = (p0f[j+½]−p0f[j−½])/dz − (ρ−ρ0) g ,   S_E = −ρ w g.
#  x is background-uniform, so x keeps the Stage-A full-variable path.
# ==========================================================================
struct Atmos
    g::Float64
    ρ0::Vector{Float64}     # padded (length Nz+2NG), z-indexed
    p0::Vector{Float64}
end

"x-BC per bcx; z-ghosts filled with the hydrostatic background (ρ0,0,0,p0) → well-balanced walls."
function apply_bc_atmos!(U,g::Grid,a::Atmos,γ,bcx)
    Nx,Nz=g.Nx,g.Nz
    @inbounds for j in 1:Nz+2NG, m in 1:NVAR, k in 1:NG
        il=NG+1; ih=NG+Nx
        if bcx===:periodic
            U[m,il-k,j]=U[m,ih-k+1,j]; U[m,ih+k,j]=U[m,il+k-1,j]
        else
            s=(bcx===:reflect && m==IMX) ? -1.0 : 1.0
            U[m,il-k,j]=s*U[m,il+k-1,j]; U[m,ih+k,j]=s*U[m,ih-k+1,j]
        end
    end
    @inbounds for i in 1:Nx+2NG, k in 1:NG
        for jg in (NG+1-k, NG+Nz+k)
            U[IRHO,i,jg]=a.ρ0[jg]; U[IMX,i,jg]=0.0; U[IMZ,i,jg]=0.0
            U[IE,i,jg]=a.p0[jg]/(γ-1)
        end
    end
    U
end

"Well-balanced dU/dt with gravity. Optional `inject(U,g,a,γ,t)` overwrites boundary ghosts
(e.g. a time-dependent pulse) right after the standard BCs."
function rhs_grav!(R,U,W,g::Grid,γ,a::Atmos,bcx,bcz; inject=nothing, t=0.0)
    Nx,Nz=g.Nx,g.Nz
    apply_bc_atmos!(U,g,a,γ,bcx)
    inject===nothing || inject(U,g,a,γ,t)
    @inbounds for j in 1:Nz+2NG, i in 1:Nx+2NG
        ρ,u,w,p=prim(U,γ,i,j); W[1,i,j]=ρ; W[2,i,j]=u; W[3,i,j]=w; W[4,i,j]=p
    end
    fill!(R,0.0)
    # --- x-direction: background uniform in x → Stage-A full-variable path ---
    @inbounds for j in NG+1:NG+Nz, i in NG:NG+Nx
        WL=ntuple(m->W[m,i,j]  +0.5*minmod(W[m,i,j]-W[m,i-1,j], W[m,i+1,j]-W[m,i,j]),4)
        WR=ntuple(m->W[m,i+1,j]-0.5*minmod(W[m,i+1,j]-W[m,i,j], W[m,i+2,j]-W[m,i+1,j]),4)
        F=hll(γ, max(WL[1],ρFLOOR),WL[2],WL[3],max(WL[4],pFLOOR),
                 max(WR[1],ρFLOOR),WR[2],WR[3],max(WR[4],pFLOOR), 1)
        for m in 1:NVAR; R[m,i,j]-=F[m]/g.dx; R[m,i+1,j]+=F[m]/g.dx; end
    end
    # --- z-direction: DEVIATION reconstruction of ρ,p (u,w have zero bg) -----
    @inbounds for j in NG:NG+Nz, i in NG+1:NG+Nx
        ρ0f=0.5*(a.ρ0[j]+a.ρ0[j+1]); p0f=0.5*(a.p0[j]+a.p0[j+1])
        dρ(jj)=W[1,i,jj]-a.ρ0[jj]; dp(jj)=W[4,i,jj]-a.p0[jj]
        ρL=ρ0f+dρ(j)  +0.5*minmod(dρ(j)-dρ(j-1),   dρ(j+1)-dρ(j))
        ρR=ρ0f+dρ(j+1)-0.5*minmod(dρ(j+1)-dρ(j),   dρ(j+2)-dρ(j+1))
        pL=p0f+dp(j)  +0.5*minmod(dp(j)-dp(j-1),   dp(j+1)-dp(j))
        pR=p0f+dp(j+1)-0.5*minmod(dp(j+1)-dp(j),   dp(j+2)-dp(j+1))
        uL=W[2,i,j]+0.5*minmod(W[2,i,j]-W[2,i,j-1], W[2,i,j+1]-W[2,i,j])
        uR=W[2,i,j+1]-0.5*minmod(W[2,i,j+1]-W[2,i,j], W[2,i,j+2]-W[2,i,j+1])
        wL=W[3,i,j]+0.5*minmod(W[3,i,j]-W[3,i,j-1], W[3,i,j+1]-W[3,i,j])
        wR=W[3,i,j+1]-0.5*minmod(W[3,i,j+1]-W[3,i,j], W[3,i,j+2]-W[3,i,j+1])
        ρL=max(ρL,ρFLOOR); ρR=max(ρR,ρFLOOR); pL=max(pL,pFLOOR); pR=max(pR,pFLOOR)
        G=hll(γ, ρL,uL,wL,pL, ρR,uR,wR,pR, 2)
        for m in 1:NVAR; R[m,i,j]-=G[m]/g.dz; R[m,i,j+1]+=G[m]/g.dz; end
    end
    # --- well-balanced gravity source ---------------------------------------
    @inbounds for j in NG+1:NG+Nz, i in NG+1:NG+Nx
        p0fp=0.5*(a.p0[j]+a.p0[j+1]); p0fm=0.5*(a.p0[j-1]+a.p0[j])
        ρ=W[1,i,j]; w=W[3,i,j]
        R[IMZ,i,j] += (p0fp-p0fm)/g.dz - (ρ-a.ρ0[j])*a.g
        R[IE ,i,j] += -ρ*w*a.g
    end
    R
end

function step_grav!(U,R,U1,W,g::Grid,γ,a::Atmos,dt,bcx,bcz; inject=nothing, t=0.0)
    rhs_grav!(R,U,W,g,γ,a,bcx,bcz; inject=inject, t=t);    @. U1 = U + dt*R
    rhs_grav!(R,U1,W,g,γ,a,bcx,bcz; inject=inject, t=t+dt); @. U = 0.5*U + 0.5*(U1 + dt*R)
    U
end

function max_dt(U,g::Grid,γ,cfl)
    Nx,Nz=g.Nx,g.Nz; dt=Inf
    @inbounds for j in NG+1:NG+Nz, i in NG+1:NG+Nx
        ρ,u,w,p=prim(U,γ,i,j); c=soundspeed(γ,ρ,p)
        dt=min(dt, cfl*g.dx/(abs(u)+c), cfl*g.dz/(abs(w)+c))
    end
    dt
end

"SSP-RK2 step in place."
function step!(U,R,U1,W,g::Grid,γ,dt,bcx,bcz)
    rhs!(R,U,W,g,γ,bcx,bcz); @. U1 = U + dt*R
    rhs!(R,U1,W,g,γ,bcx,bcz); @. U = 0.5*U + 0.5*(U1 + dt*R)
    U
end

# --------------------------------------------------------------------------
#  Stage-A validation: Sod shock tube along x (constant in z).
# --------------------------------------------------------------------------
function sod(; Nx=200, tend=0.2, cfl=0.4, γ=1.4)
    Nz=4; g=Grid(Nx,Nz,1.0,0.02)
    U=alloc(g); R=similar(U); U1=similar(U); W=similar(U)
    for j in 1:Nz+2NG, i in 1:Nx+2NG
        x=(i-NG-0.5)*g.dx
        ρ,p = x<0.5 ? (1.0,1.0) : (0.125,0.1)
        U[IRHO,i,j]=ρ; U[IMX,i,j]=0; U[IMZ,i,j]=0; U[IE,i,j]=p/(γ-1)
    end
    t=0.0; nstep=0
    while t<tend
        dt=min(max_dt(U,g,γ,cfl), tend-t)
        step!(U,R,U1,W,g,γ,dt,:outflow,:outflow); t+=dt; nstep+=1
    end
    # sample the midline
    jm=NG+2
    xs=[g.xc[i] for i in 1:Nx]
    ρs=[U[IRHO,NG+i,jm] for i in 1:Nx]
    ps=[pressure(γ,U[IRHO,NG+i,jm],U[IMX,NG+i,jm],U[IMZ,NG+i,jm],U[IE,NG+i,jm]) for i in 1:Nx]
    us=[U[IMX,NG+i,jm]/U[IRHO,NG+i,jm] for i in 1:Nx]
    (xs,ρs,ps,us,nstep,t)
end

"Stage-B validation: isothermal atmosphere ρ0=p0=exp(-z/H), g=1, H=1, held static.
Returns (peak |w|/c_iso, max relative density drift) over nstep — both →0 if well-balanced."
function wellbalance(; Nx=8, Nz=128, Lz=6.0, γ=5/3, nstep=3000, cfl=0.4)
    z0=0.0; g=Grid(Nx,Nz,1.0,Lz; z0=z0); H=1.0; grav=1.0
    zpad(j)=z0+(j-NG-0.5)*g.dz
    ρ0=[exp(-zpad(j)/H) for j in 1:Nz+2NG]
    p0=[H*grav*exp(-zpad(j)/H) for j in 1:Nz+2NG]     # c_iso² = Hg = 1
    a=Atmos(grav,ρ0,p0)
    U=alloc(g); R=similar(U); U1=similar(U); W=similar(U)
    for j in 1:Nz+2NG, i in 1:Nx+2NG
        U[IRHO,i,j]=ρ0[j]; U[IMX,i,j]=0; U[IMZ,i,j]=0; U[IE,i,j]=p0[j]/(γ-1)
    end
    dt=max_dt(U,g,γ,cfl); ciso=sqrt(H*grav); wpeak=0.0
    for n in 1:nstep
        step_grav!(U,R,U1,W,g,γ,a,dt,:periodic,:atmos)
        @inbounds for j in NG+1:NG+Nz, i in NG+1:NG+Nx
            wpeak=max(wpeak, abs(U[IMZ,i,j]/U[IRHO,i,j]))
        end
    end
    drift=0.0
    @inbounds for j in NG+1:NG+Nz, i in NG+1:NG+Nx
        drift=max(drift, abs(U[IRHO,i,j]-ρ0[j])/ρ0[j])
    end
    (wpeak/ciso, drift, nstep)
end

end # module

# ---- run the Sod validation when executed as a script ----------------------
if abspath(PROGRAM_FILE) == @__FILE__
    xs,ρs,ps,us,nstep,t = Patch2D.sod()
    @printf("Sod shock tube: %d steps to t=%.3f, %d cells\n", nstep, t, length(xs))
    @printf("  %-7s %-10s %-10s %-10s\n","x","rho","p","u")
    for i in 1:20:length(xs)
        @printf("  %-7.3f %-10.4f %-10.4f %-10.4f\n", xs[i],ρs[i],ps[i],us[i])
    end
    # crude checks vs the exact Sod solution at t=0.2
    ρL=ρs[1]; ρR=ρs[end]
    @printf("checks: rho_L=%.3f (exact 1.000), rho_R=%.3f (exact 0.125), max|u|=%.3f (exact ~0.93)\n",
            ρL, ρR, maximum(abs,us))
    okA = abs(ρL-1.0)<0.02 && abs(ρR-0.125)<0.02 && 0.7<maximum(abs,us)<1.05
    println(okA ? "STAGE A PASS: shock capturing validated." : "STAGE A CHECK FAILED — inspect.")

    println()
    machw, drift, ns = Patch2D.wellbalance()
    @printf("Well-balancing (isothermal atmosphere, %d steps): peak |w|/c_iso=%.2e, max density drift=%.2e\n",
            ns, machw, drift)
    okB = machw < 1e-3 && drift < 1e-3
    println(okB ? "STAGE B PASS: hydrostatic atmosphere held static (well-balanced)." :
                  "STAGE B CHECK FAILED — background drifts, inspect well-balancing.")
end
