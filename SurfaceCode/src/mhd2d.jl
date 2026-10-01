# ============================================================================
#  mhd2d.jl — explicit 2.5D ideal-MHD solver for the SurfaceCode pilot, Stage 4.
#
#  Extends the validated hydro core (patch2d.jl) to MHD so the actual fast-mode
#  surface wave can be propagated (hydro Stage 3 only had the acoustic proxy).
#  State: U=(ρ, ρvx, ρvy, ρvz, Bx, By, Bz, E, ψ) — 2D domain (x,z), all 3 vector
#  components (2.5D). Dedner GLM divergence cleaning (ψ). HLL for the 7 genuine
#  MHD waves + the exact 2×2 (Bn,ψ) GLM Riemann flux. MUSCL(minmod) + SSP-RK2.
#
#  Built incrementally, validated each stage:
#    4a [this file]: MHD core, validated on the Brio–Wu shock tube (γ=2).
#    4b [next]: gravity + well-balanced magneto-atmosphere (uniform background B).
#    4c [next]: fast-mode pulse + Φ(x) transmission — the actual Stage-4 result.
# ============================================================================
using Printf
module MHD2D
using Base.Threads: @threads, nthreads, threadid, maxthreadid

const NG=2
const IRHO=1; const IMX=2; const IMY=3; const IMZ=4; const IBX=5; const IBY=6; const IBZ=7; const IE=8; const IPSI=9
const NVAR=9
const ρFLOOR=1e-14; const pFLOOR=1e-12

struct Grid
    Nx::Int; Nz::Int; dx::Float64; dz::Float64; xc::Vector{Float64}; zc::Vector{Float64}
end
Grid(Nx,Nz,Lx,Lz; z0=0.0) = Grid(Nx,Nz,Lx/Nx,Lz/Nz,
    [(i-0.5)*(Lx/Nx) for i in 1:Nx],[z0+(j-0.5)*(Lz/Nz) for j in 1:Nz])
alloc(g::Grid)=zeros(NVAR, g.Nx+2NG, g.Nz+2NG)

@inline minmod(a,b)=(a*b≤0) ? 0.0 : (abs(a)<abs(b) ? a : b)

"pressure from conservative U at (i,j)"
@inline function gaspressure(γ,μ0,U,i,j)
    ρ=max(U[IRHO,i,j],ρFLOOR)
    ke=0.5*(U[IMX,i,j]^2+U[IMY,i,j]^2+U[IMZ,i,j]^2)/ρ
    me=0.5*(U[IBX,i,j]^2+U[IBY,i,j]^2+U[IBZ,i,j]^2)/μ0
    max((γ-1)*(U[IE,i,j]-ke-me), pFLOOR)
end

"primitive tuple (ρ,vx,vy,vz,Bx,By,Bz,p,ψ) at (i,j)"
@inline function prim(γ,μ0,U,i,j)
    ρ=max(U[IRHO,i,j],ρFLOOR)
    (ρ, U[IMX,i,j]/ρ, U[IMY,i,j]/ρ, U[IMZ,i,j]/ρ,
     U[IBX,i,j], U[IBY,i,j], U[IBZ,i,j], gaspressure(γ,μ0,U,i,j), U[IPSI,i,j])
end

"conservative from primitive (ρ,vx,vy,vz,Bx,By,Bz,p,ψ)"
@inline function cons(γ,μ0,ρ,vx,vy,vz,Bx,By,Bz,p,ψ)
    E=p/(γ-1)+0.5*ρ*(vx^2+vy^2+vz^2)+0.5*(Bx^2+By^2+Bz^2)/μ0
    (ρ, ρ*vx, ρ*vy, ρ*vz, Bx, By, Bz, E, ψ)
end

"fast magnetosonic speed in the normal direction (Bn = normal field component)"
@inline function cfast(γ,μ0,ρ,p,Bx,By,Bz,Bn)
    a2=γ*max(p,pFLOOR)/ρ; b2=(Bx^2+By^2+Bz^2)/(μ0*ρ); bn2=Bn^2/(μ0*ρ)
    sqrt(max(0.5*(a2+b2+sqrt(max((a2+b2)^2-4*a2*bn2,0.0))), 0.0))
end

"ideal-MHD physical flux (dir=1→x, 2→z) from a primitive state; GLM ch for Bn,ψ."
@inline function physflux(γ,μ0,ch, ρ,vx,vy,vz,Bx,By,Bz,p,ψ, dir)
    pB=0.5*(Bx^2+By^2+Bz^2)/μ0; pt=p+pB; vB=vx*Bx+vy*By+vz*Bz
    E=p/(γ-1)+0.5*ρ*(vx^2+vy^2+vz^2)+pB
    if dir==1
        vn=vx; Bn=Bx
        (ρ*vn, ρ*vn*vx+pt-Bn*Bx/μ0, ρ*vn*vy-Bn*By/μ0, ρ*vn*vz-Bn*Bz/μ0,
         ψ, vn*By-vy*Bn, vn*Bz-vz*Bn, (E+pt)*vn-Bn*vB/μ0, ch^2*Bn)
    else
        vn=vz; Bn=Bz
        (ρ*vn, ρ*vn*vx-Bn*Bx/μ0, ρ*vn*vy-Bn*By/μ0, ρ*vn*vz+pt-Bn*Bz/μ0,
         vn*Bx-vx*Bn, vn*By-vy*Bn, ψ, (E+pt)*vn-Bn*vB/μ0, ch^2*Bn)
    end
end

"HLL flux for the 7 MHD waves + exact 2×2 GLM (Bn,ψ), written DIRECTLY into Fout[:,i,j]
(no returned tuple, no closures — allocation-free hot path)."
@inline function mhdflux!(Fout, i, j, γ,μ0,ch, WL,WR, dir)
    ρL,vxL,vyL,vzL,BxL,ByL,BzL,pL,ψL=WL
    ρR,vxR,vyR,vzR,BxR,ByR,BzR,pR,ψR=WR
    vnL = dir==1 ? vxL : vzL; vnR = dir==1 ? vxR : vzR
    BnL = dir==1 ? BxL : BzL; BnR = dir==1 ? BxR : BzR
    cfL=cfast(γ,μ0,ρL,pL,BxL,ByL,BzL,BnL); cfR=cfast(γ,μ0,ρR,pR,BxR,ByR,BzR,BnR)
    SL=min(vnL-cfL, vnR-cfR); SR=max(vnL+cfL, vnR+cfR)
    FL=physflux(γ,μ0,ch, WL..., dir); FR=physflux(γ,μ0,ch, WR..., dir)
    UL=cons(γ,μ0, WL...); UR=cons(γ,μ0, WR...)
    inv=1.0/(SR-SL)
    @inbounds for m in 1:NVAR
        Fout[m,i,j] = SL≥0 ? FL[m] : SR≤0 ? FR[m] : (SR*FL[m]-SL*FR[m]+SL*SR*(UR[m]-UL[m]))*inv
    end
    # exact GLM Riemann flux for the (Bn,ψ) subsystem (linear, speed ch)
    @inbounds begin
        Fout[dir==1 ? IBX : IBZ, i, j] = 0.5*(ψL+ψR) - 0.5*ch*(BnR-BnL)
        Fout[IPSI, i, j]               = 0.5*ch^2*(BnL+BnR) - 0.5*ch*(ψR-ψL)
    end
    nothing
end

function apply_bc!(U,g::Grid,bcx,bcz)
    Nx,Nz=g.Nx,g.Nz
    @inbounds for j in 1:Nz+2NG, m in 1:NVAR, k in 1:NG
        il=NG+1; ih=NG+Nx
        if bcx===:periodic
            U[m,il-k,j]=U[m,ih-k+1,j]; U[m,ih+k,j]=U[m,il+k-1,j]
        else
            U[m,il-k,j]=U[m,il+k-1,j]; U[m,ih+k,j]=U[m,ih-k+1,j]
        end
    end
    @inbounds for i in 1:Nx+2NG, m in 1:NVAR, k in 1:NG
        jl=NG+1; jh=NG+Nz
        if bcz===:periodic
            U[m,i,jl-k]=U[m,i,jh-k+1]; U[m,i,jh+k]=U[m,i,jl+k-1]
        else
            U[m,i,jl-k]=U[m,i,jl+k-1]; U[m,i,jh+k]=U[m,i,jh-k+1]
        end
    end
    U
end

"max fast-magnetosonic wave speed (for CFL and GLM ch) — threaded reduction over columns."
function maxspeed(U,g::Grid,γ,μ0)
    part=zeros(maxthreadid())        # threadid() may exceed nthreads() (interactive pool) in Julia ≥1.9
    @threads for j in NG+1:NG+g.Nz
        s=0.0
        @inbounds for i in NG+1:NG+g.Nx
            ρ,vx,vy,vz,Bx,By,Bz,p,ψ=prim(γ,μ0,U,i,j)
            cfx=cfast(γ,μ0,ρ,p,Bx,By,Bz,Bx); cfz=cfast(γ,μ0,ρ,p,Bx,By,Bz,Bz)
            s=max(s, abs(vx)+cfx, abs(vz)+cfz)
        end
        part[threadid()]=max(part[threadid()],s)
    end
    maximum(part)
end

@inline recL(W,m,i,j,di,dj)=W[m,i,j]+0.5*minmod(W[m,i,j]-W[m,i-di,j-dj], W[m,i+di,j+dj]-W[m,i,j])
@inline recRR(W,m,i,j,di,dj)=W[m,i,j]-0.5*minmod(W[m,i,j]-W[m,i-di,j-dj], W[m,i+di,j+dj]-W[m,i,j])
# closure-free full-state reconstruction (explicit 9-tuple → no boxed loop-var capture)
@inline recLall(W,i,j,di,dj)=(recL(W,1,i,j,di,dj),recL(W,2,i,j,di,dj),recL(W,3,i,j,di,dj),
    recL(W,4,i,j,di,dj),recL(W,5,i,j,di,dj),recL(W,6,i,j,di,dj),recL(W,7,i,j,di,dj),
    recL(W,8,i,j,di,dj),recL(W,9,i,j,di,dj))
@inline recRall(W,i,j,di,dj)=(recRR(W,1,i,j,di,dj),recRR(W,2,i,j,di,dj),recRR(W,3,i,j,di,dj),
    recRR(W,4,i,j,di,dj),recRR(W,5,i,j,di,dj),recRR(W,6,i,j,di,dj),recRR(W,7,i,j,di,dj),
    recRR(W,8,i,j,di,dj),recRR(W,9,i,j,di,dj))

"Race-free threaded RHS: fill primitives, compute x/z interface flux arrays, then differences —
each stage threaded over j-columns (no scatter, barrier between stages)."
function rhs!(R,U,W,Fx,Fz,g::Grid,γ,μ0,ch,bcx,bcz)
    Nx,Nz=g.Nx,g.Nz
    apply_bc!(U,g,bcx,bcz)
    @threads for j in 1:Nz+2NG
        @inbounds for i in 1:Nx+2NG
            pr=prim(γ,μ0,U,i,j); for m in 1:NVAR; W[m,i,j]=pr[m]; end
        end
    end
    @threads for j in NG+1:NG+Nz                            # x-interfaces (Fx[·,i,j]=flux at i+½)
        @inbounds for i in NG:NG+Nx
            mhdflux!(Fx,i,j, γ,μ0,ch, recLall(W,i,j,1,0), recRall(W,i+1,j,1,0), 1)
        end
    end
    @threads for j in NG:NG+Nz                              # z-interfaces (Fz[·,i,j]=flux at j+½)
        @inbounds for i in NG+1:NG+Nx
            mhdflux!(Fz,i,j, γ,μ0,ch, recLall(W,i,j,0,1), recRall(W,i,j+1,0,1), 2)
        end
    end
    @threads for j in NG+1:NG+Nz                            # flux differences (no races)
        @inbounds for i in NG+1:NG+Nx, m in 1:NVAR
            R[m,i,j]=-(Fx[m,i,j]-Fx[m,i-1,j])/g.dx - (Fz[m,i,j]-Fz[m,i,j-1])/g.dz
        end
    end
    R
end

"SSP-RK2 step with GLM parabolic ψ-damping folded in."
function step!(U,R,U1,W,Fx,Fz,g::Grid,γ,μ0,dt,ch,bcx,bcz)
    rhs!(R,U,W,Fx,Fz,g,γ,μ0,ch,bcx,bcz);  @. U1=U+dt*R
    rhs!(R,U1,W,Fx,Fz,g,γ,μ0,ch,bcx,bcz); @. U=0.5*U+0.5*(U1+dt*R)
    α=exp(-0.18*ch*dt/min(g.dx,g.dz))
    @threads for j in NG+1:NG+g.Nz
        @inbounds for i in NG+1:NG+g.Nx; U[IPSI,i,j]*=α; end
    end
    U
end

# ==========================================================================
#  Stage 4b: gravity + WELL-BALANCED magneto-atmosphere.
#  Background field B0=(Bx0(z),0,Bz0). Two supported configurations:
#    * uniform Bx0 (science runs): no background Lorentz force, gas hydrostatic
#      balance unchanged;
#    * z-dependent Bx0 with Bz0=0 (verification test 1): magneto-hydrostatic
#      balance d/dz(p0 + Bx0²/2μ0) = -ρ0 g. (Bz0≠0 with varying Bx0 would leave
#      an unbalanced background tension force Bz0·dBx0/dz in x — rejected.)
#  Deviation reconstruction of ρ,p,Bx in z (as the hydro Stage B) + matched
#  gravity/magnetic-pressure source; v,By,Bz,ψ reconstructed normally.
# ==========================================================================
struct Atmos
    grav::Float64
    ρ0::Vector{Float64}; p0::Vector{Float64}       # padded z-arrays
    Bx0::Vector{Float64}                            # padded z-array (horizontal background field)
    Bz0::Float64                                    # uniform vertical background field
    function Atmos(grav,ρ0,p0,Bx0::AbstractVector,Bz0)
        (Bz0==0.0 || all(==(Bx0[1]),Bx0)) ||
            error("z-dependent Bx0 requires Bz0=0 (background J×B has an unbalanced x-component otherwise)")
        new(grav,ρ0,p0,Bx0,Bz0)
    end
end
"uniform-Bx0 convenience constructor (the science-run configuration)"
Atmos(grav,ρ0,p0,Bx0::Real,Bz0::Real)=Atmos(grav,ρ0,p0,fill(float(Bx0),length(ρ0)),float(Bz0))

"z-ghosts = hydrostatic background (uniform B0, v=0); x-ghosts per bcx."
function apply_bc_atmos!(U,g::Grid,a::Atmos,γ,μ0,bcx; bczlo::Symbol=:fixed)
    Nx,Nz=g.Nx,g.Nz
    @inbounds for j in 1:Nz+2NG, m in 1:NVAR, k in 1:NG
        il=NG+1; ih=NG+Nx
        if bcx===:periodic
            U[m,il-k,j]=U[m,ih-k+1,j]; U[m,ih+k,j]=U[m,il+k-1,j]
        else
            U[m,il-k,j]=U[m,il+k-1,j]; U[m,ih+k,j]=U[m,ih-k+1,j]
        end
    end
    @inbounds for i in 1:Nx+2NG, k in 1:NG
        # top ghost: fixed background (unchanged)
        jg=NG+Nz+k
        c=cons(γ,μ0, a.ρ0[jg],0.0,0.0,0.0, a.Bx0[jg],0.0,a.Bz0, a.p0[jg],0.0)
        for m in 1:NVAR; U[m,i,jg]=c[m]; end
        # bottom ghost: :fixed (default, unchanged) or :reflect — Krause et al. 2015 Sect. 2.1
        # use a reflecting lower boundary "to model the denser solar surface values"; a fixed
        # background floor instead ABSORBS the downward compression the coronal shock drives
        # into the slab, which is precisely the chromospheric signal being measured.
        jg=NG+1-k
        if bczlo===:reflect
            js=NG+k                                   # mirror partner inside the domain
            for m in 1:NVAR; U[m,i,jg]=U[m,i,js]; end
            U[IMZ,i,jg]=-U[IMZ,i,js]                  # flip the normal momentum
        else
            c=cons(γ,μ0, a.ρ0[jg],0.0,0.0,0.0, a.Bx0[jg],0.0,a.Bz0, a.p0[jg],0.0)
            for m in 1:NVAR; U[m,i,jg]=c[m]; end
        end
    end
    U
end

function rhs_grav!(R,U,W,Fx,Fz,g::Grid,γ,μ0,ch,a::Atmos,bcx; inject=nothing, t=0.0, bczlo::Symbol=:fixed)
    Nx,Nz=g.Nx,g.Nz
    apply_bc_atmos!(U,g,a,γ,μ0,bcx; bczlo=bczlo)
    inject===nothing || inject(U,g,a,γ,μ0,t)
    @threads for j in 1:Nz+2NG
        @inbounds for i in 1:Nx+2NG
            pr=prim(γ,μ0,U,i,j); for m in 1:NVAR; W[m,i,j]=pr[m]; end
        end
    end
    @threads for j in NG+1:NG+Nz                            # x-interfaces (bg uniform in x)
        @inbounds for i in NG:NG+Nx
            mhdflux!(Fx,i,j, γ,μ0,ch, recLall(W,i,j,1,0), recRall(W,i+1,j,1,0), 1)
        end
    end
    @threads for j in NG:NG+Nz                              # z-interfaces (DEVIATION recon of ρ,p,Bx)
        @inbounds for i in NG+1:NG+Nx
            ρ0f=0.5*(a.ρ0[j]+a.ρ0[j+1]); p0f=0.5*(a.p0[j]+a.p0[j+1]); Bx0f=0.5*(a.Bx0[j]+a.Bx0[j+1])
            dρj=W[1,i,j]-a.ρ0[j]; dρj1=W[1,i,j+1]-a.ρ0[j+1]
            dρm=W[1,i,j-1]-a.ρ0[j-1]; dρp=W[1,i,j+2]-a.ρ0[j+2]
            dpj=W[8,i,j]-a.p0[j]; dpj1=W[8,i,j+1]-a.p0[j+1]
            dpm=W[8,i,j-1]-a.p0[j-1]; dpp=W[8,i,j+2]-a.p0[j+2]
            dBj=W[5,i,j]-a.Bx0[j]; dBj1=W[5,i,j+1]-a.Bx0[j+1]
            dBm=W[5,i,j-1]-a.Bx0[j-1]; dBp=W[5,i,j+2]-a.Bx0[j+2]
            ρL=max(ρ0f+dρj +0.5*minmod(dρj-dρm, dρj1-dρj), ρFLOOR)
            pL=max(p0f+dpj +0.5*minmod(dpj-dpm, dpj1-dpj), pFLOOR)
            ρR=max(ρ0f+dρj1-0.5*minmod(dρj1-dρj, dρp-dρj1), ρFLOOR)
            pR=max(p0f+dpj1-0.5*minmod(dpj1-dpj, dpp-dpj1), pFLOOR)
            BxL=Bx0f+dBj +0.5*minmod(dBj-dBm, dBj1-dBj)
            BxR=Bx0f+dBj1-0.5*minmod(dBj1-dBj, dBp-dBj1)
            WL=(ρL, recL(W,2,i,j,0,1),recL(W,3,i,j,0,1),recL(W,4,i,j,0,1), BxL,
                recL(W,6,i,j,0,1),recL(W,7,i,j,0,1), pL, recL(W,9,i,j,0,1))
            WR=(ρR, recRR(W,2,i,j+1,0,1),recRR(W,3,i,j+1,0,1),recRR(W,4,i,j+1,0,1), BxR,
                recRR(W,6,i,j+1,0,1),recRR(W,7,i,j+1,0,1), pR, recRR(W,9,i,j+1,0,1))
            mhdflux!(Fz,i,j, γ,μ0,ch, WL,WR, 2)
        end
    end
    @threads for j in NG+1:NG+Nz                            # differences + well-balanced gravity source
        @inbounds for i in NG+1:NG+Nx
            for m in 1:NVAR
                R[m,i,j]=-(Fx[m,i,j]-Fx[m,i-1,j])/g.dx - (Fz[m,i,j]-Fz[m,i,j-1])/g.dz
            end
            p0fp=0.5*(a.p0[j]+a.p0[j+1]); p0fm=0.5*(a.p0[j-1]+a.p0[j])
            pB0fp=(0.5*(a.Bx0[j]+a.Bx0[j+1]))^2/(2μ0); pB0fm=(0.5*(a.Bx0[j-1]+a.Bx0[j]))^2/(2μ0)
            ρ=W[1,i,j]; vz=W[4,i,j]
            R[IMZ,i,j] += (p0fp-p0fm+pB0fp-pB0fm)/g.dz - (ρ-a.ρ0[j])*a.grav
            R[IE ,i,j] += -ρ*vz*a.grav
        end
    end
    R
end

function step_grav!(U,R,U1,W,Fx,Fz,g::Grid,γ,μ0,dt,ch,a::Atmos,bcx; inject=nothing,t=0.0,bczlo::Symbol=:fixed)
    rhs_grav!(R,U,W,Fx,Fz,g,γ,μ0,ch,a,bcx; inject=inject,t=t,bczlo=bczlo);      @. U1=U+dt*R
    rhs_grav!(R,U1,W,Fx,Fz,g,γ,μ0,ch,a,bcx; inject=inject,t=t+dt,bczlo=bczlo);  @. U=0.5*U+0.5*(U1+dt*R)
    α=exp(-0.18*ch*dt/min(g.dx,g.dz))
    @threads for j in NG+1:NG+g.Nz
        @inbounds for i in NG+1:NG+g.Nx; U[IPSI,i,j]*=α; end
    end
    U
end

"Stage-4b validation: isothermal magneto-atmosphere (uniform B0) held static."
function wellbalance(; Nx=8, Nz=128, Lz=6.0, γ=5/3, μ0=1.0, nstep=2000, cfl=0.4, β=1.0)
    z0=0.0; g=Grid(Nx,Nz,1.0,Lz; z0=z0); H=1.0; grav=1.0
    zpad(j)=z0+(j-NG-0.5)*g.dz
    ρ0=[exp(-zpad(j)/H) for j in 1:Nz+2NG]; p0=[H*grav*exp(-zpad(j)/H) for j in 1:Nz+2NG]
    Bz0=sqrt(2*μ0*p0[NG+1]/β); Bx0=0.3*Bz0                 # uniform field, plasma β at base
    a=Atmos(grav,ρ0,p0,Bx0,Bz0)
    U=alloc(g)
    for j in 1:Nz+2NG, i in 1:Nx+2NG
        c=cons(γ,μ0, ρ0[j],0.0,0.0,0.0, Bx0,0.0,Bz0, p0[j],0.0); for m in 1:NVAR; U[m,i,j]=c[m]; end
    end
    R=similar(U);U1=similar(U);W=similar(U);Fx=similar(U);Fz=similar(U)
    ch=maxspeed(U,g,γ,μ0); dt=cfl*min(g.dx,g.dz)/ch; ciso=sqrt(H*grav); vpeak=0.0
    for n in 1:nstep
        ch=maxspeed(U,g,γ,μ0); dt=cfl*min(g.dx,g.dz)/ch
        step_grav!(U,R,U1,W,Fx,Fz,g,γ,μ0,dt,ch,a,:periodic)
        @inbounds for j in NG+1:NG+Nz, i in NG+1:NG+Nx
            vpeak=max(vpeak, abs(U[IMZ,i,j]/U[IRHO,i,j]), abs(U[IMX,i,j]/U[IRHO,i,j]))
        end
    end
    drift=0.0
    @inbounds for j in NG+1:NG+Nz, i in NG+1:NG+Nx; drift=max(drift, abs(U[IRHO,i,j]-ρ0[j])/ρ0[j]); end
    (vpeak/ciso, drift, nstep)
end

# --- Stage-4a validation: Brio–Wu MHD shock tube (γ=2, μ0=1) ------------------
function briowu(; Nx=400, tend=0.1, cfl=0.4)
    γ=2.0; μ0=1.0; Nz=4; g=Grid(Nx,Nz,1.0,0.01)
    U=alloc(g)
    for j in 1:Nz+2NG, i in 1:Nx+2NG
        x=(i-NG-0.5)*g.dx
        ρ,p,By = x<0.5 ? (1.0,1.0,1.0) : (0.125,0.1,-1.0)
        c=cons(γ,μ0, ρ,0.0,0.0,0.0, 0.75,By,0.0, p,0.0)
        for m in 1:NVAR; U[m,i,j]=c[m]; end
    end
    R=similar(U);U1=similar(U);W=similar(U);Fx=similar(U);Fz=similar(U); t=0.0;n=0
    while t<tend
        ch=maxspeed(U,g,γ,μ0); dt=min(cfl*min(g.dx,g.dz)/ch, tend-t)
        step!(U,R,U1,W,Fx,Fz,g,γ,μ0,dt,ch,:outflow,:outflow); t+=dt;n+=1
    end
    jm=NG+2; xs=[g.xc[i] for i in 1:Nx]
    ρs=[U[IRHO,NG+i,jm] for i in 1:Nx]
    Bys=[U[IBY,NG+i,jm] for i in 1:Nx]
    (xs,ρs,Bys,n,t)
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    xs,ρs,Bys,n,t = MHD2D.briowu()
    @printf("Brio–Wu shock tube: %d steps to t=%.3f, %d cells\n", n, t, length(xs))
    @printf("  %-7s %-9s %-9s\n","x","rho","By")
    for i in 1:40:length(xs); @printf("  %-7.3f %-9.4f %-9.4f\n", xs[i],ρs[i],Bys[i]); end
    # Brio–Wu reference at t=0.1: ρ_L≈1.0, ρ_R≈0.125, By flips sign; peak ρ ~0.5 plateau.
    ρL=ρs[1]; ρR=ρs[end]; ByL=Bys[1]; ByR=Bys[end]
    @printf("checks: rho_L=%.3f (~1.0), rho_R=%.3f (~0.125), By_L=%.2f (~1.0), By_R=%.2f (~-1.0)\n",
            ρL,ρR,ByL,ByR)
    ok = abs(ρL-1.0)<0.03 && abs(ρR-0.125)<0.03 && ByL>0.8 && ByR<-0.8 && all(0.1 .< ρs .< 1.05)
    println(ok ? "STAGE 4a PASS: MHD shock capturing validated." : "STAGE 4a CHECK — inspect.")

    println()
    for β in (1.0, 0.1)
        local machvb, driftb, nsb = MHD2D.wellbalance(β=β)
        @printf("Magneto-atmosphere well-balancing (β=%.1f, %d steps): peak|v|/c_iso=%.2e, drift=%.2e\n",
                β, nsb, machvb, driftb)
    end
    machv,drift,_ = MHD2D.wellbalance(β=0.1)
    println(machv<1e-3 && drift<1e-3 ? "STAGE 4b PASS: magneto-atmosphere held static (well-balanced)." :
                                        "STAGE 4b CHECK — inspect.")
end
