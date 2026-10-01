# ============================================================================
#  stage4_pulse.jl — SurfaceCode pilot, Stage 4c (MHD fast mode).
#
#  The actual surface fast-mode wave: real T Tauri atmosphere + uniform dipole
#  background field at the ring latitude, an impulsive compressive pulse injected
#  at the footpoint, propagated laterally, and its energy-flux decay Φ(x) measured
#  — the MHD analog of the hydro Stage-3 result. (Ambipolar damping = Stage 5.)
#
#  Usage: julia SurfaceCode/runs/stage4_pulse.jl
# ============================================================================
include(joinpath(@__DIR__,"..","src","mhd2d.jl"))
using .MHD2D
const M = MHD2D
using Printf, DelimitedFiles

const G_grav=6.674e-11; const Msun=1.989e30; const Rsun=6.957e8
const Mstar=0.5*Msun; const Rstar=2.0*Rsun
const gsurf=G_grav*Mstar/Rstar^2
const μ0=4π*1e-7; const γ=5/3
const Bstar=0.1                                   # T (1000 G) base field

read_bg(f)=(d=readdlm(f;comments=true,comment_char='#'); p=sortperm(d[:,1]); (d[p,1],d[p,3],d[p,5]))
function interp(xq,x,y)
    xq≤x[1] && return y[1]; xq≥x[end] && return y[end]
    k=searchsortedlast(x,xq); w=(xq-x[k])/(x[k+1]-x[k]); (1-w)*y[k]+w*y[k+1]
end

function build(file; Nx=500, Lx=5e7, zbot_Hp=-4.0, ztop_Hp=3.0, cph=12, θdeg=35.0)
    zb,ρb,Pb=read_bg(file)
    ρ0s=interp(0.0,zb,ρb); P0s=interp(0.0,zb,Pb); Hp=P0s/(ρ0s*gsurf)
    zlo=zbot_Hp*Hp; Lz=(ztop_Hp-zbot_Hp)*Hp
    Nz=max(64,round(Int,Lz/(Hp/cph)))
    g=M.Grid(Nx,Nz,Lx,Lz; z0=zlo)
    zpad(j)=zlo+(j-M.NG-0.5)*g.dz
    ρ0=[max(interp(zpad(j),zb,ρb),1e-12) for j in 1:Nz+2M.NG]
    p0=[max(interp(zpad(j),zb,Pb),1e-6) for j in 1:Nz+2M.NG]
    θ=deg2rad(θdeg); Bz0=Bstar*cos(θ); Bx0=0.5*Bstar*sin(θ)     # dipole at r=R⋆
    a=M.Atmos(gsurf,ρ0,p0,Bx0,Bz0)
    cf=sqrt(γ*P0s/ρ0s + (Bx0^2+Bz0^2)/(μ0*ρ0s))                # fast speed at surface
    (g,a,Hp,cf,ρ0s)
end

const TAU_S = haskey(ENV,"TAU_S") ? parse(Float64,ENV["TAU_S"]) : 150.0   # driver duration [s]

function run_pulse(g,a,cf,Hp; M0=2.0, τ=TAU_S, cfl=0.4, textra=1.3)
    U=M.alloc(g); R=similar(U); U1=similar(U); W=similar(U); Fx=similar(U); Fz=similar(U)
    @inbounds for j in 1:g.Nz+2M.NG, i in 1:g.Nx+2M.NG
        c=M.cons(γ,μ0, a.ρ0[j],0.0,0.0,0.0, a.Bx0[j],0.0,a.Bz0, a.p0[j],0.0)
        for m in 1:M.NVAR; U[m,i,j]=c[m]; end
    end
    zpad(j)=g.zc[1]-M.NG*g.dz+(j-1)*g.dz
    zw(j)=exp(-(zpad(j)/(1.5*Hp))^2)
    function inject!(U,gg,aa,γ,μ0,t)
        (0.0≤t≤τ) || return
        vx=M0*cf*sin(π*t/τ)
        @inbounds for j in 1:gg.Nz+2M.NG, k in 1:M.NG
            il=M.NG+1-k; ρ=aa.ρ0[j]; ux=vx*zw(j)
            c=M.cons(γ,μ0, ρ,ux,0.0,0.0, aa.Bx0[j],0.0,aa.Bz0, aa.p0[j],0.0)
            for m in 1:M.NVAR; U[m,il,j]=c[m]; end
        end
    end
    tend=textra*(g.Nx*g.dx)/cf+τ
    fluxx=zeros(g.Nx); t=0.0; nstep=0
    while t<tend
        ch=M.maxspeed(U,g,γ,μ0); dt=cfl*min(g.dx,g.dz)/ch
        M.step_grav!(U,R,U1,W,Fx,Fz,g,γ,μ0,dt,ch,a,:outflow; inject=inject!,t=t); t+=dt; nstep+=1
        if nstep % 1000 == 0                                # live progress + crash guard
            vmax=0.0
            @inbounds for j in M.NG+1:M.NG+g.Nz, i in M.NG+1:M.NG+g.Nx
                vmax=max(vmax, abs(U[M.IMX,i,j]/max(U[M.IRHO,i,j],1e-14)))
            end
            @printf("      step %6d  t=%6.0f/%.0f s (%3.0f%%)  dt=%.2fs  max|vx|=%.2f km/s\n",
                    nstep, t, tend, 100t/tend, dt, vmax/1e3); flush(stdout)
            isfinite(vmax) || (println("      ** NaN/Inf — aborting this run **"); flush(stdout); break)
        end
        # Threaded over i only: each fluxx[i] is touched by exactly one thread and the
        # j-summation order is unchanged, so this is bit-identical to the serial form.
        Threads.@threads for i in 1:g.Nx
            fx=0.0
            @inbounds for j in M.NG+1:M.NG+g.Nz
                ρ,vx,vy,vz,Bx,By,Bz,p,ψ=M.prim(γ,μ0,U,M.NG+i,j)
                pB=0.5*(Bx^2+By^2+Bz^2)/μ0; pt=p+pB; vB=vx*Bx+vy*By+vz*Bz
                E=p/(γ-1)+0.5*ρ*(vx^2+vy^2+vz^2)+pB
                fx += (E+pt)*vx - Bx*vB/μ0          # MHD x energy flux (0 in the unperturbed atm.)
            end
            @inbounds fluxx[i]+= fx*g.dz*dt
        end
    end
    (collect(g.xc), fluxx, nstep)
end

function decay_length(x,f; i0f=0.15,i1f=0.9)
    n=length(x); i0=max(2,round(Int,i0f*n)); i1=round(Int,i1f*n)
    xs=x[i0:i1]; ys=log.(max.(abs.(f[i0:i1]),1e-30))
    β=hcat(ones(length(xs)),xs)\ys
    β[2]<0 ? -1/β[2] : Inf
end

if abspath(PROGRAM_FILE)==@__FILE__
    bg=joinpath(@__DIR__,"..","data","background_35deg.txt"); Lx=0.036*Rstar
    g,a,Hp,cf,ρ0s=build(bg; Nx=500, Lx=Lx)
    VA=sqrt((a.Bx0[M.NG+1]^2+a.Bz0^2)/(μ0*ρ0s))
    @printf("MHD atmosphere: Nx=%d Nz=%d, Lx=%.2e m (%.3f Rstar), Hp=%.0f km\n", g.Nx,g.Nz,Lx,Lx/Rstar,Hp/1e3)
    @printf("  surface: c_s=%.1f, V_A=%.1f, c_fast=%.1f km/s, |B0|=%.0f G (Bx0=%.0f, Bz0=%.0f G)\n",
            sqrt(γ*a.p0[M.NG+1]/ρ0s)/1e3, VA/1e3, cf/1e3, sqrt(a.Bx0[M.NG+1]^2+a.Bz0^2)*1e4, a.Bx0[M.NG+1]*1e4, a.Bz0*1e4)
    println("\nfast-mode energy-flux decay length L(Φ)  [convergence: L should be grid-independent]:")
    reslist = isempty(ARGS) ? (1,2,4) : Tuple(parse(Int,a) for a in ARGS)   # e.g. `... stage4_pulse.jl 8`
    for res in reslist
        g2,a2,Hp2,cf2,_=build(bg; Nx=500*res, Lx=Lx, cph=12*res)
        @printf("  --- res×%d: Nx=%d Nz=%d dx=%.0fkm ---\n", res,g2.Nx,g2.Nz,g2.dx/1e3); flush(stdout)
        x,flux,ns=run_pulse(g2,a2,cf2,Hp2; M0=2.0)
        L=decay_length(x,flux)
        @printf("  ==> res×%d: L(Φ)=%.4f Rstar (%.2e m) [%d steps]\n\n",
                res,L/Rstar,L,ns); flush(stdout)
    end
    @printf("\ntransit ring→pole ≈ 0.61 Rstar. If L(Φ) ≪ that (and grid-converged), the fast-mode\n")
    @printf("surface wave does not survive to the pole — the resolved test of Cranmer's channel.\n")
end
