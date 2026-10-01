# profile_step.jl — where does a production step actually spend its time?
#
# MEASURED 2026-07-25 (14-core Ryzen, Nx=2783 Nz=168, 4 threads):
#   * rhs_grav! = 75% of the step; the two serial RK broadcasts = 17.5%.
#   * Threading the broadcasts: NO GAIN (0.9x) — they already run at ~19 GB/s, i.e. at DRAM speed.
#   * Thread scaling SATURATES AT ~4 THREADS: 90.4 (2t) / 66.6 (4t) / 64.7 (8t) / 65.7 (14t) /
#     68.0 ms (16t). Same at the res8 shape (Nx=4000 Nz=672): 379 (4t) / 342 (8t) / 375 ms (14t).
#     => running with -t 14 buys nothing over -t 4..8. Use -t 4 and run jobs concurrently.
#   * 3 concurrent jobs @ 4 threads: 167 ms/step each => 1.2x AGGREGATE throughput only.
#     The box is memory-bandwidth saturated.
#   * In-cache vs DRAM per-cell cost: 46.5 ns (Nx<=256) vs 71.7 ns (Nx=2783) => ~65% compute,
#     ~35% memory. Any cache-blocking / flux-fusion rewrite is therefore capped at 1.54x.
#   * NULL RESULT: hoisting the j-invariant background terms out of the i-loops (z-flux and
#     difference loops), hoisting the loop-invariant HLL branch out of the m-loop, and replacing
#     /g.dx by *idx gave EXACTLY 1.00x — LLVM already performs all three. Do not redo these.
#   * Remaining real lever = AoS->SoA restructuring so the i-loop vectorizes (state is
#     U[NVAR,i,j], all flux math is scalar per interface). Potentially 2-4x, but a major
#     rewrite of a validated solver.
# Decomposes step_grav! into its phases at the real production grid, plus the driver-level
# diagnostic (the Phi(x) flux accumulation in paper_run.jl, which runs every step).
# Usage: julia -t 14 --project=. SurfaceCode/runs/profile_step.jl [Nx] [Nz]
include(joinpath(@__DIR__,"..","src","mhd2d.jl")); using .MHD2D; const M=MHD2D
using Printf, DelimitedFiles
using Base.Threads: @threads
const Gc=6.674e-11; const Msun=1.989e30; const Rsun=6.957e8
const Mstar=0.5*Msun; const Rstar=2.0*Rsun; const gsurf=Gc*Mstar/Rstar^2
const μ0c=4π*1e-7; const γc=5/3; const Bstar=0.1

read_bg(f)=(d=readdlm(f;comments=true,comment_char='#'); p=sortperm(d[:,1]); (d[p,1],d[p,3],d[p,5]))
function interp(xq,x,y); xq≤x[1] && return y[1]; xq≥x[end] && return y[end]
    k=searchsortedlast(x,xq); w=(xq-x[k])/(x[k+1]-x[k]); (1-w)*y[k]+w*y[k+1]; end

Nx = length(ARGS)≥1 ? parse(Int,ARGS[1]) : 2783
Nz = length(ARGS)≥2 ? parse(Int,ARGS[2]) : 168
γ=γc; μ0=μ0c
zb,ρb,Pb=read_bg(joinpath(@__DIR__,"..","data","background_35deg.txt"))
ρ0s=interp(0.,zb,ρb); P0s=interp(0.,zb,Pb); Hp=P0s/(ρ0s*gsurf)
Lx=0.1*Rstar; zlo=-4Hp; Lz=7Hp
g=M.Grid(Nx,Nz,Lx,Lz; z0=zlo); zpad(j)=zlo+(j-M.NG-0.5)*g.dz
ρ0=[max(interp(zpad(j),zb,ρb),1e-12) for j in 1:Nz+2M.NG]
p0=[max(interp(zpad(j),zb,Pb),1e-6) for j in 1:Nz+2M.NG]
θ=deg2rad(35.0); Bz0=Bstar*cos(θ); Bx0=0.5*Bstar*sin(θ); a=M.Atmos(gsurf,ρ0,p0,Bx0,Bz0)
U=M.alloc(g); R=similar(U);U1=similar(U);W=similar(U);Fx=similar(U);Fz=similar(U)
for j in 1:Nz+2M.NG,i in 1:Nx+2M.NG
    c=M.cons(γ,μ0,ρ0[j],0.,0.,0.,Bx0,0.,Bz0,p0[j],0.); for m in 1:M.NVAR; U[m,i,j]=c[m]; end; end
# seed a perturbation so the flux/limiter paths are exercised realistically
for j in M.NG+1:M.NG+Nz, i in M.NG+1:M.NG+40
    U[M.IMX,i,j] += 1e-3*U[M.IRHO,i,j]*sin(i*0.3)*exp(-((zpad(j)/(1.5Hp))^2)); end

ch=M.maxspeed(U,g,γ,μ0); dt=0.4*min(g.dx,g.dz)/ch
fluxx=zeros(Nx)

# --- the driver-level diagnostic, copied verbatim from paper_run.jl -----------
function flux_diag_orig!(fluxx,U,g,γ,μ0,dt,Nx,Nz)
    @threads for i in 1:Nx
        fx=0.0
        @inbounds for j in M.NG+1:M.NG+Nz
            ρ,vx,vy,vz,Bx,By,Bz,p,ψ=M.prim(γ,μ0,U,M.NG+i,j)
            pB=0.5*(Bx^2+By^2+Bz^2)/μ0; E=p/(γ-1)+0.5ρ*(vx^2+vy^2+vz^2)+pB
            fx+=(E+p+pB)*vx-Bx*(vx*Bx+vy*By+vz*Bz)/μ0; end
        fluxx[i]+=fx*g.dz*dt; end
end
# --- cache-friendly variant: thread over j-slabs, per-thread partial rows -----
function flux_diag_jslab!(fluxx,U,g,γ,μ0,dt,Nx,Nz,part)
    nt=size(part,2)
    fill!(part,0.0)
    @threads for tid in 1:nt
        jlo = M.NG+1 + div((tid-1)*Nz,nt); jhi = M.NG + div(tid*Nz,nt)
        @inbounds for j in jlo:jhi, i in 1:Nx
            ρ,vx,vy,vz,Bx,By,Bz,p,ψ=M.prim(γ,μ0,U,M.NG+i,j)
            pB=0.5*(Bx^2+By^2+Bz^2)/μ0; E=p/(γ-1)+0.5ρ*(vx^2+vy^2+vz^2)+pB
            part[i,tid] += (E+p+pB)*vx-Bx*(vx*Bx+vy*By+vz*Bz)/μ0
        end
    end
    @inbounds for t in 1:nt, i in 1:Nx; fluxx[i]+=part[i,t]*g.dz*dt; end
end

# --- threaded replacements for the serial broadcasts in step_grav! -----------
function axpy_thread!(U1,U,R,dt)
    @threads for j in axes(U,3)
        @inbounds for i in axes(U,2), m in 1:M.NVAR; U1[m,i,j]=U[m,i,j]+dt*R[m,i,j]; end
    end
end
function rk2_thread!(U,U1,R,dt)
    @threads for j in axes(U,3)
        @inbounds for i in axes(U,2), m in 1:M.NVAR
            U[m,i,j]=0.5*U[m,i,j]+0.5*(U1[m,i,j]+dt*R[m,i,j]); end
    end
end

bench(f,n=20) = (f(); GC.gc(); t=@elapsed(for _ in 1:n; f(); end); t/n*1e3)

M.step_grav!(U,R,U1,W,Fx,Fz,g,γ,μ0,dt,ch,a,:outflow)   # warmup
nt=Threads.nthreads(); part=zeros(Nx,nt)
flux_diag_jslab!(fluxx,U,g,γ,μ0,dt,Nx,Nz,part)          # warmup

t_step  = bench(()->M.step_grav!(U,R,U1,W,Fx,Fz,g,γ,μ0,dt,ch,a,:outflow))
t_rhs   = bench(()->M.rhs_grav!(R,U,W,Fx,Fz,g,γ,μ0,ch,a,:outflow))
t_max   = bench(()->M.maxspeed(U,g,γ,μ0))
t_bc1   = bench(()->(@. U1=U+dt*R))
t_bc2   = bench(()->(@. U=0.5*U+0.5*(U1+dt*R)))
t_bc1t  = bench(()->axpy_thread!(U1,U,R,dt))
t_bc2t  = bench(()->rk2_thread!(U,U1,R,dt))
t_diag  = bench(()->flux_diag_orig!(fluxx,U,g,γ,μ0,dt,Nx,Nz))
t_diagj = bench(()->flux_diag_jslab!(fluxx,U,g,γ,μ0,dt,Nx,Nz,part))

cells=Nx*Nz
@printf("grid %d x %d = %.2f Mcells | %d threads\n\n", Nx,Nz,cells/1e6,nt)
@printf("%-34s %9s %8s\n","phase","ms/call","%% of step")
@printf("%-34s %9.2f %8s\n","FULL step_grav! (2 RK stages)",t_step,"100%")
@printf("%-34s %9.2f %7.1f%%\n","  rhs_grav! (x2 per step)",t_rhs,200*t_rhs/t_step)
@printf("%-34s %9.2f %7.1f%%\n","  broadcast U1=U+dt*R  [SERIAL]",t_bc1,100*t_bc1/t_step)
@printf("%-34s %9.2f %7.1f%%\n","  broadcast RK2 combine [SERIAL]",t_bc2,100*t_bc2/t_step)
@printf("%-34s %9.2f %7.1f%%\n","  (residual: psi damp + bc)",t_step-2t_rhs-t_bc1-t_bc2,
        100*(t_step-2t_rhs-t_bc1-t_bc2)/t_step)
println()
@printf("%-34s %9.2f %7.1f%%   <- threaded, speedup %.1fx\n","  broadcast 1 THREADED",t_bc1t,100*t_bc1t/t_step,t_bc1/t_bc1t)
@printf("%-34s %9.2f %7.1f%%   <- threaded, speedup %.1fx\n","  broadcast 2 THREADED",t_bc2t,100*t_bc2t/t_step,t_bc2/t_bc2t)
println()
@printf("%-34s %9.2f %7.1f%%  (driver, every step)\n","maxspeed",t_max,100*t_max/t_step)
@printf("%-34s %9.2f %7.1f%%  (driver, every step)\n","Phi diagnostic (i-thread, orig)",t_diag,100*t_diag/t_step)
@printf("%-34s %9.2f %7.1f%%   <- j-slab, speedup %.1fx\n","Phi diagnostic (j-slab)",t_diagj,100*t_diagj/t_step,t_diag/t_diagj)
println()
tot_now = t_step + t_max + t_diag
tot_opt = t_step - t_bc1 - t_bc2 + t_bc1t + t_bc2t + t_max + t_diagj
@printf("driver step total now: %.2f ms   optimized: %.2f ms   => %.2fx speedup\n",
        tot_now, tot_opt, tot_now/tot_opt)
