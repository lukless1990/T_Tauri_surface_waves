# paper_run_leak.jl — paper_run.jl + a VERTICAL-LEAKAGE ledger (referee point 1, 2026-09-30).
#   Identical physics/numerics to paper_run.jl (the RK2 step is written out here only so the
#   numerical z-fluxes of BOTH stages can be read; it is the same sequence as M.step_grav!).
#   Extra output: cumulative energy and mass that crossed the TOP (z=+3Hp) and BOTTOM (z=-4Hp)
#   faces, per column x, in leak_NNN.txt (cols: x/R⋆, E_top, E_bot, M_top, M_bot) and leak_final.txt.
#   Units: time-integrated flux densities, J/m² and kg/m² (positive = leaving the domain).
#   Multiply by dx and sum over x to compare with Φ [J/m]. E is the conserved MHD total energy
#   (no potential energy); the potential-energy flux is g·z_face × mass flux (offline).
#   Runs until converged: flux-burial point (Φ/Φ0 < 7.8e-4) moves < CONV_TOL per hour, checked
#   every 30 min after the front has crossed the domain; hard cap TMAX_H (stellar hours).
# ---- original header ----
# paper_run.jl — production MHD fast-mode run for the paper figure.
#   * res×2 (dx = 50 km, cph = 24)
#   * lateral domain out to 0.1 R⋆
#   * runs until the pulse has travelled the full 0.1 R⋆ (tend = 1.05·Lx/cf + τ)
#   * driving amplitude M0 calibrated so the launched domain pulse peaks at ~10 km/s
#   * dumps a 2D v_x snapshot + the running energy flux Φ(x) every 30 min of stellar time
#     into output/paper_run/  (plus a "trigger" snapshot at 1.5τ to capture the launch).
#
# Usage: julia -t 14 SurfaceCode/paper_run.jl [M0]      (M0 default set from calibration below)
include(joinpath(@__DIR__,"mhd2d.jl")); using .MHD2D; const M=MHD2D
using Printf, DelimitedFiles
using Base.Threads: @threads
const Gc=6.674e-11; const Msun=1.989e30; const Rsun=6.957e8
const Mstar=0.5*Msun; const Rstar=2.0*Rsun; const gsurf=Gc*Mstar/Rstar^2
const μ0c=4π*1e-7; const γc=5/3
const Bstar = haskey(ENV,"BSTAR_T") ? parse(Float64,ENV["BSTAR_T"]) : 0.1   # T (0.1=1kG default, 0.01=0.1kG)

const M0_DEFAULT = 2.63         # calibrated for ~10 km/s launched-pulse peak@1.5τ (see calib_pulse.jl); overwrite via ARGS[1]
const SNAP_DT    = 1800.0       # snapshot cadence [s] = 30 min stellar time
const LXFRAC     = 0.1          # lateral domain in R⋆
const CPH        = 24           # cells per H_p  (res×2)
const TEXTRA     = 1.05         # run to 1.05×(pulse transit) so it just clears 0.1 R⋆
const SP_ON      = get(ENV,"LEFT_SPONGE","1") == "1"   # left absorbing layer (footpoint anti-reflection)
const SP_NCELL   = 40           # layer width [cells] (~0.002 R⋆; well left of the x>0.004 decay region)
const SP_RATE    = 0.1          # absorption rate [1/s]
const CONV_TOL   = haskey(ENV,"CONV_TOL") ? parse(Float64,ENV["CONV_TOL"]) : 0.01   # burial drift per hour
const TMAX_H     = haskey(ENV,"TMAX_H")   ? parse(Float64,ENV["TMAX_H"])   : 30.0   # hard cap [stellar h]
const PHI_AMB    = 7.8e-4                                                            # (v_amb/v_inj)^2

read_bg(f)=(d=readdlm(f;comments=true,comment_char='#'); p=sortperm(d[:,1]); (d[p,1],d[p,3],d[p,5]))
function interp(xq,x,y); xq≤x[1] && return y[1]; xq≥x[end] && return y[end]
    k=searchsortedlast(x,xq); w=(xq-x[k])/(x[k+1]-x[k]); (1-w)*y[k]+w*y[k+1]; end

# Operator-split left absorbing layer: over the leftmost nsp cells, damp the VELOCITY toward 0 and the
# INTERNAL ENERGY toward the hydrostatic background (p0), while leaving ρ and B UNTOUCHED. Damping ρ
# would make the layer a mass sink (draws a spurious inflow); velocity+eint damping absorbs the wave
# (both its kinetic and thermal parts) without moving mass. ξ² ramp for a smooth, low-reflection profile.
function left_sponge!(U, g, a, γ, μ0, nsp, σ0, dt)
    @inbounds for j in 1:g.Nz+2*M.NG
        eint_bg = a.p0[j]/(γ-1)
        for i in M.NG+1:M.NG+nsp
            ξ = (nsp - (i-M.NG)) / nsp
            f = σ0*ξ^2*dt; f = f/(1+f); g1 = 1-f
            ρ = max(U[M.IRHO,i,j], 1e-14)
            mx,my,mz = U[M.IMX,i,j], U[M.IMY,i,j], U[M.IMZ,i,j]
            pB = 0.5*(U[M.IBX,i,j]^2+U[M.IBY,i,j]^2+U[M.IBZ,i,j]^2)/μ0
            ke = 0.5*(mx^2+my^2+mz^2)/ρ
            eint = U[M.IE,i,j] - ke - pB
            U[M.IMX,i,j]=g1*mx; U[M.IMY,i,j]=g1*my; U[M.IMZ,i,j]=g1*mz   # v -> 0
            U[M.IE,i,j] = (g1*eint + f*eint_bg) + g1^2*ke + pB           # eint -> p0/(γ-1); new KE; ρ,B kept
        end
    end
end

# Right absorbing layer (mirror of left_sponge!): lets the pulse EXIT through x=Lx without reflecting off
# the outflow boundary, which otherwise fluctuates the flux near 0.1 R⋆. Always on — the region is
# quiescent background until the pulse arrives (v=0, eint=p0 ⇒ no-op there). Same velocity+eint damping.
function right_sponge!(U, g, a, γ, μ0, nsp, σ0, dt)
    @inbounds for j in 1:g.Nz+2*M.NG
        eint_bg = a.p0[j]/(γ-1)
        for i in M.NG+g.Nx-nsp+1 : M.NG+g.Nx
            ξ = (i - (M.NG+g.Nx-nsp)) / nsp        # 0 at inner edge, ->1 at the right boundary
            f = σ0*ξ^2*dt; f = f/(1+f); g1 = 1-f
            ρ = max(U[M.IRHO,i,j], 1e-14)
            mx,my,mz = U[M.IMX,i,j], U[M.IMY,i,j], U[M.IMZ,i,j]
            pB = 0.5*(U[M.IBX,i,j]^2+U[M.IBY,i,j]^2+U[M.IBZ,i,j]^2)/μ0
            ke = 0.5*(mx^2+my^2+mz^2)/ρ
            eint = U[M.IE,i,j] - ke - pB
            U[M.IMX,i,j]=g1*mx; U[M.IMY,i,j]=g1*my; U[M.IMZ,i,j]=g1*mz
            U[M.IE,i,j] = (g1*eint + f*eint_bg) + g1^2*ke + pB
        end
    end
end

function main()
    γ=γc; μ0=μ0c
    M0 = length(ARGS)≥1 ? parse(Float64,ARGS[1]) : M0_DEFAULT
    zb,ρb,Pb=read_bg(joinpath(@__DIR__,"background_35deg.txt"))
    ρ0s=interp(0.,zb,ρb); P0s=interp(0.,zb,Pb); Hp=P0s/(ρ0s*gsurf)
    Lx=LXFRAC*Rstar; dxt=1e5*12/CPH
    Nx=round(Int,Lx/dxt); zlo=-4Hp; Lz=7Hp; Nz=max(64,round(Int,Lz/(Hp/CPH)))
    g=M.Grid(Nx,Nz,Lx,Lz; z0=zlo); zpad(j)=zlo+(j-M.NG-0.5)*g.dz
    ρ0=[max(interp(zpad(j),zb,ρb),1e-12) for j in 1:Nz+2M.NG]
    p0=[max(interp(zpad(j),zb,Pb),1e-6) for j in 1:Nz+2M.NG]
    θ=deg2rad(35.0); Bz0=Bstar*cos(θ); Bx0=0.5*Bstar*sin(θ); a=M.Atmos(gsurf,ρ0,p0,Bx0,Bz0)
    cf=sqrt(γ*P0s/ρ0s+(Bx0^2+Bz0^2)/(μ0*ρ0s))
    U=M.alloc(g); R=similar(U);U1=similar(U);W=similar(U);Fx=similar(U);Fz=similar(U)
    for j in 1:Nz+2M.NG,i in 1:Nx+2M.NG
        c=M.cons(γ,μ0,ρ0[j],0.,0.,0.,Bx0,0.,Bz0,p0[j],0.); for m in 1:M.NVAR; U[m,i,j]=c[m]; end; end
    # driver duration: 150 s default (3× the Cranmer clump passage time — the wave-favorable choice);
    # TAU_S=50 reproduces the r_c ~ 1e-2 R⋆ clump directly (paper Sect. 2.2 / 3.4).
    τ = haskey(ENV,"TAU_S") ? parse(Float64,ENV["TAU_S"]) : 150.0
    zw(j)=exp(-(zpad(j)/(1.5Hp))^2)
    function injf(U,gg,aa,γ,μ0,t)
        (0≤t≤τ)||return; vx=M0*cf*sin(π*t/τ)
        @inbounds for j in 1:gg.Nz+2M.NG,k in 1:M.NG; il=M.NG+1-k; ρ=aa.ρ0[j]; ux=vx*zw(j)
            c=M.cons(γ,μ0,ρ,ux,0.,0.,aa.Bx0[j],0.,aa.Bz0,aa.p0[j],0.); for m in 1:M.NVAR; U[m,il,j]=c[m]; end; end
    end

    label = get(ENV,"RUN_LABEL","paper_run")            # output subdir (label runs by field, etc.)
    outdir=joinpath(@__DIR__,"..","output",label); mkpath(outdir)
    writedlm(joinpath(outdir,"grid_x.txt"), collect(g.xc)./Rstar)      # x/R⋆  (length Nx)
    writedlm(joinpath(outdir,"grid_z.txt"), [zpad(M.NG+j)/Hp for j in 1:Nz])  # z/H_p (length Nz)
    B0mag=sqrt(Bx0^2+Bz0^2)
    open(joinpath(outdir,"meta.txt"),"w") do io
        @printf(io,"M0=%.3f\nboundary_drive_kms=%.3f\ncf_kms=%.3f\nLx_Rstar=%.3f\ndx_km=%.1f\nHp_km=%.1f\nNx=%d\nNz=%d\nsnap_dt_s=%.0f\nBstar_kG=%.3f\nB0_G=%.1f\nleft_sponge=%d\nsp_ncell=%d\ntau_s=%.1f\n",
                M0, M0*cf/1e3, cf/1e3, LXFRAC, g.dx/1e3, Hp/1e3, Nx, Nz, SNAP_DT, Bstar*10, B0mag*1e4, SP_ON ? 1 : 0, SP_ON ? SP_NCELL : 0, τ)
    end
    snaptimes=Float64[]; snappeaks=Float64[]; si=0
    fluxx=zeros(Nx)
    Etop=zeros(Nx); Ebot=zeros(Nx); Mtop=zeros(Nx); Mbot=zeros(Nx)   # cumulative, leaving = positive
    jb=M.NG; jt=M.NG+Nz                                              # bottom / top boundary faces in Fz
    function acc_zflux!(w)                                           # w = dt/2 (Heun stage weight)
        @inbounds for i in 1:Nx
            Ebot[i]-=w*Fz[M.IE,M.NG+i,jb];   Etop[i]+=w*Fz[M.IE,M.NG+i,jt]
            Mbot[i]-=w*Fz[M.IRHO,M.NG+i,jb]; Mtop[i]+=w*Fz[M.IRHO,M.NG+i,jt]
        end
    end
    xc=collect(g.xc)./Rstar
    function burial()                                                # flux-burial point of current Φ
        m0=(xc.>0.002).&(xc.<0.004); Φ0=maximum(fluxx[m0])
        for i in eachindex(xc); xc[i]>0.01 && fluxx[i]/Φ0<PHI_AMB && return xc[i]; end
        return NaN
    end
    function save_leak(tag)
        writedlm(joinpath(outdir,@sprintf("leak_%03d.txt",tag)), hcat(xc,Etop,Ebot,Mtop,Mbot))
    end
    function save_snapshot(tag,t)
        vx2d=[U[M.IMX,M.NG+i,M.NG+j]/max(U[M.IRHO,M.NG+i,M.NG+j],1e-14)/1e3 for i in 1:Nx, j in 1:Nz]
        writedlm(joinpath(outdir,@sprintf("snap_vx_%03d.txt",tag)), vx2d)
        writedlm(joinpath(outdir,@sprintf("snap_flux_%03d.txt",tag)), hcat(collect(g.xc)./Rstar, copy(fluxx)))
        push!(snaptimes,t); push!(snappeaks,maximum(abs,vx2d))
        writedlm(joinpath(outdir,"snap_times.txt"), hcat(snaptimes,snappeaks))  # cols: t[s], domain peak |vx| [km/s]
        @printf("   >> snapshot %03d  t=%.0f s (%.2f h)  domain peak |vx|=%.2f km/s\n",tag,t,t/3600,snappeaks[end]); flush(stdout)
    end

    tcross=TEXTRA*(Nx*g.dx)/cf+τ                  # original paper_run end time (front has crossed)
    tend=TMAX_H*3600; t=0.0; ns=0; next_snap=SNAP_DT; trig_done=false  # trigger snap at 1.5τ, then every SNAP_DT
    bhist=Tuple{Float64,Float64}[]                 # (t, burial x) at each 30-min snapshot
    open(joinpath(outdir,"meta.txt"),"a") do io
        @printf(io,"leak_ledger=1\nconv_tol_per_h=%.4f\ntmax_h=%.1f\nt_cross_s=%.1f\n",CONV_TOL,TMAX_H,tcross); end
    @printf("paper run: M0=%.2f (drive %.1f km/s), tau=%.0f s, Nx=%d Nz=%d, Lx=%.2f Rstar, dx=%.0f km, cf=%.1f km/s\n",
            M0,M0*cf/1e3,τ,Nx,Nz,LXFRAC,g.dx/1e3,cf/1e3)
    @printf("           tend=%.0f s (%.2f h), snapshot every %.0f s (%.0f min) -> %s\n",
            tend,tend/3600,SNAP_DT,SNAP_DT/60,outdir); flush(stdout)
    while t<tend
        ch=M.maxspeed(U,g,γ,μ0); dt=0.4*min(g.dx,g.dz)/ch
        # == M.step_grav!(U,R,U1,W,Fx,Fz,g,γ,μ0,dt,ch,a,:outflow; inject=injf,t=t), written out
        #    so the boundary z-fluxes of both Heun stages enter the ledger
        M.rhs_grav!(R,U,W,Fx,Fz,g,γ,μ0,ch,a,:outflow; inject=injf,t=t);     acc_zflux!(0.5dt); @. U1=U+dt*R
        M.rhs_grav!(R,U1,W,Fx,Fz,g,γ,μ0,ch,a,:outflow; inject=injf,t=t+dt); acc_zflux!(0.5dt); @. U=0.5*U+0.5*(U1+dt*R)
        αψ=exp(-0.18*ch*dt/min(g.dx,g.dz))
        @threads for j in M.NG+1:M.NG+g.Nz
            @inbounds for i in M.NG+1:M.NG+g.Nx; U[M.IPSI,i,j]*=αψ; end
        end
        t+=dt; ns+=1
        SP_ON && t > 4τ && left_sponge!(U, g, a, γ, μ0, SP_NCELL, SP_RATE, dt)   # left: absorb after pulse clears
        SP_ON && right_sponge!(U, g, a, γ, μ0, SP_NCELL, SP_RATE, dt)             # right: let the pulse exit cleanly
        # accumulate the height/time-integrated x energy flux Φ(x)
        @threads for i in 1:Nx
            fx=0.0
            @inbounds for j in M.NG+1:M.NG+Nz
                ρ,vx,vy,vz,Bx,By,Bz,p,ψ=M.prim(γ,μ0,U,M.NG+i,j)
                pB=0.5*(Bx^2+By^2+Bz^2)/μ0; E=p/(γ-1)+0.5ρ*(vx^2+vy^2+vz^2)+pB
                fx+=(E+p+pB)*vx-Bx*(vx*Bx+vy*By+vz*Bz)/μ0; end
            fluxx[i]+=fx*g.dz*dt; end
        # trigger snapshot (captures the launched pulse) then the 30-min cadence
        if !trig_done && t≥1.5τ; si+=1; save_snapshot(si,t); trig_done=true; next_snap=SNAP_DT; end
        if t≥next_snap; si+=1; save_snapshot(si,t); save_leak(si); next_snap+=SNAP_DT
            xb=burial(); push!(bhist,(t,xb))
            @printf("   >> leak so far: top %.3e  bottom %.3e J/m (∫dx);  burial x=%.4f R⋆\n",
                    sum(Etop)*g.dx, sum(Ebot)*g.dx, xb); flush(stdout)
            # converged: past the crossing time, and burial moved < CONV_TOL per hour over the last hour
            if t>tcross && length(bhist)≥3
                (t1,x1)=bhist[end-2]; drift=abs(xb-x1)/x1/((t-t1)/3600)
                @printf("   >> burial drift %.2f%%/h (tol %.2f%%/h)\n",100drift,100CONV_TOL); flush(stdout)
                isfinite(drift) && drift<CONV_TOL && (println("   ** converged **"); break)
            end
        end
        if ns%2000==0
            vmax=0.0; @inbounds for j in M.NG+1:M.NG+Nz,i in M.NG+1:M.NG+Nx
                vmax=max(vmax,abs(U[M.IMX,i,j]/max(U[M.IRHO,i,j],1e-14))); end
            @printf("   step %6d  t=%6.0f/%.0f s (%3.0f%%)  dt=%.2fs  max|vx|=%.2f km/s\n",ns,t,tend,100t/tend,dt,vmax/1e3)
            flush(stdout); isfinite(vmax)||(println("   ** NaN/Inf — aborting **");break)
        end
    end
    # final flux + a closing snapshot at the domain edge
    si+=1; save_snapshot(si,t)
    writedlm(joinpath(outdir,"flux_final.txt"), hcat(collect(g.xc)./Rstar, fluxx))
    writedlm(joinpath(outdir,"leak_final.txt"), hcat(xc,Etop,Ebot,Mtop,Mbot))
    writedlm(joinpath(outdir,"burial_history.txt"), hcat(first.(bhist),last.(bhist)))
    @printf("done: %d steps, %d snapshots, t=%.0f s (%.2f h). data -> %s\n", ns, si, t, t/3600, outdir)
end
main()
