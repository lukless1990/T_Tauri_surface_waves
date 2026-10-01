# paper_run_cont.jl — CONTINUOUS-DRIVING variant of paper_run.jl.
#   Same atmosphere/grid/domain (res×2, 0.1 R⋆), but the footpoint is driven *indefinitely*
#   (vx = M0·cf·sin(π t/τ) for ALL t, period 2τ=300 s) instead of a single pulse. This tests
#   whether sustained accretion driving changes the transport: it should reach a STEADY STATE
#   whose time-averaged energy flux ⟨F(x)⟩ decays with the SAME length L as the single pulse.
#
#   Diagnostic: over a trailing window [tend−WIN, tend] (near-field already steady), accumulate the
#   time-averaged energy flux ⟨F(x)⟩ = (1/T_win)∫F dt.  Also dumps 2D vx snapshots every 30 min.
#
# Usage: julia -t 14 SurfaceCode/paper_run_cont.jl [M0]
include(joinpath(@__DIR__,"mhd2d.jl")); using .MHD2D; const M=MHD2D
using Printf, DelimitedFiles
using Base.Threads: @threads
const Gc=6.674e-11; const Msun=1.989e30; const Rsun=6.957e8
const Mstar=0.5*Msun; const Rstar=2.0*Rsun; const gsurf=Gc*Mstar/Rstar^2
const μ0c=4π*1e-7; const γc=5/3; const Bstar=0.1

const M0_DEFAULT = 2.63          # same footpoint drive as the pulse run (~10 km/s launched amplitude)
const SNAP_DT    = 1800.0        # 2D snapshot cadence [s]
const LXFRAC     = 0.1
const CPH        = 24            # res×2
const TEXTRA     = 1.05          # run to 1.05×(front transit) — same tend as the pulse run
const WIN        = 1800.0        # trailing steady-state averaging window [s] = 6 driving periods

read_bg(f)=(d=readdlm(f;comments=true,comment_char='#'); p=sortperm(d[:,1]); (d[p,1],d[p,3],d[p,5]))
function interp(xq,x,y); xq≤x[1] && return y[1]; xq≥x[end] && return y[end]
    k=searchsortedlast(x,xq); w=(xq-x[k])/(x[k+1]-x[k]); (1-w)*y[k]+w*y[k+1]; end

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
    # driver half-period: see paper_run.jl. TAU_S=50 is the clump-consistent value (r_c ~ 1e-2 R⋆).
    τ = haskey(ENV,"TAU_S") ? parse(Float64,ENV["TAU_S"]) : 150.0
    zw(j)=exp(-(zpad(j)/(1.5Hp))^2)
    function injf(U,gg,aa,γ,μ0,t)                       # CONTINUOUS: no (0≤t≤τ) gate
        vx=M0*cf*sin(π*t/τ)
        @inbounds for j in 1:gg.Nz+2M.NG,k in 1:M.NG; il=M.NG+1-k; ρ=aa.ρ0[j]; ux=vx*zw(j)
            c=M.cons(γ,μ0,ρ,ux,0.,0.,aa.Bx0[j],0.,aa.Bz0,aa.p0[j],0.); for m in 1:M.NVAR; U[m,il,j]=c[m]; end; end
    end

    label = get(ENV,"RUN_LABEL","paper_run_cont")
    outdir=joinpath(@__DIR__,"..","output",label); mkpath(outdir)
    writedlm(joinpath(outdir,"grid_x.txt"), collect(g.xc)./Rstar)
    writedlm(joinpath(outdir,"grid_z.txt"), [zpad(M.NG+j)/Hp for j in 1:Nz])
    tend=TEXTRA*(Nx*g.dx)/cf+τ; twin0=tend-WIN
    open(joinpath(outdir,"meta.txt"),"w") do io
        @printf(io,"mode=continuous\nM0=%.3f\nboundary_drive_kms=%.3f\ncf_kms=%.3f\ndrive_period_s=%.1f\nLx_Rstar=%.3f\ndx_km=%.1f\nHp_km=%.1f\nNx=%d\nNz=%d\ntend_s=%.0f\nwin0_s=%.0f\n",
                M0,M0*cf/1e3,cf/1e3,2τ,LXFRAC,g.dx/1e3,Hp/1e3,Nx,Nz,tend,twin0)
    end
    snaptimes=Float64[]; snappeaks=Float64[]; si=0
    Favg=zeros(Nx); Twin=0.0                            # trailing-window time-averaged energy flux
    function save_snapshot(tag,t)
        vx2d=[U[M.IMX,M.NG+i,M.NG+j]/max(U[M.IRHO,M.NG+i,M.NG+j],1e-14)/1e3 for i in 1:Nx, j in 1:Nz]
        writedlm(joinpath(outdir,@sprintf("snap_vx_%03d.txt",tag)), vx2d)
        push!(snaptimes,t); push!(snappeaks,maximum(abs,vx2d))
        writedlm(joinpath(outdir,"snap_times.txt"), hcat(snaptimes,snappeaks))
        @printf("   >> snapshot %03d  t=%.0f s (%.2f h)  domain peak |vx|=%.2f km/s\n",tag,t,t/3600,snappeaks[end]); flush(stdout)
    end

    t=0.0; ns=0; next_snap=SNAP_DT; trig_done=false      # early snap at 1.5τ, then every SNAP_DT (match pulse run)
    @printf("CONTINUOUS run: M0=%.2f (drive %.1f km/s, period %.0f s), Nx=%d Nz=%d, Lx=%.2f Rstar, dx=%.0f km, cf=%.1f km/s\n",
            M0,M0*cf/1e3,2τ,Nx,Nz,LXFRAC,g.dx/1e3,cf/1e3)
    @printf("           tend=%.0f s (%.2f h); steady window [%.0f,%.0f] s -> %s\n",tend,tend/3600,twin0,tend,outdir); flush(stdout)
    while t<tend
        ch=M.maxspeed(U,g,γ,μ0); dt=0.4*min(g.dx,g.dz)/ch
        M.step_grav!(U,R,U1,W,Fx,Fz,g,γ,μ0,dt,ch,a,:outflow; inject=injf,t=t); t+=dt; ns+=1
        if t≥twin0                                       # accumulate steady-state time-averaged flux
            @threads for i in 1:Nx
                fx=0.0
                @inbounds for j in M.NG+1:M.NG+Nz
                    ρ,vx,vy,vz,Bx,By,Bz,p,ψ=M.prim(γ,μ0,U,M.NG+i,j)
                    pB=0.5*(Bx^2+By^2+Bz^2)/μ0; E=p/(γ-1)+0.5ρ*(vx^2+vy^2+vz^2)+pB
                    fx+=(E+p+pB)*vx-Bx*(vx*Bx+vy*By+vz*Bz)/μ0; end
                Favg[i]+=fx*g.dz*dt; end
            Twin+=dt
        end
        if !trig_done && t≥1.5τ; si+=1; save_snapshot(si,t); trig_done=true; end
        if t≥next_snap; si+=1; save_snapshot(si,t); next_snap+=SNAP_DT; end
        if ns%2000==0
            vmax=0.0; @inbounds for j in M.NG+1:M.NG+Nz,i in M.NG+1:M.NG+Nx
                vmax=max(vmax,abs(U[M.IMX,i,j]/max(U[M.IRHO,i,j],1e-14))); end
            @printf("   step %6d  t=%6.0f/%.0f s (%3.0f%%)  dt=%.2fs  max|vx|=%.2f km/s%s\n",
                    ns,t,tend,100t/tend,dt,vmax/1e3, t≥twin0 ? "  [avg]" : "")
            flush(stdout); isfinite(vmax)||(println("   ** NaN/Inf — aborting **");break)
        end
    end
    Favg ./= max(Twin,1e-30)                              # -> time-averaged flux ⟨F(x)⟩
    writedlm(joinpath(outdir,"flux_steady.txt"), hcat(collect(g.xc)./Rstar, Favg))
    si+=1; save_snapshot(si,t)
    @printf("done: %d steps, %d snapshots, t=%.0f s (%.2f h), window=%.0f s. data -> %s\n", ns, si, t, t/3600, Twin, outdir)
end
main()
