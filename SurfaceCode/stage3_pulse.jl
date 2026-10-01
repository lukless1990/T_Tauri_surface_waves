# ============================================================================
#  stage3_pulse.jl — SurfaceCode pilot, Stage 3 (hydro).
#
#  Ingest the real exported T Tauri near-surface atmosphere, inject an impulsive
#  compressive pulse at the footpoint (x=0), propagate it laterally, and measure
#  how its amplitude decays with distance x. Comparing a STRONG (shock-forming)
#  vs a WEAK (linear) pulse isolates the *shock-dissipation* decay length — the
#  first term Cranmer omits. (MHD fast mode + ambipolar come in Stages 4–5.)
#
#  Usage: julia SurfaceCode/stage3_pulse.jl [Lx_frac] [M0]
# ============================================================================
include(joinpath(@__DIR__, "patch2d.jl"))
using .Patch2D: NG, Grid, Atmos, alloc, step_grav!, max_dt, IRHO, IMX, IMZ, IE, prim, pressure
using Printf, DelimitedFiles

const G_grav = 6.674e-11; const Msun = 1.989e30; const Rsun = 6.957e8
const Mstar = 0.5*Msun; const Rstar = 2.0*Rsun
const gsurf = G_grav*Mstar/Rstar^2        # ≈34.3 m/s² at τ=1 (const over the thin patch)
const γ = 5/3

# --- read the background_*.txt exported by export_background.jl --------------
function read_background(file)
    data = readdlm(file; comments=true, comment_char='#')
    z = data[:,1]; ρ = data[:,3]; T = data[:,4]; P = data[:,5]   # SI
    p = sortperm(z)
    (z[p], ρ[p], T[p], P[p])
end
lininterp(xq, x, y) = begin
    xq ≤ x[1]   && return y[1]
    xq ≥ x[end] && return y[end]
    k = searchsortedlast(x, xq); w = (xq-x[k])/(x[k+1]-x[k])
    (1-w)*y[k] + w*y[k+1]
end

# --- build the 2D atmosphere (background replicated in x) --------------------
function build_atmos(file; Nx=500, Lx=5e7, zbot_Hp=-4.0, ztop_Hp=6.0, cph=12)
    zb, ρb, Tb, Pb = read_background(file)
    # photospheric scale height (z≈0) sets the vertical range/resolution
    ρ0s = lininterp(0.0, zb, ρb); P0s = lininterp(0.0, zb, Pb)
    Hp = P0s/(ρ0s*gsurf)
    zlo = zbot_Hp*Hp; zhi = ztop_Hp*Hp; Lz = zhi-zlo
    Nz = max(64, round(Int, Lz/(Hp/cph)))         # cph cells per surface scale height
    g = Grid(Nx, Nz, Lx, Lz; z0=zlo)
    zpad(j) = zlo + (j-NG-0.5)*g.dz
    ρ0 = [max(lininterp(zpad(j), zb, ρb), 1e-12) for j in 1:Nz+2NG]
    p0 = [max(lininterp(zpad(j), zb, Pb), 1e-6) for j in 1:Nz+2NG]
    a = Atmos(gsurf, ρ0, p0)
    csurf = sqrt(γ*P0s/ρ0s)
    (g, a, Hp, csurf, ρ0s, P0s)
end

# --- pulse: compressive velocity at the x-lo ghost, localized near z≈0 --------
const TAU_S = haskey(ENV,"TAU_S") ? parse(Float64,ENV["TAU_S"]) : 150.0   # driver duration [s]

function run_pulse(g::Grid, a::Atmos, csurf, Hp; M0=2.0, τ=TAU_S, cfl=0.4, textra=1.2)
    U=alloc(g); R=similar(U); U1=similar(U); W=similar(U)
    @inbounds for j in 1:g.Nz+2NG, i in 1:g.Nx+2NG
        U[IRHO,i,j]=a.ρ0[j]; U[IMX,i,j]=0; U[IMZ,i,j]=0; U[IE,i,j]=a.p0[j]/(γ-1)
    end
    zpad(j)= g.zc[1] - (NG-0)*g.dz + (j-1)*g.dz    # padded-index → z (j=NG+1 ↔ zc[1])
    zw(j)  = exp(-(zpad(j)/(1.5*Hp))^2)            # localize pulse near z≈0 (photosphere)
    function inject!(U, gg, aa, γ, t)
        (0.0 ≤ t ≤ τ) || return
        u = M0*csurf*sin(π*t/τ)
        @inbounds for j in 1:gg.Nz+2NG, k in 1:NG
            il=NG+1-k
            ρ=aa.ρ0[j]; ux=u*zw(j)
            U[IRHO,il,j]=ρ; U[IMX,il,j]=ρ*ux; U[IMZ,il,j]=0
            U[IE,il,j]=aa.p0[j]/(γ-1)+0.5*ρ*ux^2
        end
    end
    tend=textra*(g.Nx*g.dx)/csurf + τ
    ampx=zeros(g.Nx)                                # peak |u_x| per x-column
    fluxx=zeros(g.Nx)                               # ∫∫ (E+p)u_x dz dt — height/time-integrated energy flux
    t=0.0; nstep=0; dt=max_dt(U,g,γ,cfl)
    while t<tend
        dt=max_dt(U,g,γ,cfl)                        # adaptive: shock raises |u|+c above the static value
        step_grav!(U,R,U1,W,g,γ,a,dt,:outflow,:atmos; inject=inject!, t=t); t+=dt; nstep+=1
        @inbounds for i in 1:g.Nx
            m=0.0; fx=0.0
            for j in NG+1:NG+g.Nz
                ρ=max(U[IRHO,NG+i,j],1e-14); ux=U[IMX,NG+i,j]/ρ
                p=pressure(γ,ρ,U[IMX,NG+i,j],U[IMZ,NG+i,j],U[IE,NG+i,j])
                m=max(m, abs(ux))
                fx += (U[IE,NG+i,j]+p)*ux            # x energy flux (0 in the unperturbed atmosphere)
            end
            ampx[i]=max(ampx[i], m)
            fluxx[i]+= fx*g.dz*dt                    # accumulate energy transported past x
        end
    end
    (collect(g.xc), ampx, fluxx, nstep, dt)
end

# --- decay-length fit: log(amp) ~ a - x/L over the range where the pulse is clean
function decay_length(x, amp; i0frac=0.15, i1frac=0.9)
    n=length(x); i0=max(2,round(Int,i0frac*n)); i1=round(Int,i1frac*n)
    xs=x[i0:i1]; ys=log.(max.(amp[i0:i1], 1e-30))
    X=hcat(ones(length(xs)), xs); β=X\ys           # least squares
    slope=β[2]
    L = slope<0 ? -1/slope : Inf
    (L, exp(β[1]))
end

# ---------------------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    Lxfrac = length(ARGS)≥1 ? parse(Float64,ARGS[1]) : 0.036
    bgfile = joinpath(@__DIR__, "background_35deg.txt")
    Lx = Lxfrac*Rstar
    g, a, Hp, csurf, ρ0s, P0s = build_atmos(bgfile; Nx=500, Lx=Lx)
    @printf("atmosphere: Nx=%d Nz=%d, Lx=%.2e m (%.3f Rstar), Hp=%.0f km, c_s(surf)=%.1f km/s, rho_s=%.2e\n",
            g.Nx, g.Nz, Lx, Lxfrac, Hp/1e3, csurf/1e3, ρ0s)
    @printf("  %-16s %-22s %-22s\n","pulse","L(peak-amp)","L(energy-flux Φ)")
    for (lbl,M0) in (("strong (shock)",2.0), ("weak (linear)",0.05))
        x, amp, flux, nstep, dt = run_pulse(g, a, csurf, Hp; M0=M0)
        La,_ = decay_length(x, amp); Lf,_ = decay_length(x, flux)
        @printf("  %-16s %8.4f Rstar (%.2e m)   %8.4f Rstar (%.2e m)   [%d steps]\n",
                lbl, La/Rstar, La, Lf/Rstar, Lf, nstep)
    end
    # convergence check: does L(Φ) for the strong pulse change with resolution?
    println("\nconvergence of L(Φ), strong pulse (numerical diffusion → L grows with resolution):")
    for res in (1,2)
        g2,a2,Hp2,cs2,_,_ = build_atmos(bgfile; Nx=500*res, Lx=Lx, cph=12*res)
        x,_,flux,ns,_ = run_pulse(g2,a2,cs2,Hp2; M0=2.0)
        Lf,_ = decay_length(x,flux)
        @printf("  res×%d (Nx=%d Nz=%d, dx=%.0f km): L(Φ)=%.4f Rstar (%.2e m) [%d steps]\n",
                res, g2.Nx, g2.Nz, g2.dx/1e3, Lf/Rstar, Lf, ns)
    end
    println("(if L(Φ) grows markedly with resolution, the decay is numerical, not physical — needs")
    println(" further refinement; the physical shock length is the converged value.)")
end
