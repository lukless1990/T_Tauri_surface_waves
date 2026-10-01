# ============================================================================
#  solar_moreton.jl — SurfaceCode/SolarRuns: a flare-triggered Moreton wave on
#  the quiet Sun, run with the validated MHD2D solver core (../src/mhd2d.jl).
#
#  Scenario (blast-wave variant; all values + references in solar_parameters.pdf):
#  an impulsive flare pressure pulse in the low corona launches a fast-mode shock;
#  its lower flank sweeps the chromosphere — the Moreton front (Uchida 1968;
#  Vršnak & Cliver 2008 blast scenario; setup after Krause et al. 2015).
#  Uniform, vertical quiet-Sun B0 (default 3 G): lateral propagation is a
#  PERPENDICULAR fast wave, ambient V_f = sqrt(c_s²+V_A²) ≈ 550-800 km/s — the
#  speed observed Moreton fronts decelerate toward (Warmuth et al. 2004).
#  Two backgrounds via MODE (see below): gravity-free Krause slab (production)
#  or full hydrostatic VAL-C stratification (high-res variant).
#
#  Numerics beyond the parent core (all in this file, core untouched):
#  per-stage bounded-state protection (relative ρ floor + velocity ceiling +
#  eint clamp on flagged cells — HLL/MUSCL has no positivity guarantee at the
#  slab/TR interface), sacrificial guard + post-pulse relaxer columns at the
#  flare site, mass-conserving boundary sponges (paper_run.jl design).
#
#  Diagnostics written to SolarRuns/run_<label>/:
#    background.txt              z, ρ0, p0, T0, c_s, V_A, V_f
#    snap_vz_###.txt / snap_vx_###.txt / snap_dpp_###.txt   2D fields (x-subsampled)
#    trace_chromo_vz.txt         v_z(x,t) at z=Z_HA  — the "Hα proxy" stack
#    trace_corona_dpp.txt        Δp/p0(x,t) at z=Z_EIT — the coronal (EUV) front
#    trace_times.txt, snap_times.txt, meta.txt
#
#  Usage:  julia -t 14 SurfaceCode/SolarRuns/solar_moreton.jl
#     env: MODE (krause|valc), A_PULSE, B0_G, THETA_B_DEG, Z_FL_MM, TAU_FL_S,
#          LX_MM, LZ_MM, DX_KM, DZ_KM, TEND_S, RUN_LABEL, ...(see getf calls)
# ============================================================================
include(joinpath(@__DIR__,"..","src","mhd2d.jl")); using .MHD2D; const M=MHD2D
using Printf, DelimitedFiles

# --- physical constants ------------------------------------------------------
const kB=1.380649e-23; const mH=1.6735575e-27
const μ0c=4π*1e-7;     const γc=5/3
const gsun=274.0                      # m/s² (surface; -8% at z=30 Mm — neglected)
const Rsun=6.957e8

# --- adopted solar-atmosphere parameters (references in solar_parameters.pdf) -
const T_PHOT   = 6420.0               # K   VAL-C kinetic T at tau_500=1 (NOT T_eff)
const P_PHOT   = 1.172e4              # Pa  gas pressure at tau_500=1 (VAL-C)
const T_MINV   = 4170.0               # K   temperature minimum (VAL-C: 4170 K)
const Z_MINV   = 0.515e6              # m   height of the temperature minimum (VAL-C)
const T_CHROM  = 7000.0               # K   upper-chromospheric plateau (VAL-C ~7160 @2 Mm)
const Z_CHTOP  = 2.0e6                # m   top of the chromospheric rise
const Z_TRTOP  = parse(Float64,get(ENV,"Z_TRTOP_MM","3.5"))*1e6  # m  top of the numerically
                                           # BROADENED transition region (Lionello-style): compact-
                                           # support smoothstep from Z_CHTOP to Z_TRTOP. The real
                                           # ~0.28 Mm TR at 2.27-2.54 Mm is unresolvable at dz>=100 km,
                                           # and tanh tails would leak MK temperatures into the
                                           # chromosphere — smoothstep keeps z<Z_CHTOP exactly VAL-C.
const T_COR    = 1.5e6                # K   quiet-corona temperature
const GRAD_SUB = 1.0e-2               # K/m sub-photospheric gradient (anchor layer)

# --- run controls -------------------------------------------------------------
# MODE selects the background:
#   "krause" (default) — the published Moreton-simulation recipe (Krause et al. 2015):
#       GRAVITY-FREE isobaric corona (T=1.6 MK, p=2.65e-3 Pa, n~1.2e8 cm^-3) with a dense
#       chromospheric slab (T=1e4 K, density contrast ~160) below Z_SLAB. Robust at
#       dz>=200 km because the background is an exact static solution with no thin
#       hydrostatic transition region.
#   "valc"  — the full hydrostatic VAL-C stratification (photosphere/T-min/chromosphere/
#       TR/corona). Physically complete but the TR density scale length is ~0.2-0.3 Mm
#       for ANY resolvable T(z) path (hydrostatics fixes the contrast), so it requires
#       dz <~ 75 km to survive the blast impact — the expensive high-res variant.
const MODE = get(ENV,"MODE","krause")
const KR   = MODE=="krause"
getf(k,d)=parse(Float64,get(ENV,k,string(d)))
const T_COR_K = 1.6e6                      # K   Krause corona temperature
const P_COR_K = 2.65e-3                    # Pa  Krause coronal/slab pressure (isobaric)
const T_SLAB  = 1.0e4                      # K   Krause chromosphere-slab temperature
const Z_SLAB  = getf("Z_SLAB_MM", 5.0)*1e6 # m   slab top (Krause: 5 Mm)
const W_SLAB  = 1.0e6                      # m   slab-edge smoothstep half-width
const A_PULSE = getf("A_PULSE", KR ? 100.0 : 30.0) # peak Δp/p0 (Krause15 observable-Moreton
                                           #   threshold ~100; VC08 blast fiducial 20-40)
const B0_G    = getf("B0_G", 3.0)          # coronal field [G] (Krause15 nominal 3 G;
                                           #   West+11 quiet Sun ~3 G)
const θB      = deg2rad(getf("THETA_B_DEG", 0.0))  # 0 = vertical field
const Z_FL    = getf("Z_FL_MM", KR ? 35.0 : 10.0)*1e6  # pulse centre height (Krause15: 35 Mm)
const X_FL    = getf("X_FL_MM", KR ? 40.0 : 30.0)*1e6  # pulse centre (lateral)
const σ_FL    = getf("SIG_FL_MM", 5.0)*1e6 # lateral Gaussian σ (kernel <~10 Mm, VC08)
const T_FL_MAX = getf("T_FL_MAX_MK", 40.0)*1e6  # K  max flare-region temperature (Krause15
                                           # Sect. 2.1: T<=40 MK, Aschwanden 2005); above it the
                                           # PULSE IS MASS-LOADED instead of heated further
const σ_FZ    = getf("SIG_FZ_MM", KR ? 5.0 : 2.5)*1e6  # vertical σ (elliptical in valc mode:
                                           # keeps the in-situ tail off the TR, which in
                                           # PRESSURE rivals the low corona)
const τ_FL    = getf("TAU_FL_S", KR ? 50.0 : 60.0)     # impulsive phase (Krause: 40-50 s)
const Lx      = getf("LX_MM", KR ? 250.0 : 200.0)*1e6
const Lz      = getf("LZ_MM", KR ? 60.0 : 32.0)*1e6    # domain top (above zlo)
const zlo     = KR ? 0.0 : -0.6e6          # valc: shallow sub-photospheric anchor
const dxr     = getf("DX_KM", KR ? 250.0 : 200.0)*1e3
const dzr     = getf("DZ_KM", KR ? 250.0 : 100.0)*1e3
const TEND    = getf("TEND_S", 420.0)
const SNAP_DT = getf("SNAP_DT_S", 30.0)
const TRACE_DT= 1.0                        # cadence of the x-trace stacks [s]
const XSUB    = 2                          # x-subsampling of 2D snapshots
const Z_HA    = getf("Z_HA_MM", KR ? 4.0 : 1.2)*1e6   # "Hα proxy" trace height (upper slab /
                                           #   chromosphere — where the skirt swing is largest)
const Z_EIT   = getf("Z_EIT_MM", KR ? 15.0 : 8.0)*1e6 # coronal-front trace height
const SP_MM   = getf("SP_MM", 0.0)          # sponge width [Mm]; 0 => legacy 50-CELL width.
                                           # A cell-specified sponge changes PHYSICAL thickness
                                           # with resolution: at dz=500 km the 50-cell top layer
                                           # is 25 Mm and swallows the upper half of a plug at
                                           # z_fl=35 Mm. Specify Mm for resolution studies.
const SP_NC   = SP_MM > 0 ? max(4, round(Int, SP_MM*1e6/min(dxr,dzr))) : 50   # sponge width [cells]
const SP_RATE = 0.5                        # sponge damping rate [1/s]
const VMAX2   = (3.0e6)^2                  # velocity ceiling [m/s]² of the protection block
const BCZLO   = Symbol(get(ENV,"BCZLO","fixed"))   # :fixed (legacy) | :reflect (Krause15 floor)

# --- temperature / mean-molecular-weight profiles ------------------------------
"temperature vs height: gravity-free slab+corona (krause) or piecewise VAL-C (valc)"
function Tprof(z)
    if KR
        u = clamp((z-(Z_SLAB-W_SLAB))/(2W_SLAB), 0.0, 1.0)
        return T_SLAB + (T_COR_K-T_SLAB)*u*u*(3.0-2.0*u)
    end
    if z < 0.0
        return T_PHOT - GRAD_SUB*z                     # rises inward
    elseif z < Z_MINV
        return T_PHOT + (T_MINV-T_PHOT)*z/Z_MINV
    elseif z < Z_CHTOP
        return T_MINV + (T_CHROM-T_MINV)*(z-Z_MINV)/(Z_CHTOP-Z_MINV)
    else
        u = clamp((z-Z_CHTOP)/(Z_TRTOP-Z_CHTOP), 0.0, 1.0)
        return T_CHROM + (T_COR-T_CHROM)*u*u*(3.0-2.0*u)
    end
end
"mean molecular weight: fully ionized in krause mode (their slab is ionized plasma);
neutral (1.25) below ~10⁴ K -> ionized H+He (0.6) in valc mode"
μmw(T) = KR ? 0.6 : 0.6 + (1.25-0.6)*0.5*(1.0-tanh((T-2.0e4)/2.0e4))

"background column: krause = isobaric gravity-free (exact static state at any
resolution); valc = hydrostatic dp/dz = -p μ m_H g/(k_B T) on a 1 km fine grid"
function build_atmosphere()
    zf = collect(range(zlo, zlo+Lz+4dzr; step=1e3))    # 1 km fine grid, past the top ghosts
    if KR
        p = fill(P_COR_K, length(zf))
        ρ = [P_COR_K*μmw(Tprof(zz))*mH/(kB*Tprof(zz)) for zz in zf]
        return (zf, ρ, p)
    end
    j0 = findmin(abs.(zf))[2]                          # index of z≈0
    lnp = zeros(length(zf)); lnp[j0]=log(P_PHOT)
    invH(z) = μmw(Tprof(z))*mH*gsun/(kB*Tprof(z))
    for j in j0+1:length(zf);  lnp[j]=lnp[j-1]-0.5*(invH(zf[j-1])+invH(zf[j]))*(zf[j]-zf[j-1]); end
    for j in j0-1:-1:1;        lnp[j]=lnp[j+1]+0.5*(invH(zf[j+1])+invH(zf[j]))*(zf[j+1]-zf[j]); end
    p  = exp.(lnp)
    ρ  = [p[j]*μmw(Tprof(zf[j]))*mH/(kB*Tprof(zf[j])) for j in eachindex(zf)]
    (zf, ρ, p)
end
function interp(xq,x,y); xq≤x[1] && return y[1]; xq≥x[end] && return y[end]
    k=searchsortedlast(x,xq); w=(xq-x[k])/(x[k+1]-x[k]); (1-w)*y[k]+w*y[k+1]; end

# --- operator-split sponge: damp v and e_int toward the hydrostatic background,
#     leave ρ and B untouched (mass-conserving; same design as paper_run.jl) ----
@inline function sponge_cell!(U,a,γ,μ0,i,j,f)
    g1=1-f
    ρ=max(U[M.IRHO,i,j],1e-14)
    mx,my,mz=U[M.IMX,i,j],U[M.IMY,i,j],U[M.IMZ,i,j]
    pB=0.5*(U[M.IBX,i,j]^2+U[M.IBY,i,j]^2+U[M.IBZ,i,j]^2)/μ0
    ke=0.5*(mx^2+my^2+mz^2)/ρ
    eint=U[M.IE,i,j]-ke-pB
    U[M.IMX,i,j]=g1*mx; U[M.IMY,i,j]=g1*my; U[M.IMZ,i,j]=g1*mz
    U[M.IE,i,j]=(g1*eint+f*a.p0[j]/(γ-1))+g1^2*ke+pB
end
function sponges!(U,g,a,γ,μ0,dt)
    Nx,Nz=g.Nx,g.Nz
    @inbounds for j in 1:Nz+2M.NG, k in 1:SP_NC          # left & right x-layers
        ξ=(SP_NC-k+1)/SP_NC; f=SP_RATE*ξ^2*dt; f=f/(1+f)
        sponge_cell!(U,a,γ,μ0,M.NG+k,j,f); sponge_cell!(U,a,γ,μ0,M.NG+Nx-k+1,j,f)
    end
    @inbounds for k in 1:SP_NC, i in M.NG+1:M.NG+Nx      # top z-layer
        ξ=k/SP_NC; f=SP_RATE*ξ^2*dt; f=f/(1+f)
        sponge_cell!(U,a,γ,μ0,i,M.NG+Nz-SP_NC+k,f)
    end
end

function main()
    γ=γc; μ0=μ0c
    zfine,ρfine,pfine = build_atmosphere()
    Nx=round(Int,Lx/dxr); Nz=round(Int,Lz/dzr)
    g=M.Grid(Nx,Nz,Lx,Lz; z0=zlo); zpad(j)=zlo+(j-M.NG-0.5)*g.dz
    ρ0=[interp(zpad(j),zfine,ρfine) for j in 1:Nz+2M.NG]
    p0=[interp(zpad(j),zfine,pfine) for j in 1:Nz+2M.NG]
    B0=B0_G*1e-4; Bx0=B0*sin(θB); Bz0=B0*cos(θB)
    a=M.Atmos(KR ? 0.0 : gsun, ρ0,p0,Bx0,Bz0)   # krause mode is gravity-free (isobaric exact state)

    # --- background diagnostics + table -------------------------------------
    label=get(ENV,"RUN_LABEL",@sprintf("A%.0f_B%.0fG",A_PULSE,B0_G))
    outdir=joinpath(@__DIR__,"run_"*label); mkpath(outdir)
    open(joinpath(outdir,"background.txt"),"w") do io
        println(io,"# z[m]  rho[kg/m3]  p[Pa]  T[K]  c_s[m/s]  V_A[m/s]  V_fast_perp[m/s]")
        for j in M.NG+1:M.NG+Nz
            z=zpad(j); T=Tprof(z); cs=sqrt(γ*p0[j]/ρ0[j]); va=B0/sqrt(μ0*ρ0[j])
            @printf(io,"%.4e  %.6e  %.6e  %.1f  %.4e  %.4e  %.4e\n",
                    z,ρ0[j],p0[j],T,cs,va,sqrt(cs^2+va^2))
        end
    end
    jcor=findfirst(j->zpad(j)>10e6, M.NG+1:M.NG+Nz)+M.NG
    csc=sqrt(γ*p0[jcor]/ρ0[jcor]); vac=B0/sqrt(μ0*ρ0[jcor]); vfc=sqrt(csc^2+vac^2)
    ncor=ρ0[jcor]/(1.2*mH)   # rough total H density for the printout
    @printf("solar atmosphere: p_phot=%.2e Pa, corona @10Mm: rho=%.2e kg/m3 (n_H~%.1e cm^-3), c_s=%.0f km/s, V_A=%.0f km/s, V_f=%.0f km/s\n",
            P_PHOT,ρ0[jcor],ncor*1e-6,csc/1e3,vac/1e3,vfc/1e3)

    # --- initial state = static background -----------------------------------
    U=M.alloc(g); R=similar(U);U1=similar(U);W=similar(U);Fx=similar(U);Fz=similar(U)
    for j in 1:Nz+2M.NG, i in 1:Nx+2M.NG
        c=M.cons(γ,μ0,ρ0[j],0.,0.,0.,Bx0,0.,Bz0,p0[j],0.); for m in 1:M.NVAR; U[m,i,j]=c[m]; end
    end

    # --- driver: operator-split heating, Gaussian in (x,z). Two DRIVE modes:
    #   "pulse" (default) — impulsive flare blast: sin ramp over τ_FL, total central
    #       deposit = A·p0(z_fl)/(γ-1) (peak overpressure ≈ A). Cranmer/blast scenario.
    #   "piston" — SUSTAINED, laterally EXPANDING driver over τ_DRIVE: mimics the
    #       continuous lateral expansion of a CME flank (Chen 2002; Vršnak & Cliver 2008),
    #       which keeps pumping the front during propagation rather than launching it once.
    #       Steady deposition rate = A·p0/(γ-1)/τ_FL held over τ_DRIVE (cos ramps at both
    #       ends), with the driven lateral extent growing as σ_x(t)=σ_FL(1+V_pist·t/σ_FL)
    #       so the driving edge tracks outward with the front. This is the "extreme,
    #       continuously driven" companion to the impulsive run — the referee test for
    #       whether the short visibility range is driver-limited or dissipation-limited.
    p0fl=interp(Z_FL,zfine,pfine)
    zgate(z)= get(ENV,"ZGATE","0")=="0" ? 1.0 : 0.5*(1.0+tanh((z-4.0e6)/1.0e6))
    DRIVE   = get(ENV,"DRIVE","pulse")
    τ_DRIVE = getf("TAU_DRIVE_S", 180.0)           # sustained-driving duration [s] (CME accel. phase)
    V_PIST  = getf("V_PIST_KMS", 800.0)*1e3        # lateral expansion speed of the piston edge [m/s]
    q0     = A_PULSE*p0fl/(γ-1) * (π/(2τ_FL))      # ∫₀^τ q0 sin(πt/τ) dt = A·p0/(γ-1) (pulse)
    qdot_p = A_PULSE*p0fl/(γ-1) / τ_FL             # steady deposition rate (piston)
    # --- Krause et al. (2015) flare plug (DRIVE=krause_ic | krause_piston) -----------
    # Their Sect. 2.1 prescription, which our energy-deposition driver does NOT reproduce:
    # the pressure in the flare region is raised by a factor F=1+A*w at CONSTANT DENSITY,
    # "but, if for a given pressure increment the maximum temperature is exceeded, the
    # density is increased to maintain the temperature below the threshold". So their
    # Dp/p ~ 1e4 blast is a HOT, DENSE plug (T~40 MK, n~400 n_u), not a pure heating
    # deposit -- which is also why it stays numerically tame: mass loading holds V_A
    # (and hence the local fast speed) down. T/T_u = (p/p0)/(rho/rho0) needs no mu.
    R_T = T_FL_MAX/T_COR_K
    function krause_plug!(U)
        @inbounds for j in M.NG+1:M.NG+Nz
            gz=exp(-(zpad(j)-Z_FL)^2/(2σ_FZ^2))*zgate(zpad(j))
            gz<1e-6 && continue
            for i in M.NG+1:M.NG+Nx
                x=(i-M.NG-0.5)*g.dx
                w=gz*exp(-(x-X_FL)^2/(2σ_FL^2))
                w<1e-6 && continue
                F  = 1.0 + A_PULSE*w              # p/p0 in the plug
                TT = min(F, R_T)                  # T/T_u, capped at T_FL_MAX
                ρn = ρ0[j]*F/TT                   # density raised to hold T <= T_max
                ρc = max(U[M.IRHO,i,j],1e-14)
                c  = M.cons(γ,μ0, ρn,
                            U[M.IMX,i,j]/ρc, U[M.IMY,i,j]/ρc, U[M.IMZ,i,j]/ρc,
                            U[M.IBX,i,j], U[M.IBY,i,j], U[M.IBZ,i,j],
                            p0[j]*F, U[M.IPSI,i,j])
                for m in 1:M.NVAR; U[m,i,j]=c[m]; end
            end
        end
    end
    function flare_heat!(U,t,dt)
        if DRIVE=="krause_ic"
            return                                 # applied once, before the time loop
        elseif DRIVE=="krause_piston"
            t≤τ_DRIVE && krause_plug!(U)           # held for τ_DRIVE, then released
            return
        elseif DRIVE=="piston"
            t>τ_DRIVE && return
            ramp = 0.5*(1-cos(π*min(t,τ_FL)/τ_FL)) * 0.5*(1-cos(π*min(τ_DRIVE-t,τ_FL)/τ_FL))
            σx   = σ_FL + V_PIST*t
            amp  = qdot_p*ramp*dt
            @inbounds for j in M.NG+1:M.NG+Nz
                gz=exp(-(zpad(j)-Z_FL)^2/(2σ_FZ^2))*zgate(zpad(j))
                for i in M.NG+1:M.NG+Nx
                    x=(i-M.NG-0.5)*g.dx
                    U[M.IE,i,j]+=amp*gz*exp(-(x-X_FL)^2/(2σx^2))
                end
            end
        else
            t>τ_FL && return
            amp=q0*sin(π*t/τ_FL)*dt
            @inbounds for j in M.NG+1:M.NG+Nz
                gz=exp(-(zpad(j)-Z_FL)^2/(2σ_FZ^2))*zgate(zpad(j))
                for i in M.NG+1:M.NG+Nx
                    x=(i-M.NG-0.5)*g.dx
                    U[M.IE,i,j]+=amp*gz*exp(-(x-X_FL)^2/(2σ_FL^2))
                end
            end
        end
    end
    # post-pulse flare-site relaxer: once the front has left (t>3τ_fl), gently damp the
    # velocity in a Gaussian column around the flare site (z<10 Mm) so the cavity
    # refill / fall-back onto the TR cannot steepen without limit. The escaping wave is
    # ~V_f·t away by then; this only removes the flare-site aftermath, not the front.
    # post-pulse flare-site relaxer: once the front has detached, gently damp the flare
    # column so the site cannot keep ringing and pumping slosh into the domain for the
    # whole run (in ideal MHD the hot wake has no physical relaxation channel). The zone
    # covers the wake column: below ~10 Mm in valc mode; the slab-top-to-above-pulse
    # column in krause mode. Moreton kinematics are read at x-x_fl >~ 30 Mm.
    site_zrelax(z) = KR ? 0.5*(1.0-tanh((z-(Z_FL+2σ_FZ))/2.0e6)) :
                          0.5*(1.0-tanh((z-10.0e6)/1.0e6))
    t_wake = ((DRIVE=="piston"||DRIVE=="krause_piston") ? τ_DRIVE : τ_FL) + 10.0  # wait for the driver
    function site_relax!(U,t,dt)
        t<t_wake && return
        @inbounds for j in M.NG+1:M.NG+Nz
            gz=site_zrelax(zpad(j))
            gz<1e-3 && continue
            for i in M.NG+1:M.NG+Nx
                x=(i-M.NG-0.5)*g.dx
                f=0.10*dt*gz*exp(-(x-X_FL)^2/(2*(2σ_FL)^2)); f=f/(1+f)
                f<1e-6 && continue
                sponge_cell!(U,a,γ,μ0,i,j,f)
            end
        end
    end
    # sacrificial TR guard column: the blast's DIRECT overhead impact on the (numerically
    # thin) transition region has no relaxation channel in ideal MHD and evacuates TR cells
    # however the pulse is shaped; gently relax the TR layer in a narrow column under the
    # flare at all times. Moreton kinematics are read at x-x_fl >~ 30 Mm, far outside.
    zgd_c = KR ? Z_SLAB : 2.7e6          # guarded interface: slab edge / transition region
    zgd_w = KR ? 1.5e6  : 1.2e6
    function tr_guard!(U,dt)
        @inbounds for j in M.NG+1:M.NG+Nz
            z=zpad(j); (zgd_c-3zgd_w<z<zgd_c+3zgd_w) || continue
            gzz=exp(-((z-zgd_c)/zgd_w)^2)
            for i in M.NG+1:M.NG+Nx
                x=(i-M.NG-0.5)*g.dx
                f=0.15*dt*gzz*exp(-(x-X_FL)^2/(2*(1.5σ_FL)^2)); f=f/(1+f)
                f<1e-7 && continue
                sponge_cell!(U,a,γ,μ0,i,j,f)
            end
        end
    end

    # --- output machinery -----------------------------------------------------
    jha  = argmin(abs.([zpad(j) for j in M.NG+1:M.NG+Nz].-Z_HA ))+M.NG
    jeit = argmin(abs.([zpad(j) for j in M.NG+1:M.NG+Nz].-Z_EIT))+M.NG
    xs=collect(g.xc)./1e6                       # Mm
    writedlm(joinpath(outdir,"grid_x_Mm.txt"),xs[1:XSUB:end])
    writedlm(joinpath(outdir,"grid_z_Mm.txt"),[zpad(M.NG+j)/1e6 for j in 1:Nz])
    open(joinpath(outdir,"meta.txt"),"w") do io
        @printf(io,"bczlo=%s\nsp_nc=%d\nmode=%s\ndrive=%s\nA_pulse=%.1f\nB0_G=%.2f\ntheta_B_deg=%.1f\nz_fl_Mm=%.1f\nx_fl_Mm=%.1f\nsig_fl_Mm=%.1f\ntau_fl_s=%.0f\nLx_Mm=%.0f\nLz_Mm=%.0f\ndx_km=%.0f\ndz_km=%.0f\nNx=%d\nNz=%d\ntend_s=%.0f\ncf_cor_kms=%.1f\nz_ha_Mm=%.2f\nz_eit_Mm=%.2f\nxsub=%d\n",
                String(BCZLO),SP_NC,MODE,DRIVE,A_PULSE,B0_G,rad2deg(θB),Z_FL/1e6,X_FL/1e6,σ_FL/1e6,τ_FL,Lx/1e6,Lz/1e6,
                g.dx/1e3,g.dz/1e3,Nx,Nz,TEND,vfc/1e3,zpad(jha)/1e6,zpad(jeit)/1e6,XSUB)
    end
    trvz=Vector{Vector{Float64}}(); trdp=Vector{Vector{Float64}}(); trt=Float64[]
    snapt=Float64[]; si=0
    function save_snapshot(t)
        si+=1
        vz2d=[U[M.IMZ,M.NG+i,M.NG+j]/max(U[M.IRHO,M.NG+i,M.NG+j],1e-14)/1e3 for i in 1:XSUB:Nx, j in 1:Nz]
        vx2d=[U[M.IMX,M.NG+i,M.NG+j]/max(U[M.IRHO,M.NG+i,M.NG+j],1e-14)/1e3 for i in 1:XSUB:Nx, j in 1:Nz]
        dpp =[M.gaspressure(γ,μ0,U,M.NG+i,M.NG+j)/p0[M.NG+j]-1.0 for i in 1:XSUB:Nx, j in 1:Nz]
        writedlm(joinpath(outdir,@sprintf("snap_vz_%03d.txt",si)),vz2d)
        writedlm(joinpath(outdir,@sprintf("snap_vx_%03d.txt",si)),vx2d)
        writedlm(joinpath(outdir,@sprintf("snap_dpp_%03d.txt",si)),dpp)
        push!(snapt,t); writedlm(joinpath(outdir,"snap_times.txt"),snapt)
        @printf("   >> snapshot %03d  t=%5.0f s   peak|vz|=%.2f km/s  peak dp/p=%.2f\n",
                si,t,maximum(abs,vz2d),maximum(dpp)); flush(stdout)
    end

    # bounded-state protection (PLUTO/Athena-style failsafe): HLL+MUSCL has no positivity
    # guarantee, and a blast-driven rarefaction can empty a thin interface cell
    # (rho<0 -> v=m/rho_floor runaway that poisons its neighbours). In flagged cells:
    # rho floored at 1e-3*rho0 (physical wake evacuation is ~1e-1*rho0), |v| capped at
    # VMAX (physical max in this problem ~1500 km/s), eint clamped to [1e-3,100]*p0
    # (blast peak is A). Hit counts reported so silent intervention cannot pass unnoticed.
    # CRITICAL: applied to BOTH RK stages — the intermediate state U1 is where positivity
    # first fails, and stage-2 fluxes computed from a bad U1 poison the full step.
    nfloor=0
    function protect!(A)
        @inbounds for j in M.NG+1:M.NG+Nz, i in M.NG+1:M.NG+Nx
            rfl=1e-3*a.ρ0[j]
            ρc=A[M.IRHO,i,j]
            flag = ρc<rfl
            flag && (A[M.IRHO,i,j]=rfl; ρc=rfl)
            v2=(A[M.IMX,i,j]^2+A[M.IMY,i,j]^2+A[M.IMZ,i,j]^2)/ρc^2
            if flag || v2>VMAX2
                pB=0.5*(A[M.IBX,i,j]^2+A[M.IBY,i,j]^2+A[M.IBZ,i,j]^2)/μ0
                ke=0.5*ρc*v2
                eint=A[M.IE,i,j]-ke-pB
                if v2>VMAX2
                    sc=sqrt(VMAX2/v2)
                    A[M.IMX,i,j]*=sc; A[M.IMY,i,j]*=sc; A[M.IMZ,i,j]*=sc
                    ke*=VMAX2/v2
                end
                eint=clamp(eint, 1e-3*a.p0[j]/(γ-1), 100.0*a.p0[j]/(γ-1))
                A[M.IE,i,j]=eint+ke+pB; nfloor+=1
            end
        end
    end
    # local SSP-RK2 step = M.step_grav! with per-stage protection (solver core untouched)
    function step_protected!(dt,ch)
        M.rhs_grav!(R,U,W,Fx,Fz,g,γ,μ0,ch,a,:outflow; bczlo=BCZLO);  @. U1=U+dt*R
        protect!(U1)
        M.rhs_grav!(R,U1,W,Fx,Fz,g,γ,μ0,ch,a,:outflow; bczlo=BCZLO); @. U=0.5*U+0.5*(U1+dt*R)
        protect!(U)
        α=exp(-0.18*ch*dt/min(g.dx,g.dz))
        @inbounds for j in M.NG+1:M.NG+Nz, i in M.NG+1:M.NG+Nx; U[M.IPSI,i,j]*=α; end
    end

    @printf("Moreton run '%s' [%s]: %dx%d cells (dx=%.0f km, dz=%.0f km), A=%.0f, B0=%.1f G, tend=%.0f s\n",
            label,MODE,Nx,Nz,g.dx/1e3,g.dz/1e3,A_PULSE,B0_G,TEND); flush(stdout)
    (DRIVE=="krause_ic" || DRIVE=="krause_piston") && krause_plug!(U)   # t=0 plug
    t=0.0; ns=0; next_snap=SNAP_DT; next_tr=0.0
    while t<TEND
        ch=M.maxspeed(U,g,γ,μ0); dt=min(0.3*min(g.dx,g.dz)/ch, TEND-t)
        if ch>3e7   # true runaway despite protection: report and abort
            imn=0;jmn=0;rmn=Inf; imx=0;jmx=0;smx=0.0
            for j in M.NG+1:M.NG+Nz, i in M.NG+1:M.NG+Nx
                r=U[M.IRHO,i,j]/a.ρ0[j]; r<rmn && (rmn=r;imn=i;jmn=j)
                sp=abs(U[M.IMX,i,j]/max(U[M.IRHO,i,j],1e-14))+abs(U[M.IMZ,i,j]/max(U[M.IRHO,i,j],1e-14))
                sp>smx && (smx=sp;imx=i;jmx=j)
            end
            @printf("   ** runaway step %d t=%.2f ch=%.2e: min rho/rho0=%.2e at x=%.1f z=%.2f; max|v|=%.2e at x=%.1f z=%.2f **\n",
                    ns,t,ch,rmn,(imn-M.NG-0.5)*g.dx/1e6,zpad(jmn)/1e6,
                    smx,(imx-M.NG-0.5)*g.dx/1e6,zpad(jmx)/1e6); flush(stdout); break
        end
        step_protected!(dt,ch)
        flare_heat!(U,t,dt)
        site_relax!(U,t,dt)          # post-pulse wake damping (both modes)
        tr_guard!(U,dt)              # interface guard column under the flare (both modes)
        sponges!(U,g,a,γ,μ0,dt)
        t+=dt; ns+=1
        if t≥next_tr
            push!(trt,t)
            push!(trvz,[U[M.IMZ,M.NG+i,jha]/max(U[M.IRHO,M.NG+i,jha],1e-14)/1e3 for i in 1:Nx])
            push!(trdp,[M.gaspressure(γ,μ0,U,M.NG+i,jeit)/p0[jeit]-1.0 for i in 1:Nx])
            next_tr+=TRACE_DT
        end
        t≥next_snap && (save_snapshot(t); next_snap+=SNAP_DT)
        if !isfinite(sum(@view U[M.IE,M.NG+1:M.NG+Nx,M.NG+1:M.NG+Nz]))
            found=false
            for j in M.NG+1:M.NG+Nz, i in M.NG+1:M.NG+Nx
                if !all(m->isfinite(U[m,i,j]),1:M.NVAR)
                    @printf("   ** NaN at step %d t=%.2f dt=%.4f: first cell x=%.1f Mm z=%.2f Mm **\n",
                            ns,t,dt,(i-M.NG-0.5)*g.dx/1e6,zpad(j)/1e6); found=true; break
                end
            end
            found || println("   ** NaN in E-sum but no single non-finite cell?? **")
            flush(stdout); break
        end
        if ns%500==0
            @printf("   step %6d  t=%6.1f/%.0f s (%3.0f%%)  dt=%.3f s\n",ns,t,TEND,100t/TEND,dt)
            flush(stdout)
        end
    end
    save_snapshot(t)
    writedlm(joinpath(outdir,"trace_times.txt"),trt)
    writedlm(joinpath(outdir,"trace_chromo_vz.txt"),reduce(hcat,trvz)')   # rows=time, cols=x
    writedlm(joinpath(outdir,"trace_corona_dpp.txt"),reduce(hcat,trdp)')
    writedlm(joinpath(outdir,"trace_x_Mm.txt"),xs)
    @printf("done: %d steps to t=%.0f s, %d snapshots, %d floor-cell hits -> %s\n",ns,t,si,nfloor,outdir)
end
main()
