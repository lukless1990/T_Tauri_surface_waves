# test1_stratified_wave.jl — SurfaceCode verification test 1 (paper appendix):
# linear wave in a gravitationally stratified atmosphere vs the EXACT analytic solution.
#
# Single-fluid ideal-MHD limit of the stratified-atmosphere test of
# Popescu Braileanu et al. 2019, A&A 627, A25, Sect. 5 (Mancha3D-2F verification):
#   isothermal layer ρ0,p0 ∝ e^{-z/H}; horizontal field Bx0(z)=B00 e^{-z/2H} (same
#   magnetic and gas scale height ⇒ c_s AND V_A constant with height);
#   magneto-hydrostatic balance (p00 + B00²/2μ0)/H = ρ00·g.
#   Vertically propagating fast wave — EXACT solution (no WKB; u=e^{z/2H}F reduces the
#   wave equation ∂²u/∂t² = a²u'' − (a²/H)u' to F̈ = a²F'' − (a²/4H²)F):
#       u_z(z,t) = V e^{z/2H} cos(ωt−kz),   ω² = a²(k² + 1/4H²),  a² = c_s²+V_A²
#   with polarization (convention e^{i(ωt−kz)}, physical = Re):
#       Bx1 = (k B00 V/ω) cos φ                                    (flat amplitude)
#       ρ1  = ρ00 e^{-z/2H} V [ (k/ω) cos φ + (1/2Hω) sin φ ]
#       p1  = p00 e^{-z/2H} V [ (γk/ω) cos φ + ((1−γ/2)/Hω) sin φ ]
# Driven at the bottom ghost cells (full polarization, smoothstep ramp over 4 periods),
# absorbed at the top by the production-style sponge; compared in the stationary state.
#
# What this verifies (none of Sod / Brio–Wu / static well-balancing touch it):
#   * wave propagation ON the deviation-reconstructed stratified background,
#     including the new z-dependent-Bx0 magneto-hydrostatic path in rhs_grav!,
#   * the e^{z/2H} amplitude growth (wave-energy conservation in stratification),
#   * the absorbing sponge (a standing-wave ripple in the error = reflection).
#
# Amplitude: the ladder drives V=1e-4 c_s. (The paper drives 1e-3 c0 but takes its
# <2% against a LINEARIZED twin run; our code is fully nonlinear, and at 1e-3 the
# physical second-order steepening — Fubini 2k harmonic, σ/2 with the e^{z/2H}
# growth — floors ε at ~3%, independent of resolution. At 1e-4 that floor drops
# to ~0.3% and the comparison is truncation-limited, as intended.)
# A separate 1e-3 run then MEASURES the 2k harmonic and checks it against the
# Fubini prediction — the paper's Fig.-7 nonlinear demonstration, and a direct
# check of the steepening mechanism the science result rests on.
#
# PASS: machine-precision static balance of the new background; ε(u_z) (normalized
# error, their Eq. 48; window z∈[0.5,1.75]H) decreasing monotonically along the
# resolution ladder and < 2% at the finest rung. The raw pointwise ε of the
# production TVD scheme is limiter-distortion-dominated (first order at wave
# extrema) — the 2% criterion is met at high resolution, while at production-like
# resolution the physics-relevant fidelity is reported per wavelength band:
# carrier amplitude ratio (its deficit matches the independent test-2 dissipation
# floor) and phase error.
# The dampB=false rerun quantifies the reflection of the PRODUCTION sponge
# (which damps v+eint only, leaving B untouched) — reported, not gated.
#
# Usage: julia -t 4 SurfaceCode/test1_stratified_wave.jl
include(joinpath(@__DIR__,"mhd2d.jl")); using .MHD2D; const M=MHD2D
using Printf, DelimitedFiles

# Top absorbing layer, production-style (velocity+eint damped toward the background,
# ρ untouched — mass-conserving). dampB additionally relaxes B → background: for this
# x-uniform test that is ∇·B-safe (∂x(f·Bx)=0) and removes the frozen Bx1 residue the
# v-only sponge leaves behind (the paper uses a PML for the same purpose).
function top_sponge!(U,g,a,γ,μ0,nsp,σ0,dt; dampB::Bool=true)
    @inbounds for j in M.NG+g.Nz-nsp+1:M.NG+g.Nz
        ξ=(j-(M.NG+g.Nz-nsp))/nsp
        f=σ0*ξ^2*dt; f=f/(1+f); g1=1-f
        eint_bg=a.p0[j]/(γ-1)
        for i in 1:g.Nx+2M.NG
            ρ=max(U[M.IRHO,i,j],1e-14)
            mx,my,mz=U[M.IMX,i,j],U[M.IMY,i,j],U[M.IMZ,i,j]
            pBold=0.5*(U[M.IBX,i,j]^2+U[M.IBY,i,j]^2+U[M.IBZ,i,j]^2)/μ0
            ke=0.5*(mx^2+my^2+mz^2)/ρ
            eint=U[M.IE,i,j]-ke-pBold
            U[M.IMX,i,j]=g1*mx; U[M.IMY,i,j]=g1*my; U[M.IMZ,i,j]=g1*mz
            if dampB
                U[M.IBX,i,j]=g1*U[M.IBX,i,j]+f*a.Bx0[j]
                U[M.IBY,i,j]*=g1
                U[M.IBZ,i,j]=g1*U[M.IBZ,i,j]+f*a.Bz0
            end
            pBnew=0.5*(U[M.IBX,i,j]^2+U[M.IBY,i,j]^2+U[M.IBZ,i,j]^2)/μ0
            U[M.IE,i,j]=(g1*eint+f*eint_bg)+g1^2*ke+pBnew
        end
    end
end

function run_case(; B00=1.0, cpl=64, amp=1e-4, dampB=true, dump=nothing)
    γ=5/3; μ0=1.0; H=1.0; ρ00=1.0; p00=1.0
    pB00=B00^2/(2μ0); grav=(p00+pB00)/(ρ00*H)          # MHS: total-pressure scale height = H
    cs2=γ*p00/ρ00; vA2=B00^2/(μ0*ρ00); a2=cs2+vA2
    λ=H/4; kw=2π/λ; ω=sqrt(a2*(kw^2+1/(4H^2))); Pw=2π/ω
    V=amp*sqrt(cs2)                                    # driving amplitude (see header)
    Lz=3.5H; Nz=round(Int,Lz*cpl/λ); dz=Lz/Nz; Nx=4; Lx=Nx*dz
    g=M.Grid(Nx,Nz,Lx,Lz; z0=0.0); zpad(j)=(j-M.NG-0.5)*dz
    ρ0=[ρ00*exp(-zpad(j)/H)   for j in 1:Nz+2M.NG]
    p0=[p00*exp(-zpad(j)/H)   for j in 1:Nz+2M.NG]
    Bx0=[B00*exp(-zpad(j)/2H) for j in 1:Nz+2M.NG]
    atm=M.Atmos(grav,ρ0,p0,Bx0,0.0)
    U=M.alloc(g); R=similar(U);U1=similar(U);W=similar(U);Fx=similar(U);Fz=similar(U)
    setbg!() = (for j in 1:Nz+2M.NG, i in 1:Nx+2M.NG
        c=M.cons(γ,μ0,ρ0[j],0.,0.,0.,Bx0[j],0.,0.,p0[j],0.)
        for m in 1:M.NVAR; U[m,i,j]=c[m]; end
    end)

    # -- static pre-check: the NEW z-dependent-Bx0 background must be held to machine precision --
    setbg!(); vpk=0.0
    for n in 1:400
        ch=M.maxspeed(U,g,γ,μ0); dt=0.4*min(g.dx,g.dz)/ch
        M.step_grav!(U,R,U1,W,Fx,Fz,g,γ,μ0,dt,ch,atm,:periodic)
        @inbounds for j in M.NG+1:M.NG+Nz, i in M.NG+1:M.NG+Nx
            vpk=max(vpk, abs(U[M.IMZ,i,j]/U[M.IRHO,i,j]), abs(U[M.IMX,i,j]/U[M.IRHO,i,j]))
        end
    end
    drift=0.0
    @inbounds for j in M.NG+1:M.NG+Nz, i in M.NG+1:M.NG+Nx
        drift=max(drift, abs(U[M.IRHO,i,j]-ρ0[j])/ρ0[j])
    end

    # -- driver: full-polarization analytic state in the bottom ghosts, smoothstep ramp --
    tramp=4Pw
    function inj!(U,gg,aa,γi,μ0i,t)
        s=clamp(t/tramp,0.0,1.0); s=s*s*(3-2s)
        @inbounds for kk in 1:M.NG
            jg=M.NG+1-kk; z=zpad(jg); φ=ω*t-kw*z
            ez2=exp(z/2H); emz2=exp(-z/2H); emz=exp(-z/H)
            uz  = s*V*ez2*cos(φ)
            Bx1 = s*(kw*B00*V/ω)*cos(φ)
            ρ1  = s*ρ00*emz2*V*((kw/ω)*cos(φ)+(1/(2H*ω))*sin(φ))
            p1  = s*p00*emz2*V*((γi*kw/ω)*cos(φ)+((1-γi/2)/(H*ω))*sin(φ))
            c=M.cons(γi,μ0i, ρ00*emz+ρ1, 0.0,0.0,uz, B00*emz2+Bx1,0.0,0.0, p00*emz+p1, 0.0)
            for i in 1:gg.Nx+2M.NG, m in 1:M.NVAR; U[m,i,jg]=c[m]; end
        end
    end

    # -- driven run to the stationary state --
    setbg!()
    tcross=Lz/sqrt(a2); tend=3.2*tcross
    nsp=round(Int,1.5H/dz); σ0=10*sqrt(a2)/H           # sponge: top 1.5H, ~e⁻³ per pass
    t=0.0; ns=0
    while t<tend
        ch=M.maxspeed(U,g,γ,μ0); dt=min(0.4*min(g.dx,g.dz)/ch, tend-t)
        M.step_grav!(U,R,U1,W,Fx,Fz,g,γ,μ0,dt,ch,atm,:periodic; inject=inj!,t=t)
        top_sponge!(U,g,atm,γ,μ0,nsp,σ0,dt; dampB=dampB)
        t+=dt; ns+=1
    end

    # -- normalized error (their Eq. 48) in the window, against the exact solution --
    i0=M.NG+1; num=0.0; den=0.0; numB=0.0; denB=0.0; xspread=0.0
    rows=Tuple{Float64,Float64,Float64,Float64,Float64}[]
    @inbounds for j in M.NG+1:M.NG+Nz
        z=zpad(j); (0.5H≤z≤1.75H) || continue
        φ=ω*t-kw*z
        ua=V*exp(z/2H)*cos(φ);        un=U[M.IMZ,i0,j]/U[M.IRHO,i0,j]
        Ba=(kw*B00*V/ω)*cos(φ);       Bn=U[M.IBX,i0,j]-Bx0[j]
        num+=(un-ua)^2; den+=ua^2
        if B00>0; numB+=(Bn-Ba)^2; denB+=Ba^2; end
        umin=Inf; umax=-Inf
        for i in M.NG+1:M.NG+Nx
            u=U[M.IMZ,i,j]/U[M.IRHO,i,j]; umin=min(umin,u); umax=max(umax,u)
        end
        xspread=max(xspread,umax-umin)
        push!(rows,(z,un,ua,Bn,Ba))
    end
    εu=sqrt(num/den); εB= B00>0 ? sqrt(numB/denB) : NaN
    # per-band carrier fidelity + 2k harmonic: project u/(V e^{z/2H}) on the exact phase
    rmin=Inf; rmax=-Inf; dφmax=0.0
    bandz=Float64[]; bandr=Float64[]; bandA2=Float64[]
    let per=round(Int,λ/dz)
        nb=length(rows)÷per
        for b in 0:nb-1
            Ac=0.0; As=0.0; A2c=0.0; A2s=0.0; zm=0.0
            for m in 1:per
                z,un,_,_,_=rows[b*per+m]; φ=ω*t-kw*z; gn=un/(V*exp(z/2H)); zm+=z/per
                Ac+=2gn*cos(φ)/per; As+=2gn*sin(φ)/per
                A2c+=2gn*cos(2φ)/per; A2s+=2gn*sin(2φ)/per
            end
            rb=hypot(Ac,As); δ=atan(-As,Ac)
            rmin=min(rmin,rb); rmax=max(rmax,rb); dφmax=max(dφmax,abs(δ))
            push!(bandz,zm); push!(bandr,rb); push!(bandA2,hypot(A2c,A2s))
        end
    end
    if dump!==nothing
        writedlm(dump, [collect(r) for r in rows])     # z, uz_num, uz_exact, Bx1_num, Bx1_exact
    end
    (εu=εu, εB=εB, drift=drift, vpk_static=vpk/sqrt(cs2), xspread=xspread/V, ns=ns, Nz=Nz,
     rmin=rmin, rmax=rmax, dφmax=dφmax, bandz=bandz, bandr=bandr, bandA2=bandA2, V=V, H=H, kw=kw, a=sqrt(a2), γ=γ)
end

outdir=joinpath(@__DIR__,"..","output","wavetests"); mkpath(outdir)
println("TEST 1 — stratified-atmosphere wave vs exact solution")
println("  (single-fluid limit of Popescu Braileanu et al. 2019 Sect. 5)")
println("  PASS: static balance ~1e-16; eps(u_z) monotone-decreasing along the ladder, <2% at the finest rung\n")
allpass=true
conv=[]                                                 # (B00, cpl, eps_u, eps_B, rmin, rmax, dphmax)
for (label,B00) in (("hydro (B=0, acoustic-gravity)",0.0), ("MHD (β=2, fast mode)",1.0))
    @printf("%s:\n", label)
    prev=Inf
    for cpl in (32,64,128,256,512)
        dump = cpl==512 ? joinpath(outdir, B00==0 ? "test1_hydro.txt" : "test1_mhd.txt") : nothing
        r=run_case(B00=B00, cpl=cpl, dampB=true, dump=dump)
        ok = r.εu<prev && r.drift<1e-10 && r.vpk_static<1e-10 && (cpl<512 || r.εu<0.02)
        global allpass &= ok
        @printf("  cpl=%3d (Nz=%5d): eps(u_z)=%.4f  eps(Bx1)=%s  band amp=[%.3f,%.3f]  |dphase|<%.3f  static drift=%.1e  [%s]\n",
                cpl, r.Nz, r.εu, isnan(r.εB) ? "  --  " : @sprintf("%.4f",r.εB),
                r.rmin, r.rmax, r.dφmax, r.drift, ok ? "PASS" : "FAIL")
        flush(stdout)
        push!(conv,(B00,cpl,r.εu,r.εB,r.rmin,r.rmax,r.dφmax))
        prev=r.εu
    end
end
writedlm(joinpath(outdir,"test1_convergence.txt"), [collect(c) for c in conv])
# informational: reflection of the PRODUCTION sponge (v+eint only, B untouched)
r=run_case(B00=1.0, cpl=128, dampB=false)
@printf("\nproduction sponge (no B damping), MHD cpl=128: eps(u_z)=%.4f  (excess over the dampB run ≈ sponge reflection)\n", r.εu)

# nonlinear steepening check (paper's Fig.-7 analog): 1e-3 c_s drive, hydro, 512 c/λ.
# Fubini: 2nd-harmonic/carrier = σ/2, σ(z) = (γ+1)/2 · (V/c_s) · k · 2H(e^{z/2H}−1).
rn=run_case(B00=0.0, cpl=512, amp=1e-3, dump=joinpath(outdir,"test1_nonlinear.txt"))
println("\nsteepening check (V=1e-3 c_s, hydro, 512 c/λ): measured 2k harmonic vs Fubini σ/2:")
for (zb,rb,A2) in zip(rn.bandz,rn.bandr,rn.bandA2)
    σ=(rn.γ+1)/2*(rn.V/rn.a)*rn.kw*2rn.H*(exp(zb/2rn.H)-1)
    @printf("  z=%.2fH: measured %.4f   Fubini %.4f\n", zb, A2/rb, σ/2)
end
println(allpass ? "\nTEST 1 PASS: stratified wave converges to the exact solution (<2% at the finest rung)." :
                  "\nTEST 1 FAIL — inspect.")
