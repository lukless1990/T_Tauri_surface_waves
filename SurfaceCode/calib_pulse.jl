# calib_pulse.jl — find the driving amplitude M0 that launches a v_x ≈ 10 km/s pulse INTO the domain.
# The boundary drive is M0*cf; after transmission through the ghost interface + vertical dispersion the
# domain peak is smaller (M0=2 gave ~6.3 km/s). Sweep M0, report the domain-peak |vx| (a) at the trigger
# snapshot t=1.5τ (the definition used for the paper's injected amplitude) and (b) the global max over the
# launch phase. Small domain / short time — cheap. Usage: julia -t 14 SurfaceCode/calib_pulse.jl [M0 M0 ...]
include(joinpath(@__DIR__,"mhd2d.jl")); using .MHD2D; const M=MHD2D
using Printf, DelimitedFiles
const Gc=6.674e-11; const Msun=1.989e30; const Rsun=6.957e8
const Mstar=0.5*Msun; const Rstar=2.0*Rsun; const gsurf=Gc*Mstar/Rstar^2
const μ0c=4π*1e-7; const γc=5/3
const Bstar = haskey(ENV,"BSTAR_T") ? parse(Float64,ENV["BSTAR_T"]) : 0.1   # T (0.1=1kG default, 0.01=0.1kG)

read_bg(f)=(d=readdlm(f;comments=true,comment_char='#'); p=sortperm(d[:,1]); (d[p,1],d[p,3],d[p,5]))
function interp(xq,x,y); xq≤x[1] && return y[1]; xq≥x[end] && return y[end]
    k=searchsortedlast(x,xq); w=(xq-x[k])/(x[k+1]-x[k]); (1-w)*y[k]+w*y[k+1]; end

const TAU_S = haskey(ENV,"TAU_S") ? parse(Float64,ENV["TAU_S"]) : 150.0   # driver duration [s]

function peak_for(M0; Lxfrac=0.03, cph=24, τ=TAU_S, tend=800.0)
    γ=γc; μ0=μ0c
    zb,ρb,Pb=read_bg(joinpath(@__DIR__,"background_35deg.txt"))
    ρ0s=interp(0.,zb,ρb); P0s=interp(0.,zb,Pb); Hp=P0s/(ρ0s*gsurf)
    Lx=Lxfrac*Rstar; dxt=1e5*12/cph
    Nx=round(Int,Lx/dxt); zlo=-4Hp; Lz=7Hp; Nz=max(64,round(Int,Lz/(Hp/cph)))
    g=M.Grid(Nx,Nz,Lx,Lz; z0=zlo); zpad(j)=zlo+(j-M.NG-0.5)*g.dz
    ρ0=[max(interp(zpad(j),zb,ρb),1e-12) for j in 1:Nz+2M.NG]
    p0=[max(interp(zpad(j),zb,Pb),1e-6) for j in 1:Nz+2M.NG]
    θ=deg2rad(35.0); Bz0=Bstar*cos(θ); Bx0=0.5*Bstar*sin(θ); a=M.Atmos(gsurf,ρ0,p0,Bx0,Bz0)
    cf=sqrt(γ*P0s/ρ0s+(Bx0^2+Bz0^2)/(μ0*ρ0s))
    U=M.alloc(g); R=similar(U);U1=similar(U);W=similar(U);Fx=similar(U);Fz=similar(U)
    for j in 1:Nz+2M.NG,i in 1:Nx+2M.NG
        c=M.cons(γ,μ0,ρ0[j],0.,0.,0.,Bx0,0.,Bz0,p0[j],0.); for m in 1:M.NVAR; U[m,i,j]=c[m]; end; end
    zw(j)=exp(-(zpad(j)/(1.5Hp))^2)
    function injf(U,gg,aa,γ,μ0,t)
        (0≤t≤τ)||return; vx=M0*cf*sin(π*t/τ)
        @inbounds for j in 1:gg.Nz+2M.NG,k in 1:M.NG; il=M.NG+1-k; ρ=aa.ρ0[j]; ux=vx*zw(j)
            c=M.cons(γ,μ0,ρ,ux,0.,0.,aa.Bx0[j],0.,aa.Bz0,aa.p0[j],0.); for m in 1:M.NVAR; U[m,il,j]=c[m]; end; end
    end
    domainmax()=begin mx=0.0
        @inbounds for j in M.NG+1:M.NG+Nz,i in M.NG+1:M.NG+Nx
            mx=max(mx,abs(U[M.IMX,i,j]/max(U[M.IRHO,i,j],1e-14))); end; mx; end
    t=0.0; gmax=0.0; trigmax=-1.0; ttrig=1.5τ; got=false
    while t<tend
        ch=M.maxspeed(U,g,γ,μ0); dt=0.4*min(g.dx,g.dz)/ch
        M.step_grav!(U,R,U1,W,Fx,Fz,g,γ,μ0,dt,ch,a,:outflow; inject=injf,t=t); t+=dt
        dm=domainmax(); gmax=max(gmax,dm)
        if !got && t≥ttrig; trigmax=dm; got=true; end
        isfinite(gmax)||(return (cf,NaN,NaN))
    end
    (cf, trigmax, gmax)
end

M0list = isempty(ARGS) ? [2.0,2.8,3.2,3.6,4.0] : [parse(Float64,a) for a in ARGS]
@printf("calibrating injected v_x peak vs M0  (target domain peak = 10 km/s, tau = %.0f s)\n", TAU_S)
@printf("%6s  %12s  %14s  %12s\n","M0","boundary","peak@1.5τ","global max"); flush(stdout)
for M0 in M0list
    cf,trig,gm=peak_for(M0)
    @printf("%6.2f  %9.1f km/s  %10.2f km/s  %8.2f km/s\n", M0, M0*cf/1e3, trig/1e3, gm/1e3); flush(stdout)
end
