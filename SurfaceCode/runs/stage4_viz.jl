# stage4_viz.jl — run one MHD fast-mode pulse and dump Φ(x) + 2D v_x snapshots for plotting.
# (Everything inside main() — top-level loops in Julia are type-unstable and slow.)
include(joinpath(@__DIR__,"..","src","mhd2d.jl")); using .MHD2D; const M=MHD2D
using Printf, DelimitedFiles
using Base.Threads: @threads
const Gc=6.674e-11; const Msun=1.989e30; const Rsun=6.957e8
const Mstar=0.5*Msun; const Rstar=2.0*Rsun; const gsurf=Gc*Mstar/Rstar^2
const μ0c=4π*1e-7; const γc=5/3; const Bstar=0.1

read_bg(f)=(d=readdlm(f;comments=true,comment_char='#'); p=sortperm(d[:,1]); (d[p,1],d[p,3],d[p,5]))
function interp(xq,x,y); xq≤x[1] && return y[1]; xq≥x[end] && return y[end]
    k=searchsortedlast(x,xq); w=(xq-x[k])/(x[k+1]-x[k]); (1-w)*y[k]+w*y[k+1]; end

function main()
    γ=γc; μ0=μ0c
    Lxfrac = length(ARGS)≥1 ? parse(Float64,ARGS[1]) : 0.036   # x-domain in R⋆
    cph    = length(ARGS)≥2 ? parse(Int,ARGS[2])     : 12      # cells per H_p (12=res×1, 24=res×2)
    zb,ρb,Pb=read_bg(joinpath(@__DIR__,"..","data","background_35deg.txt"))
    ρ0s=interp(0.,zb,ρb); P0s=interp(0.,zb,Pb); Hp=P0s/(ρ0s*gsurf)
    Lx=Lxfrac*Rstar; dxt=1e5*12/cph                            # dx: 100 km (res×1) / 50 km (res×2)
    Nx=round(Int,Lx/dxt); zlo=-4Hp; Lz=7Hp; Nz=max(64,round(Int,Lz/(Hp/cph)))
    g=M.Grid(Nx,Nz,Lx,Lz; z0=zlo); zpad(j)=zlo+(j-M.NG-0.5)*g.dz
    ρ0=[max(interp(zpad(j),zb,ρb),1e-12) for j in 1:Nz+2M.NG]
    p0=[max(interp(zpad(j),zb,Pb),1e-6) for j in 1:Nz+2M.NG]
    θ=deg2rad(35.0); Bz0=Bstar*cos(θ); Bx0=0.5*Bstar*sin(θ); a=M.Atmos(gsurf,ρ0,p0,Bx0,Bz0)
    cf=sqrt(γ*P0s/ρ0s+(Bx0^2+Bz0^2)/(μ0*ρ0s))
    U=M.alloc(g); R=similar(U);U1=similar(U);W=similar(U);Fx=similar(U);Fz=similar(U)
    for j in 1:Nz+2M.NG,i in 1:Nx+2M.NG
        c=M.cons(γ,μ0,ρ0[j],0.,0.,0.,Bx0,0.,Bz0,p0[j],0.); for m in 1:M.NVAR; U[m,i,j]=c[m]; end; end
    M0=2.0; τ=150.0; zw(j)=exp(-(zpad(j)/(1.5Hp))^2)
    function injf(U,gg,aa,γ,μ0,t)
        (0≤t≤τ)||return; vx=M0*cf*sin(π*t/τ)
        @inbounds for j in 1:gg.Nz+2M.NG,k in 1:M.NG; il=M.NG+1-k; ρ=aa.ρ0[j]; ux=vx*zw(j)
            c=M.cons(γ,μ0,ρ,ux,0.,0.,aa.Bx0[j],0.,aa.Bz0,aa.p0[j],0.); for m in 1:M.NVAR; U[m,il,j]=c[m]; end; end
    end
    outdir=joinpath(@__DIR__,"..","..","output","plots","surfacecode"); mkpath(outdir)
    writedlm(joinpath(outdir,"viz_x.txt"), collect(g.xc)./Rstar)
    writedlm(joinpath(outdir,"viz_z.txt"), [zpad(M.NG+j)/Hp for j in 1:Nz])
    save2d(tag)=writedlm(joinpath(outdir,"viz_vx_$tag.txt"),
        [U[M.IMX,M.NG+i,M.NG+j]/max(U[M.IRHO,M.NG+i,M.NG+j],1e-14)/1e3 for i in 1:Nx, j in 1:Nz])
    tend=1.3*(Nx*g.dx)/cf+τ; snaps=[1.5*τ, 0.11*tend, 0.25*tend]; si=1  # trigger, then mid-decay (still visible)
    snaptimes=Float64[]
    fluxx=zeros(Nx); t=0.0; ns=0
    @printf("viz run: Nx=%d Nz=%d, Lx=%.3f Rstar, dx=%.0f km, cf=%.1f km/s\n",Nx,Nz,Lxfrac,g.dx/1e3,cf/1e3); flush(stdout)
    while t<tend
        ch=M.maxspeed(U,g,γ,μ0); dt=0.4*min(g.dx,g.dz)/ch
        M.step_grav!(U,R,U1,W,Fx,Fz,g,γ,μ0,dt,ch,a,:outflow; inject=injf,t=t); t+=dt; ns+=1
        ns%2000==0 && (@printf("   step %6d  t=%.0f/%.0f s (%2.0f%%)  dt=%.2fs\n",ns,t,tend,100t/tend,dt); flush(stdout))
        @threads for i in 1:Nx                         # threaded flux accumulation (no race: each i owns fluxx[i])
            fx=0.0
            @inbounds for j in M.NG+1:M.NG+Nz
                ρ,vx,vy,vz,Bx,By,Bz,p,ψ=M.prim(γ,μ0,U,M.NG+i,j)
                pB=0.5*(Bx^2+By^2+Bz^2)/μ0; E=p/(γ-1)+0.5ρ*(vx^2+vy^2+vz^2)+pB
                fx+=(E+p+pB)*vx-Bx*(vx*Bx+vy*By+vz*Bz)/μ0; end
            fluxx[i]+=fx*g.dz*dt; end
        if si≤length(snaps) && t≥snaps[si]; save2d(si); push!(snaptimes,t); si+=1; end
    end
    writedlm(joinpath(outdir,"viz_snaptimes.txt"), snaptimes)
    writedlm(joinpath(outdir,"viz_flux.txt"), hcat(collect(g.xc)./Rstar, fluxx))
    @printf("wrote viz data to %s  (Nx=%d Nz=%d, %d steps, cf=%.1f km/s)\n", outdir, Nx, Nz, ns, cf/1e3)
end
main()
