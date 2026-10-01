# run_shocktests.jl — dump the two solver-validation shock tubes for the appendix figures.
#   Sod (hydro, Patch2D): x, rho, u, p at t=0.2 (Nx=200).  Exact solution added in the plotter.
#   Brio-Wu (MHD, MHD2D):  x, rho, vx, p, By at t=0.1, at Nx=400 (fiducial) and Nx=1600 (reference).
include(joinpath(@__DIR__,"patch2d.jl")); using .Patch2D
include(joinpath(@__DIR__,"mhd2d.jl"));   using .MHD2D
using DelimitedFiles, Printf
const P=Patch2D; const M=MHD2D

outdir=joinpath(@__DIR__,"..","output","shocktests"); mkpath(outdir)

# ---- Sod (hydro) ----
xs,ρs,ps,us,ns,t = P.sod(Nx=200)
writedlm(joinpath(outdir,"sod.txt"), hcat(xs,ρs,us,ps))     # x, rho, u, p
@printf("Sod:    Nx=%d, %d steps to t=%.3f\n", length(xs), ns, t)

# ---- Brio-Wu (MHD): full fields at two resolutions ----
function briowu_full(Nx; tend=0.1, cfl=0.4)
    γ=2.0; μ0=1.0; Nz=4; g=M.Grid(Nx,Nz,1.0,0.01)
    U=M.alloc(g)
    for j in 1:Nz+2M.NG, i in 1:Nx+2M.NG
        x=(i-M.NG-0.5)*g.dx
        ρ,p,By = x<0.5 ? (1.0,1.0,1.0) : (0.125,0.1,-1.0)
        c=M.cons(γ,μ0, ρ,0.0,0.0,0.0, 0.75,By,0.0, p,0.0)
        for m in 1:M.NVAR; U[m,i,j]=c[m]; end
    end
    R=similar(U);U1=similar(U);W=similar(U);Fx=similar(U);Fz=similar(U); t=0.0; n=0
    while t<tend
        ch=M.maxspeed(U,g,γ,μ0); dt=min(cfl*min(g.dx,g.dz)/ch, tend-t)
        M.step!(U,R,U1,W,Fx,Fz,g,γ,μ0,dt,ch,:outflow,:outflow); t+=dt; n+=1
    end
    jm=M.NG+2; xs=Float64[]; ρa=Float64[]; va=Float64[]; pa=Float64[]; Bya=Float64[]
    for i in 1:Nx
        ρ,vx,vy,vz,Bx,By,Bz,p,ψ = M.prim(γ,μ0,U,M.NG+i,jm)
        push!(xs,g.xc[i]); push!(ρa,ρ); push!(va,vx); push!(pa,p); push!(Bya,By)
    end
    (xs,ρa,va,pa,Bya,n,t)
end
x4,ρ4,v4,p4,By4,n4,t4 = briowu_full(400)
xR,ρR,vR,pR,ByR,nR,tR = briowu_full(1600)
writedlm(joinpath(outdir,"briowu_400.txt"), hcat(x4,ρ4,v4,p4,By4))
writedlm(joinpath(outdir,"briowu_ref.txt"), hcat(xR,ρR,vR,pR,ByR))
@printf("Brio-Wu: Nx=400 (%d steps) + Nx=1600 reference (%d steps), t=%.3f\n", n4, nR, t4)
println("wrote -> ", outdir)
