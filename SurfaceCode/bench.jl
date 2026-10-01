# quick allocation + timing probe for the MHD step (run: julia -t N SurfaceCode/bench.jl)
include(joinpath(@__DIR__,"mhd2d.jl")); using .MHD2D; const M=MHD2D; using Printf
γ=5/3; μ0=1.0; N=128; g=M.Grid(N,N,1.0,4.0)
zpad(j)=(j-M.NG-0.5)*g.dz
ρ0=[exp(-zpad(j)) for j in 1:N+2M.NG]; p0=[exp(-zpad(j)) for j in 1:N+2M.NG]
a=M.Atmos(1.0,ρ0,p0,0.3,1.0)
U=M.alloc(g)
for j in 1:N+2M.NG, i in 1:N+2M.NG
    c=M.cons(γ,μ0,ρ0[j],0.,0.,0.,0.3,0.,1.0,p0[j],0.); for m in 1:M.NVAR; U[m,i,j]=c[m]; end
end
R=similar(U);U1=similar(U);W=similar(U);Fx=similar(U);Fz=similar(U)
ch=M.maxspeed(U,g,γ,μ0); dt=0.1*min(g.dx,g.dz)/ch
M.step_grav!(U,R,U1,W,Fx,Fz,g,γ,μ0,dt,ch,a,:periodic)      # warmup/compile
nb=@allocated for _ in 1:20; M.step_grav!(U,R,U1,W,Fx,Fz,g,γ,μ0,dt,ch,a,:periodic); end
t=@elapsed for _ in 1:100; M.step_grav!(U,R,U1,W,Fx,Fz,g,γ,μ0,dt,ch,a,:periodic); end
@printf("grid %d² | %d threads | %.2f MB/20steps (%.1f KB/step) | %.2f ms/step | %.1f Mcell-updates/s\n",
        N, Threads.nthreads(), nb/1e6, nb/20/1e3, t/100*1e3, N*N*2/(t/100)/1e6)
