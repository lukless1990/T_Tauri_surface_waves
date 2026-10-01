# ============================================================================
#  export_background.jl  —  SurfaceCode pilot, task 1 (see README.md §11, §13)
#
#  Dump the near-surface background stratification of the Stellar2D fiducial
#  T Tauri star at the accretion-ring latitude, for ingestion by the explicit
#  surface-wave patch code. Columns: height, ρ, T, P, ionization fraction,
#  neutral mass fraction, n_e, dipole B(r,θ_ring), V_A, c_s, H_p.
#
#  The sub-photospheric part (τ≳1) is the REAL mapped MLT envelope. Above the
#  τ=1 boundary the Stellar2D domain stops, so the atmosphere is extended
#  upward with an isothermal-hydrostatic column at the photospheric T — a
#  PLACEHOLDER chromosphere for the pilot (flagged in the output header;
#  replace with a proper chromospheric model for production).
#
#  Usage:
#    julia --project=. SurfaceCode/export_background.jl [theta_ring_deg] [n_Hp_up]
#      theta_ring_deg : ring colatitude (default 35)
#      n_Hp_up        : scale heights of chromosphere to append above τ=1 (default 10)
#  Writes: SurfaceCode/background_<theta>deg.txt
# ============================================================================

using Stellar2D, JLD2, Printf
using Stellar2D: Constants, SahaEOS, saha_state, rho_from_PT

const mH = Constants.m_H
const kB = Constants.k_B
const μ0 = Constants.μ0
const G  = Constants.G

θring_deg = length(ARGS) ≥ 1 ? parse(Float64, ARGS[1]) : 35.0
n_Hp_up   = length(ARGS) ≥ 2 ? parse(Float64, ARGS[2]) : 10.0
θ = deg2rad(θring_deg)

# Metal electron floor. The H+He Saha EOS omits metals; at T_eff~4000 K the real electron supply is
# dominated by easily-ionized alkalis/Ca (x_e ~ 1e-4), NOT by H (which gives x_ion ~ 1e-7 here). Since
# ambipolar damping ∝ 1/n_i, using the H+He value ALONE overestimates the damping — the conservative
# direction for a "the wave damps" thesis, so it must be corrected. We add x_metal electrons per H
# nucleus above T_metal (metals recombine below that). Set STELLAR2D_XMETAL=0 to recover raw H+He Saha.
x_metal = parse(Float64, get(ENV, "STELLAR2D_XMETAL", "1e-4"))
T_metal = 2500.0   # K; below this the alkalis largely recombine / molecules form

# --- load the cached envelope + EOS (built by the run scripts) --------------
cachef = joinpath(@__DIR__, "..", "output", "eosenv_cache.jld2")
isfile(cachef) || error("missing $cachef — run a Stellar2D run once to build it, or point to your cache")
@load cachef eos env ROUT
Rstar = 2.0 * Constants.R_sun
Mstar = 0.5 * Constants.M_sun
Bstar = 0.1            # Tesla = 1000 G, matches the runs / Cranmer's base field
GM    = G * Mstar

# Saha EOS for the ionization state (env/eos in the cache is the TABULATED eos;
# saha_state needs the underlying SahaEOS — rebuild it with the run composition).
p_saha = Params(; Mstar=Mstar, Rstar=Rstar)
sah = SahaEOS(p_saha)
Xh = sah.X

# --- dipole field at the ring latitude (initial2.jl seed_dipole!) -----------
# B_r = B⋆ (R/r)³ cosθ · ... ; at the surface r=R: B_r=B⋆cosθ, B_θ=½B⋆sinθ.
kdip = Bstar * Rstar^3 / 2
Br(r) = 2kdip * cos(θ) / r^3
Bθ(r) = kdip * sin(θ) / r^3

# --- ionization / neutral fraction from Saha number densities ---------------
function ion_state(ρ, T)
    nH2, nHI, nHII, nHeI, nHeII, nHeIII, ne = saha_state(sah, ρ, T)
    nH = Xh * ρ / mH                                   # total H nuclei
    ne_metal = (T > T_metal) ? x_metal * nH : 0.0      # alkali/Ca electrons (not in H+He Saha)
    ne_eff = ne + ne_metal                             # total free electrons
    x_ion = ne_eff / max(nH, eps())                    # electrons per H nucleus
    ρn = (nHI + 2nH2) * mH + nHeI * 4mH                # neutral mass density (metals trace ⇒ ignored)
    ρi = nHII * mH + (nHeII + nHeIII) * 4mH + ne_metal*mH
    ξn = ρn / max(ρn + ρi, eps())                      # neutral mass fraction
    (x_ion, ξn, ne_eff)
end

# --- assemble the vertical profile ------------------------------------------
# sub-photosphere: the real envelope, from its top (τ≈1) down a few scale heights.
r_env = env.r; ρ_env = env.ρ; T_env = env.Tk; P_env = env.P
r_surf = r_env[end]                                    # τ≈1 boundary radius
# keep the outer ~12 pressure scale heights of the envelope (near-surface only)
Hp_surf = P_env[end] / (ρ_env[end] * GM / r_surf^2)
z_bot   = -12 * Hp_surf
imin = findfirst(r -> r - r_surf ≥ z_bot, r_env)
imin === nothing && (imin = 1)
zsub  = r_env[imin:end] .- r_surf
ρsub  = ρ_env[imin:end]; Tsub = T_env[imin:end]; Psub = P_env[imin:end]

# chromosphere placeholder: isothermal-hydrostatic column at the top T, appended above τ=1.
Ttop = T_env[end]; Ptop = P_env[end]
Hp_top = Ptop / (ρ_env[end] * GM / r_surf^2)
nup = 120
zup = collect(range(zsub[end] + Hp_top/20, zsub[end] + n_Hp_up*Hp_top; length=nup))[2:end]
# integrate dP/dz = -ρ g upward with T fixed = Ttop (ρ from the EOS at (P,Ttop))
function chromosphere(zup, z0, P0, Ttop)
    zout = Float64[]; ρout = Float64[]; Pout = Float64[]
    Pc = P0; zc = z0
    for zt in zup
        dz = zt - zc
        ρc = rho_from_PT(eos, Pc, Ttop)
        g  = GM / (r_surf + zc)^2
        Pc = Pc * exp(-dz * ρc * g / Pc)               # local exponential step (T fixed)
        zc = zt
        push!(zout, zt); push!(ρout, rho_from_PT(eos, Pc, Ttop)); push!(Pout, Pc)
    end
    (zout, ρout, Pout)
end
zup_full, ρup, Pup = chromosphere(zup, zsub[end], Ptop, Ttop)

zall = vcat(zsub, zup_full)
ρall = vcat(ρsub, ρup)
Tall = vcat(Tsub, fill(Ttop, length(zup_full)))
Pall = vcat(Psub, Pup)

# --- write the table --------------------------------------------------------
outf = joinpath(@__DIR__, @sprintf("background_%02ddeg.txt", round(Int, θring_deg)))
open(outf, "w") do io
    @printf(io, "# Stellar2D near-surface background for the SurfaceCode pilot\n")
    @printf(io, "# star: M=%.3e kg (0.5 Msun), R=%.3e m (2 Rsun), B_star=%.1f G, theta_ring=%.1f deg\n",
            Mstar, Rstar, Bstar*1e4, θring_deg)
    @printf(io, "# composition X=%.3f; tau=1 surface at r=%.6e m (r/Rstar=%.4f)\n",
            Xh, r_surf, r_surf/Rstar)
    @printf(io, "# z=0 at tau=1 surface; z<0 sub-photospheric (REAL MLT envelope), z>0 chromosphere\n")
    @printf(io, "#   (z>0 is an ISOTHERMAL-HYDROSTATIC PLACEHOLDER at T=%.0f K — replace for production)\n", Ttop)
    @printf(io, "# ionization: H+He Saha + metal e- floor x_metal=%.1e per H above T=%.0f K (STELLAR2D_XMETAL to change);\n", x_metal, T_metal)
    @printf(io, "#   without the metal floor the H-only x_ion~1e-7 at the surface OVERestimates ambipolar damping.\n")
    @printf(io, "# columns: z[m]  r/Rstar  rho[kg/m3]  T[K]  P[Pa]  x_ion(ne/nH)  xi_neutral  n_e[m-3]  B_r[T]  B_th[T]  |B|[T]  V_A[m/s]  c_s[m/s]  H_p[m]\n")
    for k in eachindex(zall)
        r = r_surf + zall[k]; ρ = ρall[k]; T = Tall[k]; P = Pall[k]
        xion, ξn, ne = ion_state(ρ, T)
        br = Br(r); bt = Bθ(r); bmag = hypot(br, bt)
        VA = bmag / sqrt(μ0 * ρ)
        cs = sqrt(5/3 * P / ρ)                          # approx (γ=5/3); surface code recomputes
        Hp = P / (ρ * GM / r^2)
        @printf(io, "% .6e % .6f % .6e % .6e % .6e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e\n",
                zall[k], r/Rstar, ρ, T, P, xion, ξn, ne, br, bt, bmag, VA, cs, Hp)
    end
end

# --- console summary --------------------------------------------------------
kτ = length(zsub)                                       # index of the τ=1 surface
@printf("wrote %s  (%d points: %d sub-photospheric + %d chromosphere placeholder)\n",
        outf, length(zall), kτ, length(zup_full))
xs, ξs, nes = ion_state(ρall[kτ], Tall[kτ])
@printf("at tau=1 (z=0): rho=%.2e kg/m3, T=%.0f K, x_ion=%.2e, xi_neutral=%.4f, V_A=%.1f km/s, H_p=%.0f km\n",
        ρall[kτ], Tall[kτ], xs, ξs, Br(r_surf)/sqrt(μ0*ρall[kτ])/1e3, Pall[kτ]/(ρall[kτ]*GM/r_surf^2)/1e3)
@printf("vertical span: z=%.2e .. %.2e m  (%.1f H_p below to %.1f H_p above tau=1)\n",
        zall[1], zall[end], -zall[1]/Hp_surf, zall[end]/Hp_top)
