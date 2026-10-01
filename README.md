# SurfaceCode

> **Note:** This code is research software under active development. It has been validated for the setups described in the paper, but it may still contain errors or limitations, and interfaces may change. Please check results carefully, and feel free to report issues.

Explicit 2.5D ideal-MHD patch code for the near-surface layers of a classical T Tauri star,
used to test whether accretion-generated fast-mode waves can travel along the stellar surface
from the accretion ring to the polar wind base (the transport assumption of the
accretion-powered stellar wind model of Cranmer 2008).

Companion code to: L. Gehrig, *Can accretion-driven waves survive the journey along the
surface of T Tauri stars?* (submitted to A&A).

## Numerical method

- Finite-volume ideal MHD with gravity, all three vector components retained, ∂_y ≡ 0.
- HLL fluxes; second-order MUSCL reconstruction with the minmod limiter; SSP-RK2 time stepping.
- Hyperbolic–parabolic (GLM) divergence cleaning.
- Well-balanced: deviations from the magnetohydrostatic background are reconstructed, so an
  unperturbed column stays static to machine precision.
- Threaded with `Threads.@threads` and allocation-free in the time loop. The kernel is
  memory-bandwidth bound and saturates at about 4–8 threads.

## Requirements

- Julia ≥ 1.10. The solver and run scripts use only the standard library
  (`Printf`, `DelimitedFiles`, `Base.Threads`).
- Python 3 with `numpy`, `matplotlib`, `scipy` for the plotting and figure scripts.

## Layout

```
SurfaceCode/
├── src/            solver source
├── runs/           run, test and calibration scripts (Julia + shell)
├── data/           stellar background stratification (input to every T Tauri run)
├── analysis/       plotting of run output (Python)
├── Paper/figscripts/   scripts that produce the paper figures
└── SolarRuns/      solar Moreton-wave benchmark (App. B), own README
```

All scripts are run from the repository root. Run output is written to `output/` (git-ignored);
the solar runs write to `SurfaceCode/SolarRuns/run_<label>/`.

### `src/`: solver

| File | Content |
|---|---|
| `mhd2d.jl` | MHD solver core (module `MHD2D`): HLL fluxes, MUSCL/minmod, SSP-RK2, GLM cleaning, well-balanced gravity source and boundary conditions |
| `patch2d.jl` | hydrodynamic core (module `Patch2D`), used for the Sod tube and the hydrodynamic runs |

### `runs/`: runs, tests, calibration

| File | Purpose |
|---|---|
| `run_shocktests.jl` | Sod and Brio–Wu shock tubes (App. A) |
| `test1_stratified_wave.jl` | driven wave in a stratified magneto-atmosphere vs. the exact solution (App. A) |
| `test2_numdiss.jl` | numerical dissipation floor on ideal eigenmodes (App. A) |
| `calib_pulse.jl` | calibrates the boundary amplitude `M0` to a launched pulse peak of 10 km/s |
| `paper_run.jl` | production single-pulse run: lateral energy flux Φ(x) and its decay (Sects. 4.1–4.2) |
| `paper_run_cont.jl` | continuously driven run (Sect. 4.3) |
| `paper_run_leak.jl` | `paper_run.jl` plus a ledger of the energy and mass crossing the top and bottom boundaries, integrated until the burial point converges; with `M0 = 0` it is the unperturbed control run |
| `stage3_pulse.jl` | hydrodynamic twin and weak-pulse control |
| `stage4_pulse.jl` | resolution ladder (Δx = 100 → 13 km) |
| `run_tau50_campaign.sh` | runs the full paper campaign with the exact settings used |
| `stage4_viz.jl`, `run_res8.sh` | visualization run and the finest ladder rung |
| `profile_step.jl`, `bench.jl` | performance profiling |

### `data/`: stellar background

`background_35deg.txt` and `background_55deg.txt` hold the one-dimensional, hydrostatic
near-surface stratification of the fiducial classical T Tauri star
(M = 0.5 M☉, R = 2 R☉, T_eff = 4000 K, X = 0.70). Every T Tauri run reads one of them and
interpolates it onto its vertical grid; this is the "stellar background" of paper Sect. 3.

**Where they come from.**
- Below the photosphere (z < 0): a mixing-length envelope model of this star
  (Böhm-Vitense 1958, α_MLT = 1.8), integrated in hydrostatic equilibrium with an
  ionization-aware equation of state (coupled Saha equilibria of H and He and H₂ dissociation)
  and Rosseland-mean opacities (Bell & Lin 1994).
- Above the photosphere (z > 0): an isothermal (T = 4000 K) hydrostatic layer, attached as a
  placeholder chromosphere.
- The electron density includes a floor from thermally ionized metals of 10⁻⁴ electrons per H
  nucleus above 2500 K. Without it, H + He Saha ionization alone gives x_ion ≈ 10⁻⁷ at the
  surface and overestimates ion–neutral damping. This only enters the ionization columns.
- The magnetic columns are those of an aligned dipole with a polar field of 1 kG, evaluated at
  the accretion-ring colatitude θ_ring = 35° or 55°.

**What they are used for.**
- The runs take z, ρ and p (columns 1, 3, 5) to build the hydrostatic background that the
  well-balanced scheme holds static; the pressure scale height at z = 0 (H_p ≈ 840 km) sets
  the vertical grid and domain ([−4, +3] H_p).
- The run scripts set the magnetic field themselves (uniform field of the dipole at θ_ring = 35°,
  strength from `BSTAR_T`), so the B, V_A and c_s columns are for reference only.
- The ionization columns (x_ion, ξ_n, n_e) are the input to the analytic ion–neutral damping
  estimates (paper Sect. 2.3); the simulations themselves are ideal MHD.
- The two files have identical ρ, T and p and differ only in the field columns. The paper
  uses `background_35deg.txt`.

**Format.** Plain text, `#` comment header, 394 rows from z = −10.0 to +8.4 Mm, 14 columns
(SI units): `z[m]  r/R⋆  ρ[kg/m³]  T[K]  P[Pa]  x_ion(n_e/n_H)  ξ_neutral  n_e[m⁻³]  B_r[T]  B_θ[T]  |B|[T]  V_A[m/s]  c_s[m/s]  H_p[m]`.

### `analysis/`, `Paper/figscripts/`, `SolarRuns/`

`analysis/plot_*.py` make diagnostic plots of runs in `output/` (written to `analysis/plots/`).
`Paper/figscripts/make_paper_figs.py` and `make_solar_fig.py` reproduce the paper figures.
`SolarRuns/` adapts the same solver to the 2006 December 6 solar Moreton wave (App. B).

## Output

### Where

| Script | Output directory |
|---|---|
| `paper_run.jl`, `paper_run_leak.jl`, `paper_run_cont.jl` | `output/<RUN_LABEL>/` (labels: default `paper_run` / `paper_run_cont`, set with `RUN_LABEL`) |
| `run_shocktests.jl` | `output/shocktests/` |
| `test1_stratified_wave.jl`, `test2_numdiss.jl` | `output/wavetests/` |
| `stage4_viz.jl` | `output/plots/surfacecode/` |
| `calib_pulse.jl`, `stage3_pulse.jl`, `stage4_pulse.jl` | no files: results (calibrated amplitudes, fitted decay lengths) are printed to stdout; redirect to a log |
| `SolarRuns/solar_moreton.jl` | `SurfaceCode/SolarRuns/run_<label>/`, described in `SolarRuns/README.md` |
| `analysis/plot_*.py` | figures in `SurfaceCode/analysis/plots/` |

All runs also print flushed progress to stdout (step, time, time step, peak |v_x|).

### Format

All files are plain text, tab-separated, readable with `numpy.loadtxt` or Julia's
`readdlm`. Units are SI unless the column says otherwise. "Per unit length" means per metre
in the invariant y direction, since the model is 2.5D.

**T Tauri pulse runs** (`paper_run.jl`, `paper_run_leak.jl`):

| File | Content |
|---|---|
| `meta.txt` | `key=value` lines: run parameters (`M0`, boundary drive and fast speed [km/s], `Lx_Rstar`, `dx_km`, `Hp_km`, `Nx`, `Nz`, `snap_dt_s`, `Bstar_kG`, `B0_G`, absorbing-layer settings, `tau_s`) |
| `grid_x.txt` | cell-centre x / R⋆, length Nx |
| `grid_z.txt` | cell-centre z / H_p, length Nz |
| `snap_vx_NNN.txt` | lateral velocity v_x [km/s] on the full grid: an Nx × Nz matrix, row i ↔ `grid_x[i]`, column j ↔ `grid_z[j]` |
| `snap_flux_NNN.txt` | two columns: x / R⋆, running energy flux Φ(x) [J/m] accumulated up to that snapshot |
| `snap_times.txt` | one row per snapshot (row NNN ↔ files `_NNN`): time [s], domain peak \|v_x\| [km/s] |
| `flux_final.txt` | as `snap_flux_NNN.txt`, at the end of the run |

Snapshots are written at t = 1.5 τ (the launched pulse), then every 30 min of stellar time,
and once at the end. Φ(x) is the lateral total-energy flux, integrated over the full height
of the domain and over time:
Φ(x) = ∫∫ [(E + p⋆) v_x − B_x (v·B)/μ₀] dz dt.

**Vertical energy ledger** (`paper_run_leak.jl` only):

| File | Content |
|---|---|
| `leak_NNN.txt` | five columns per x cell: x / R⋆, E_top, E_bot [J/m²], M_top, M_bot [kg/m²]: cumulative energy and mass that crossed the top (z = +3 H_p) and bottom (z = −4 H_p) boundaries above/below that cell, positive = leaving the domain. Multiply by Δx and sum over x for the totals per unit length [J/m, kg/m] |
| `leak_final.txt` | as `leak_NNN.txt`, at the end of the run |
| `burial_history.txt` | time [s], flux-burial point x / R⋆ (first x > 0.01 R⋆ with Φ/Φ₀ < 7.8×10⁻⁴) at each 30-min snapshot |
| `meta.txt` | additionally `leak_ledger`, `conv_tol_per_h`, `tmax_h`, `t_cross_s` |

E is the conserved MHD total energy without gravitational potential energy; the potential-energy
flux through a boundary at height z_b is g z_b × the mass flux. `leak_final.txt` and
`burial_history.txt` are only written when a run ends by itself (convergence or `TMAX_H`), not
when it is interrupted; the 30-min `leak_NNN` / `snap_flux_NNN` files are always available.

**Continuously driven run** (`paper_run_cont.jl`): `meta.txt`, `grid_x.txt`, `grid_z.txt`,
`snap_vx_NNN.txt` and `snap_times.txt` as above, plus `flux_steady.txt`: x / R⋆ and the
time-averaged, height-integrated lateral energy flux ⟨F(x)⟩ [W/m] over the trailing steady window.

**Validation tests** (dimensionless code units):

| File | Columns |
|---|---|
| `shocktests/sod.txt` | x, ρ, u, p (N_x = 200, t = 0.2) |
| `shocktests/briowu_400.txt`, `briowu_ref.txt` | x, ρ, v_x, p, B_y (N_x = 400 and the N_x = 1600 reference, t = 0.1) |
| `wavetests/test1_mhd.txt`, `test1_hydro.txt`, `test1_nonlinear.txt` | z, u_z numerical, u_z exact, B_x1 numerical, B_x1 exact (512 cells per wavelength) |
| `wavetests/test1_convergence.txt` | B₀, cells per wavelength, error norms ε(u_z) and ε(B_x), minimum and maximum over height of the numerical/exact amplitude ratio (fitted per wavelength band), maximum phase error [rad] |
| `wavetests/test2_numdiss.txt` | mode (`fast`/`alfven`), cells per wavelength, amplitude decay per period, numerical energy-flux e-folding length L_num / λ |

## Quick start

```bash
# validation
julia -t 4 SurfaceCode/runs/run_shocktests.jl
julia -t 4 SurfaceCode/runs/test1_stratified_wave.jl

# production single pulse, B_star = 1 kG, driver duration 50 s, M0 = 2.90
TAU_S=50 BSTAR_T=0.1 RUN_LABEL=paper_run_tau50 julia -t 8 SurfaceCode/runs/paper_run.jl 2.90

# same with the vertical energy ledger, run until the burial point converges
TAU_S=50 BSTAR_T=0.1 RUN_LABEL=paper_run_tau50_leak julia -t 8 SurfaceCode/runs/paper_run_leak.jl 2.90

# unperturbed control (no pulse), 15 h of stellar time
TAU_S=50 BSTAR_T=0.1 TMAX_H=15 CONV_TOL=0 RUN_LABEL=paper_run_tau50_quiet \
    julia -t 8 SurfaceCode/runs/paper_run_leak.jl 0.0
```

`BSTAR_T` is the dipole polar field in tesla (0.1 = 1 kG, 0.01 = 0.1 kG). For the remaining runs
and their exact arguments, see `SurfaceCode/runs/run_tau50_campaign.sh`. A production run takes
roughly 1 h on an 8-core workstation.

## Parameters

The production scripts (`runs/paper_run.jl`, `paper_run_leak.jl`, `paper_run_cont.jl`) take their
parameters from three places: environment variables and a command-line argument (for what is
varied between runs), `const` definitions at the top of the script, and a few values set inside
`main()`. The table names the constant or variable to look for.

| Parameter | Paper value | Where to change it |
|---|---|---|
| Polar dipole field B⋆ | 1 kG / 0.1 kG | env `BSTAR_T` in tesla (0.1 = 1 kG, 0.01 = 0.1 kG; default 0.1). In `paper_run_cont.jl` it is `const Bstar` |
| Driver amplitude M0 (boundary v_x = M0·c_f) | 2.90 (1 kG), 3.45 (0.1 kG) for τ = 50 s | first command-line argument; default `M0_DEFAULT = 2.63`, which is the calibration for τ = 150 s. Recalibrate with `calib_pulse.jl` whenever τ, B⋆ or the background change |
| Driver duration τ | 50 s | env `TAU_S`. **The default is 150 s**, so set `TAU_S=50` to reproduce the paper |
| Driver vertical extent | Gaussian, width 1.5 H_p | `zw(j)` in `main()` |
| Continuous-driving period | 2τ | follows from `TAU_S` (`paper_run_cont.jl`) |
| Lateral domain size | 0.1 R⋆ | `const LXFRAC` |
| Vertical domain | z ∈ [−4, +3] H_p | `zlo = -4Hp; Lz = 7Hp` in `main()` |
| Resolution | 24 cells per H_p (Δz ≈ 35 km), Δx = 50 km | `const CPH`: it sets both Δz = H_p/CPH and Δx = 1200 km/CPH |
| CFL number | 0.4 | `dt = 0.4*min(g.dx,g.dz)/ch` in the time loop |
| Run length | front transit + τ (paper_run, cont); until converged (leak) | `const TEXTRA` (multiple of the transit time); in `paper_run_leak.jl` env `TMAX_H` (hard cap, stellar hours, default 30) and `CONV_TOL` (burial drift per hour that ends the run, default 0.01; 0 disables) |
| Ambient burial level | Φ/Φ₀ = 7.8×10⁻⁴ = (0.28/10)² | `const PHI_AMB` (`paper_run_leak.jl`, convergence check only) |
| Snapshot cadence | 30 min | `const SNAP_DT` [s] |
| Steady-state averaging window | 30 min | `const WIN` [s] (`paper_run_cont.jl`) |
| Lateral absorbing layers | 40 cells, rate 0.1 s⁻¹ | `const SP_NCELL`, `SP_RATE`; env `LEFT_SPONGE=0` turns off the left layer |
| Ring colatitude θ_ring | 35° | `θ = deg2rad(35.0)` in `main()`; also pick the matching `data/background_*deg.txt` (`read_bg(...)` in `main()`) |
| Star (M⋆, R⋆, gravity) | 0.5 M☉, 2 R☉ | `const Mstar`, `Rstar` at the top of each script. The background file belongs to this star, so a different star needs a new background file in `data/` |
| Adiabatic index γ | 5/3 | `const γc` |
| Output directory name | — | env `RUN_LABEL` |
| Threads | 8 | `julia -t N`; the kernel saturates at about 4–8 threads |

Solver settings (Riemann solver, limiter, GLM cleaning: c_h equal to the maximum fast speed, with ψ
damped by exp(−0.18 c_h Δt / min(Δx, Δz)) each step) live in `src/mhd2d.jl` and are shared by all runs, including the
solar benchmark.

Other scripts:
- `calib_pulse.jl`: the amplitudes to test are command-line arguments (`calib_pulse.jl 2.8 2.9 3.0`); `TAU_S` and `BSTAR_T` as above.
- `stage3_pulse.jl`: first argument is the lateral domain in R⋆ (default 0.036); `TAU_S` as above.
- `stage4_pulse.jl`: the arguments are the resolution multipliers of the ladder (`stage4_pulse.jl 1 2 4 8`); `TAU_S` as above.
- `SolarRuns/solar_moreton.jl`: environment variables, listed in `SolarRuns/README.md`.

## License

MIT, see `LICENSE`.
