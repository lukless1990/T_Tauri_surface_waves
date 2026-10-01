# SurfaceCode

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

## License

MIT, see `LICENSE`.
