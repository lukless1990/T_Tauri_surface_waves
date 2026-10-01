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

All scripts are run from the repository root. Run output is written to `output/` (git-ignored).

| File | Purpose |
|---|---|
| `SurfaceCode/mhd2d.jl` | MHD solver core (module `MHD2D`) |
| `SurfaceCode/patch2d.jl` | hydrodynamic core, used by the stage-3 hydro runs |
| `SurfaceCode/background_35deg.txt`, `background_55deg.txt` | near-surface stellar backgrounds (ring colatitude 35° / 55°) |
| `SurfaceCode/export_background.jl` | generates the background files; needs the separate Stellar2D package, so its output is provided here |
| `SurfaceCode/run_shocktests.jl`, `plot_shocktests.py` | Sod and Brio–Wu shock tubes (App. A) |
| `SurfaceCode/test1_stratified_wave.jl`, `plot_wavetests.py` | driven wave in a stratified magneto-atmosphere vs. the exact solution (App. A) |
| `SurfaceCode/test2_numdiss.jl` | numerical dissipation floor on ideal eigenmodes (App. A) |
| `SurfaceCode/calib_pulse.jl` | calibrates the boundary amplitude `M0` to a launched pulse peak of 10 km/s |
| `SurfaceCode/paper_run.jl` | production single-pulse run: decay of the lateral energy flux Φ(x) (Sects. 4.1–4.2) |
| `SurfaceCode/paper_run_cont.jl` | continuously driven run (Sect. 4.3) |
| `SurfaceCode/paper_run_leak.jl` | `paper_run.jl` plus a ledger of the energy and mass crossing the top and bottom boundaries, integrated until the burial point converges (also used with `M0 = 0` as the unperturbed control) |
| `SurfaceCode/stage3_pulse.jl` | hydrodynamic twin and weak-pulse control |
| `SurfaceCode/stage4_pulse.jl` | resolution ladder (Δx = 100 → 13 km) |
| `SurfaceCode/run_tau50_campaign.sh` | runs the full paper campaign with the exact settings used |
| `SurfaceCode/plot_*.py`, `stage4_viz.jl`, `profile_step.jl`, `bench.jl` | diagnostics, visualization and profiling |
| `SurfaceCode/Paper/figscripts/` | scripts that produce the paper figures from `output/` |
| `SurfaceCode/SolarRuns/` | solar Moreton-wave benchmark of the same solver (App. B); see its own README |

## Quick start

```bash
# validation
julia -t 4 SurfaceCode/run_shocktests.jl
julia -t 4 SurfaceCode/test1_stratified_wave.jl

# production single pulse, B_star = 1 kG, driver duration 50 s, M0 = 2.90
TAU_S=50 BSTAR_T=0.1 RUN_LABEL=paper_run_tau50 julia -t 8 SurfaceCode/paper_run.jl 2.90

# same with the vertical energy ledger, run until the burial point converges
TAU_S=50 BSTAR_T=0.1 RUN_LABEL=paper_run_tau50_leak julia -t 8 SurfaceCode/paper_run_leak.jl 2.90

# unperturbed control (no pulse), 15 h of stellar time
TAU_S=50 BSTAR_T=0.1 TMAX_H=15 CONV_TOL=0 RUN_LABEL=paper_run_tau50_quiet \
    julia -t 8 SurfaceCode/paper_run_leak.jl 0.0
```

`BSTAR_T` is the dipole polar field in tesla (0.1 = 1 kG, 0.01 = 0.1 kG). For the remaining runs
and their exact arguments, see `SurfaceCode/run_tau50_campaign.sh`. A production run takes
roughly 1 h on an 8-core workstation.

## License

MIT, see `LICENSE`.
