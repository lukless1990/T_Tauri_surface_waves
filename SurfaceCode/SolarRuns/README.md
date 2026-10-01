# SolarRuns — a flare-triggered Moreton wave with the SurfaceCode solver

Adaptation of the SurfaceCode pilot (`../mhd2d.jl`, validated explicit 2.5D MHD) from the
T Tauri surface-wave problem to the **Sun**: a flare pressure pulse in the low corona
launches a fast-mode shock whose lower flank sweeps the chromosphere — the **Moreton
wave** scenario (Uchida 1968, blast-wave variant per Vršnak & Cliver 2008).

All adopted parameter values and their literature references: **`solar_parameters.pdf`**.

Two background modes (`MODE` env):

- **`krause` (default, production)** — the published Moreton-simulation recipe (Krause et
  al. 2015): gravity-free isobaric 1.6 MK corona (p = 2.65e-3 Pa, n ≈ 1.2e8 cm⁻³) with a
  dense chromospheric slab (T = 1e4 K, contrast 160) below 5 Mm; pulse Δp/p = 50 at
  z = 35 Mm for 50 s. Exact static background at any resolution ⇒ robust; our uniform
  250 km cells are finer than the published run's finest AMR cell (390 km).
- **`valc` (high-res variant)** — full hydrostatic VAL-C stratification (photosphere,
  T-minimum, chromosphere, TR, corona). Hydrostatics fixes the chromosphere→corona
  density contrast at ~1e5, so the TR density scale length is ~0.2–0.3 Mm for *any*
  T(z) — needs dz ≲ 75 km to survive the blast impact (multi-hour run; dz ≥ 100 km
  attempts die by TR-cell positivity loss at the flare column — measured, see git log).

## Files

| file | what |
|---|---|
| `solar_moreton.jl` | background builder + run script (uses `../mhd2d.jl`) |
| `plot_solar.py` | figures from a run directory |
| `solar_parameters.tex/.pdf` | adopted values, driver calibration, references |
| `run_<label>/` | run output: background table, 2D snapshots, distance–time traces |
| `plots_<label>/` | figures for a run |

## Usage

```bash
julia -t 14 SurfaceCode/SolarRuns/solar_moreton.jl        # production defaults
python3 SurfaceCode/SolarRuns/plot_solar.py SurfaceCode/SolarRuns/run_<label>
```

Key env knobs: `A_PULSE` (flare overpressure ratio, default 25), `B0_G` (quiet-Sun field,
default 5 G), `THETA_B_DEG`, `LX_MM/LZ_MM/DX_KM/DZ_KM`, `TEND_S`, `RUN_LABEL`.

## Diagnostics

- `trace_chromo_vz.txt` — v_z(x,t) at z=1.5 Mm: the Moreton ("Hα") distance–time stack.
- `trace_corona_dpp.txt` — Δp/p(x,t) at z=8 Mm: the coronal (EUV-wave) front.
- `snap_{vz,vx,dpp}_###.txt` — 2D fields every 30 s (x-subsampled by `xsub` in `meta.txt`).
- `background.txt` — the stratification and the V_fast(z) waveguide profile.

## Design notes

- Uniform, near-vertical B0 → lateral propagation is a perpendicular fast wave,
  V_f = √(c_s²+V_A²); the V_f(z) minimum at the temperature-minimum layer + monotonic rise
  into the corona is the refraction waveguide the Moreton mechanism needs.
- Stability stack (all in the run script; solver core untouched): per-stage bounded-state
  protection (relative ρ floor + velocity ceiling + e_int clamp — HLL/MUSCL has no
  positivity guarantee and the RK2 *intermediate* stage is where the slab/TR interface
  first fails under the blast), a sacrificial guard column at the interface under the
  flare, and a post-pulse flare-site relaxer (in ideal MHD the hot wake has no physical
  relaxation channel and otherwise rings for the whole run). None touch the escaping
  front — kinematics are read at x−x_fl > 30 Mm (see `solar_parameters.pdf`).
- Sponges (velocity+e_int → hydrostatic background; mass-conserving, same design as the
  pilot's `paper_run.jl`) at top and both lateral boundaries.

## Production result (run_production: A=100, B0=3 G, 250 km cells)

Compared against the observed 2006 Dec 6 Moreton wave (Balasubramaniam et al. 2010,
ApJ 723, 587 — PDF in `RefMoreton/`; same event simulated by Krause et al. 2015),
constant-deceleration fit per their Eq. 1 (`compare_observed.py`, fig5):

| | simulation | observed |
|---|---|---|
| mean front speed | **857 km/s** | ~850 km/s |
| fit V0 / deceleration | 891 km/s / −0.42 km/s² | ~1125 km/s / −0.87 km/s² |
| final tracked speed | 813 km/s | ~700–800 km/s |
| chromo. swing (far field) | 2.0 km/s | 1–4 km/s |
| chromo. signal range | 31–126 Mm | ~120–550 Mm |
| coronal front tracked to | 204 Mm (dp/p<3%) | EUV counterpart ~R_sun |

Propagation-phase kinematics reproduced; visibility range shorter (minimum-observable
pulse, 2D slab, 210 Mm domain vs an extreme X6.5 event). See solar_parameters.pdf §6.

## Referee test — driver-limited, not dissipation-limited (run_extreme)

Objection: does the impulsive run's short range (204 Mm vs observed ~550 Mm) mean the
solver over-dissipates (which would undermine the T Tauri decay-length result)? Tested by
adding a sustained CME-piston driver (`DRIVE=piston`: steady deposition over τ_drive=180 s
with a laterally expanding edge, V_pist=800 km/s) in a 600 Mm domain — SAME solver,
dissipation physics, and 250 km resolution, only the driver and box change.

| run | mean speed | front reached | driver |
|---|---|---|---|
| impulsive blast (A=100, 210 Mm box) | 857 km/s | 204 Mm | one-shot |
| sustained CME-piston (A=150, 600 Mm box) | 833 km/s | **553 Mm** | continuous |
| observed 2006 Dec 6 | ~850 km/s | ~550 Mm | CME |

The sustained-driver front reaches 553 Mm ≈ the observed ~550 Mm at the same mean speed:
**the short impulsive range is driver-limited, not dissipation-limited.** Cost: the extreme
piston over-drives the chromospheric swing (~130 km/s vs obs 1–4) — the two runs bracket
reality (impulsive = right amplitude, under-range; piston = right range, over-driven).
`compare_drivers.py` → fig6_driver_test.png; full argument in solar_parameters.pdf §6.

### Mid-range driver: the range-vs-swing trade-off scan

Sought a single sustained driver hitting BOTH the observed ~550 Mm range and the
1–4 km/s chromospheric swing. Scan (`tradeoff.py` → fig7):

| run | A (width) | mean | range | swing |
|---|---|---|---|---|
| mid20 | 20 | 726 | 439 Mm | 27 km/s |
| widegentle | 25 (σ=15) | 802 | 553 Mm | 78 km/s |
| mid50 | 50 | 778 | 548 Mm | 63 km/s |
| extreme | 150 | 831 | 553 Mm | 133 km/s |
| **observed** | — | ~850 | ~550 Mm | 1–4 km/s |

Range saturates ~550 Mm for any A≳25 (sustained-driving + threshold, not amplitude);
swing scales with deposited energy — widening the driver made it *worse*. The sweet-spot
corner (550 Mm AND few-km/s swing) is unreachable with one knob **because the 2D slab omits
geometric spreading** (§5 simpl. iv): in 3D the arc-front's local amplitude falls as ~r^−1/2,
so a real front reaching 550 Mm is intrinsically weak there. This is the *conservative*
direction — with spreading, amplitude would decay faster — which strengthens the
anti-over-dissipation argument rather than weakening it. Full discussion: solar_parameters.pdf §6.1.
