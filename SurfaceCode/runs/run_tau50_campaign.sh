#!/bin/bash
# run_tau50_campaign.sh — re-run the whole paper campaign at the clump-consistent driver
# duration tau = 50 s (r_c ~ 1e-2 R*, v_ff ~ 300 km/s; see paper Sect. 2.2).
#
# The 1 kG single pulse is already done (output/paper_run_tau50, M0=2.90 -> 10.10 km/s).
# This script does the rest, cheapest first so failures surface early:
#   1. calibrate M0 for 0.1 kG at tau=50            (~4 min)
#   2. stage3 hydro twin + weak control pulse       (~10 min)
#   3. 0.1 kG single pulse, production grid         (~35 min)
#   4. continuous driving, production grid          (~32 min)
#   5. stage4 convergence ladder res x1,2,4,8       (~2.4 h)
#
# tau=150 data is left in place (output/paper_run_{1kG,0p1kG,cont}) — it is the comparison
# data for the driver-duration null reported in Sect. 4.3.
set -u
cd "$(dirname "$0")/../.."   # repository root (scripts are run as SurfaceCode/runs/...)
SC=SurfaceCode/runs
T=50
NT=14
TARGET=10.10          # launched-peak target [km/s]: match the 1 kG tau=50 run exactly
stamp() { date +%H:%M:%S; }
die() { echo "!!! FAILED at $1 ($(stamp)) — campaign aborted"; tail -20 "$2"; exit 1; }

echo "=== tau=50 campaign starting $(stamp) ==="

# ---- 1. calibrate the 0.1 kG boundary amplitude at tau=50 -------------------
echo ">>> [1/5] 0.1 kG calibration  $(stamp)"
TAU_S=$T BSTAR_T=0.01 julia -t $NT --project=. $SC/calib_pulse.jl 3.45 3.70 3.80 3.90 \
    > $SC/logs_tau50_calib_0p1kG.log 2>&1 || die "0.1 kG calibration" $SC/logs_tau50_calib_0p1kG.log
cat $SC/logs_tau50_calib_0p1kG.log

# pick the M0 whose launched peak (col 4) is closest to TARGET
M0=$(awk -v t=$TARGET '$3=="km/s" && $5=="km/s" {d=($4-t); if(d<0)d=-d; if(b==""||d<b){b=d;m=$1}} END{print m}' \
     $SC/logs_tau50_calib_0p1kG.log)
[ -z "$M0" ] && die "M0 parse (no calibration rows)" $SC/logs_tau50_calib_0p1kG.log
echo ">>> selected M0=$M0 for 0.1 kG (target ${TARGET} km/s)"

# ---- 2. stage3 hydro: strong shock twin + weak linear control ---------------
echo ">>> [2/5] stage3 hydro (strong + weak control)  $(stamp)"
TAU_S=$T julia -t $NT --project=. $SC/stage3_pulse.jl \
    > $SC/logs_tau50_stage3.log 2>&1 || die "stage3 hydro" $SC/logs_tau50_stage3.log
grep -E "pulse|strong|weak|res×" $SC/logs_tau50_stage3.log

# ---- 3. 0.1 kG production single pulse -------------------------------------
echo ">>> [3/5] 0.1 kG production pulse, M0=$M0  $(stamp)"
TAU_S=$T BSTAR_T=0.01 RUN_LABEL=paper_run_0p1kG_tau50 \
    julia -t $NT --project=. $SC/paper_run.jl $M0 \
    > $SC/paper_run_0p1kG_tau50.log 2>&1 || die "0.1 kG production" $SC/paper_run_0p1kG_tau50.log
tail -2 $SC/paper_run_0p1kG_tau50.log

# ---- 4. continuous driving --------------------------------------------------
echo ">>> [4/5] continuous driving, M0=2.90  $(stamp)"
TAU_S=$T RUN_LABEL=paper_run_cont_tau50 \
    julia -t $NT --project=. $SC/paper_run_cont.jl 2.90 \
    > $SC/paper_run_cont_tau50.log 2>&1 || die "continuous" $SC/paper_run_cont_tau50.log
tail -3 $SC/paper_run_cont_tau50.log

# ---- 5. convergence ladder --------------------------------------------------
# NOTE: at tau=50 the launch wavelength is ~830 km, so res x1 (dx=100 km) resolves it with
# only ~8 cells/lambda — where test2_numdiss measured L_num ~ 1.4 lambda. Expect the coarse
# rung to be grid-dominated; judge convergence on res x2,4,8.
echo ">>> [5/5] stage4 convergence ladder res x1,2,4,8  $(stamp)"
TAU_S=$T julia -t $NT --project=. $SC/stage4_pulse.jl 1 2 4 8 \
    > $SC/logs_tau50_ladder.log 2>&1 || die "stage4 ladder" $SC/logs_tau50_ladder.log
grep "==>" $SC/logs_tau50_ladder.log

echo "=== tau=50 campaign COMPLETE $(stamp) ==="
