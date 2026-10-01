#!/bin/bash
# Wait for the extended viz run to finish, then start the res×8 MHD convergence run.
cd "$(dirname "$0")/../.."   # repository root (scripts are run as SurfaceCode/runs/...)
while pgrep -f stage4_viz >/dev/null 2>&1; do sleep 30; done
sleep 5
echo ">>> viz done; starting res×8 $(date +%H:%M:%S)"
julia -t 14 SurfaceCode/runs/stage4_pulse.jl 8 > SurfaceCode/runs/stage4_res8.log 2>&1
echo ">>> res×8 done $(date +%H:%M:%S)"
