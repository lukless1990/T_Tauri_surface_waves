#!/usr/bin/env python3
"""front_range.py — front position from the LATEST snapshot of a run (live-monitor helper).

Usage: python3 front_range.py run_<label> [snap_index]

Reports, as distance from the flare site:
  coronal  3%  : rightmost |dp/p| > 0.03 at z_eit  (our standard coronal-front threshold)
  coronal 13%  : rightmost |dp/p| > 0.13           (proxy for Krause's 8% COMPRESSION ratio,
                                                    drho/rho ~= (dp/p)/gamma for gamma=5/3)
  chromo       : rightmost |v_z| > 0.3 km/s at z_ha
Cells within x_fl+30 Mm are excluded (flare-site guard zone), as in compare_observed.py.
"""
import sys, glob, os
import numpy as np

run = (sys.argv[1] if len(sys.argv) > 1 else "run_krfull").rstrip("/")
snaps = sorted(glob.glob(os.path.join(run, "snap_dpp_*.txt")))
if not snaps:
    print(f"{run}: no snapshots yet"); sys.exit(0)
si = int(sys.argv[2]) if len(sys.argv) > 2 else len(snaps)
x  = np.loadtxt(os.path.join(run, "grid_x_Mm.txt"))
z  = np.loadtxt(os.path.join(run, "grid_z_Mm.txt"))
meta = dict(l.strip().split("=") for l in open(os.path.join(run, "meta.txt")))
xfl  = float(meta["x_fl_Mm"]); zha = float(meta["z_ha_Mm"]); zeit = float(meta["z_eit_Mm"])
ts = np.atleast_1d(np.loadtxt(os.path.join(run, "snap_times.txt")))
d  = np.loadtxt(os.path.join(run, f"snap_dpp_{si:03d}.txt"))
v  = np.loadtxt(os.path.join(run, f"snap_vz_{si:03d}.txt"))
je, jh = np.argmin(abs(z - zeit)), np.argmin(abs(z - zha))
m = x > xfl + 30

def far(sig, thr):
    i = np.where((np.abs(sig) > thr) & m)[0]
    return x[i[-1]] - xfl if len(i) else float("nan")

dc, vc = d[:, je], v[:, jh]
f3 = far(dc, 0.03)
# Front amplitude = max |dp/p| in the leading 20 Mm BEHIND the front. NOT the window maximum:
# that is pinned at the guard-zone edge (x_fl+30 Mm) by the flare-site wake and says nothing
# about the propagating front.
lead = (x - xfl > f3 - 20) & (x - xfl <= f3) if np.isfinite(f3) else np.zeros_like(x, bool)
amp = np.abs(dc[lead]).max() if lead.any() else float("nan")
print(f"{os.path.basename(run)} snap {si:03d} t={ts[si-1]:5.0f}s | "
      f"coronal 3%={f3:6.1f} Mm  13%={far(dc,0.13):6.1f} Mm | "
      f"chromo={far(vc,0.30):6.1f} Mm | dp/p at front={amp:6.3f}")
