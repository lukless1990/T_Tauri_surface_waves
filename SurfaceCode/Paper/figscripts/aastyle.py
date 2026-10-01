# aastyle.py — shared matplotlib style for the A&A paper figures.
import matplotlib as mpl

COL    = 3.464   # 88 mm A&A single-column width [in]
TWOCOL = 7.087   # 180 mm A&A two-column width  [in]

# colorblind-safe palette (Okabe-Ito derived)
C = dict(blue="#0173b2", orange="#de8f05", green="#029e73",
         red="#d55e00", purple="#cc78bc", gray="#555555")

def setup():
    mpl.use("Agg")
    mpl.rcParams.update({
        "font.size": 8, "axes.labelsize": 8.5, "axes.titlesize": 8.5,
        "legend.fontsize": 6.8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
        "font.family": "serif", "mathtext.fontset": "stix",
        "font.serif": ["STIXGeneral", "DejaVu Serif"],
        "axes.linewidth": 0.6, "lines.linewidth": 1.1,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.minor.size": 1.8, "ytick.minor.size": 1.8,
        "legend.framealpha": 0.9, "legend.edgecolor": "0.85",
        "savefig.dpi": 300, "pdf.fonttype": 42,
        "figure.constrained_layout.use": True,
    })
