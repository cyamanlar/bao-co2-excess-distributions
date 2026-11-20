"""Shared, quiet matplotlib style for all figures."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SURFACE = "#fcfcfb"
NEUTRAL = "#c9c8bf"   # observed data (histograms)

# one fixed colour + line style per distribution (identity is never colour-only)
COLOR = {"Normal": "#2a78d6", "Lognormal": "#eb6834", "Gamma": "#1baf7a",
         "Weibull": "#eda100", "Gumbel": "#e87ba4", "Exponential": "#52514e"}
LINE = {"Normal": "-", "Lognormal": "--", "Gamma": "-", "Weibull": "-.",
        "Gumbel": ":", "Exponential": (0, (1, 1))}


def apply():
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE, "savefig.dpi": 150,
        "savefig.bbox": "tight",
        "font.family": "sans-serif", "font.size": 10,
        "text.color": INK, "axes.labelcolor": INK2,
        "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
        "axes.axisbelow": True,
        "xtick.color": MUTED, "ytick.color": MUTED,
        "axes.titlesize": 11, "axes.titleweight": "bold",
        "axes.titlelocation": "left", "legend.frameon": False,
        "lines.linewidth": 2,
    })
