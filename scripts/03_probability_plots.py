"""Part C - Q-Q and P-P plots (Cunnane plotting position)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

import bao  # noqa: E402
import gof  # noqa: E402
import plotstyle as ps  # noqa: E402

ps.apply()
_, x = bao.excess_sample()
fits = gof.fit_all(x)
pp = {k: gof.probability_plot(x, d) for k, (d, _) in fits.items()}

pd.DataFrame([{"distribution": k, "r_QQ": v["r_qq"], "r_PP": v["r_pp"],
               "max_theoretical_quantile": v["q"].max()} for k, v in pp.items()]
             ).round(4).to_csv(bao.RESULTS / "C_probability_plot_correlation.csv", index=False)

for kind in ("qq", "pp"):
    fig, axs = plt.subplots(2, 3, figsize=(9, 6.2))
    for ax, (k, v) in zip(axs.flat, pp.items()):
        if kind == "qq":
            xv, yv = v["q"], v["x"]
            lo, hi = min(0, xv.min()), max(xv.max(), yv.max()) * 1.03
            r = v["r_qq"]
        else:
            xv, yv, lo, hi, r = v["F"], v["p"], 0, 1, v["r_pp"]
        ax.plot([lo, hi], [lo, hi], color=ps.INK, lw=0.8)
        ax.scatter(xv, yv, s=8, color=ps.COLOR[k], edgecolor="none", alpha=0.85, zorder=3)
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
        ax.set_title(f"{k}  (r = {r:.3f})", fontsize=10)
        ax.tick_params(labelsize=8)
    axs.flat[-1].axis("off")
    if kind == "qq":
        fig.supxlabel(r"Theoretical quantile $F^{-1}(p_i)$ (ppm)")
        fig.supylabel(r"Observed quantile $x_{(i)}$ (ppm)")
        name = "fig3_qq_plots.png"
    else:
        fig.supxlabel(r"Theoretical probability $F(x_{(i)})$")
        fig.supylabel(r"Empirical probability $p_i = (i-0.4)/(N+0.2)$")
        name = "fig4_pp_plots.png"
    fig.tight_layout(); fig.savefig(bao.FIGURES / name); plt.close(fig)

print(pd.read_csv(bao.RESULTS / "C_probability_plot_correlation.csv").to_string(index=False))
