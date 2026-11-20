"""Part B - Fit five candidate distributions and superpose them on the histogram."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import bao  # noqa: E402
import gof  # noqa: E402
import plotstyle as ps  # noqa: E402

ps.apply()
_, x = bao.excess_sample()
fits = gof.fit_all(x)

rows = []
for k, (d, p) in fits.items():
    names = list(p)
    rows.append({"distribution": k,
                 "param1": names[0], "value1": p[names[0]],
                 "param2": names[1], "value2": p[names[1]],
                 "mean": d.mean(), "variance": d.var(),
                 "skewness": float(d.stats(moments="s")),
                 "P(X<0)": d.cdf(0)})
rows.append({"distribution": "Sample", "mean": x.mean(), "variance": x.var(ddof=1),
             "skewness": pd.Series(x).skew(), "P(X<0)": 0.0})
tab = pd.DataFrame(rows)
tab.round(4).to_csv(bao.RESULTS / "B_parameters.csv", index=False)

w = 2.0
fig, ax = plt.subplots(figsize=(9, 4.6))
ax.hist(x, bins=np.arange(0, 48, w), color=ps.NEUTRAL, edgecolor=ps.SURFACE, lw=0.8,
        label=f"Observed (N = {len(x)})")
xx = np.linspace(0.05, 46, 600)
for k, (d, _) in fits.items():
    ax.plot(xx, len(x) * w * d.pdf(xx), color=ps.COLOR[k], ls=ps.LINE[k], label=k)
ax.set_xlim(0, 46); ax.set_ylim(0, 110)
ax.set_xlabel("x = CO$_2$ excess above the 5 ppm threshold (ppm)")
ax.set_ylabel("Days per 2-ppm class")
ax.set_title("Histogram of the excess with fitted densities (N·w·f(x))")
ax.legend(ncol=2)
fig.tight_layout(); fig.savefig(bao.FIGURES / "fig2_histogram_fitted_pdfs.png"); plt.close(fig)
print(tab.round(4).to_string(index=False))
