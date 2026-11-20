"""Part A - Daily CO2, background curve, polluted days and the excess variable."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

import bao  # noqa: E402
import plotstyle as ps  # noqa: E402

ps.apply()
bao.RESULTS.mkdir(exist_ok=True); bao.FIGURES.mkdir(exist_ok=True)

samples = bao.load_samples()
daily, coef = bao.build_daily()
x = daily.loc[daily["event"], "excess"].to_numpy()
y = np.log(x)

daily.round(4).to_csv(bao.RESULTS / "A_daily.csv", index=False)

summary = {
    "samples_total": len(samples),
    "samples_valid": int(samples["valid"].sum()),
    "sampling_days": len(daily),
    "first_day": daily["date"].iloc[0], "last_day": daily["date"].iloc[-1],
    "days_used_in_curve_fit": int(daily["used_in_fit"].sum()),
    "curve_a0_ppm": coef[0], "curve_trend_ppm_per_yr": coef[1],
    "curve_seasonal_amplitude_ppm": float(np.hypot(coef[2], coef[3])),
    "event_days_N": len(x), "event_probability": len(x) / len(daily),
    "x_mean": x.mean(), "x_median": np.median(x), "x_std": x.std(ddof=1),
    "x_cv": x.std(ddof=1) / x.mean(), "x_skewness": stats.skew(x, bias=False),
    "x_min": x.min(), "x_max": x.max(),
    "x_q25": np.quantile(x, .25), "x_q75": np.quantile(x, .75),
    "lnx_mean": y.mean(), "lnx_std": y.std(ddof=1), "lnx_skewness": stats.skew(y, bias=False),
}
pd.Series(summary).to_csv(bao.RESULTS / "A_summary.csv", header=["value"])

# ---- figure 1: series with background curve + histogram of the anomaly
t = daily["t"].to_numpy()
fig, (a1, a2) = plt.subplots(2, 1, figsize=(9, 6.2), gridspec_kw={"height_ratios": [1.3, 1]})
ev = daily["event"].to_numpy()
a1.scatter(t[~ev], daily["co2"][~ev], s=4, color=ps.NEUTRAL, label="Daily mean CO$_2$")
a1.scatter(t[ev], daily["co2"][ev], s=6, color=ps.COLOR["Lognormal"], label=r"Polluted day ($\Delta$ > 5 ppm)")
a1.plot(t, daily["background"], color=ps.INK, lw=1.4, label="Background curve b(t)")
a1.set_ylabel("CO$_2$ (ppm)"); a1.set_title("BAO tower (300 m), daily mean CO$_2$")
a1.legend(loc="upper left", markerscale=2)
bins = np.arange(-10, 52, 1)
d = daily["delta"].to_numpy()
a2.hist(d[~ev], bins=bins, color=ps.NEUTRAL, edgecolor=ps.SURFACE, lw=0.4, label="Other days")
a2.hist(d[ev], bins=bins, color=ps.COLOR["Lognormal"], edgecolor=ps.SURFACE, lw=0.4, label="Polluted days")
a2.axvline(bao.THRESHOLD, color=ps.INK, lw=1, ls="--")
a2.set_xlabel(r"$\Delta$ = C − b(t)  (ppm)"); a2.set_ylabel("Number of days")
a2.set_title("Anomaly from the background curve"); a2.legend()
fig.tight_layout(); fig.savefig(bao.FIGURES / "fig1_series_and_anomaly.png"); plt.close(fig)

print(pd.Series(summary).to_string())
