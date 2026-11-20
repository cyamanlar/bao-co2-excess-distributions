"""Part E - Model selection (AIC/BIC), parameter uncertainty, tail probabilities
and seasonal dependence of the excess."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

import bao  # noqa: E402
import gof  # noqa: E402
import plotstyle as ps  # noqa: E402

ps.apply()
dates, x = bao.excess_sample()
daily, _ = bao.build_daily()
N, n_days = len(x), len(daily)
p_event = N / n_days
fits = gof.fit_all(x)

# ---- AIC / BIC with maximum-likelihood fits (+ one-parameter exponential)
ab = []
for k in gof.DISTRIBUTIONS + ["Exponential"]:
    r = gof.mle_aic_bic(k, x)
    ab.append({"distribution": k, "k": r["k"], "loglik": r["loglik"], "AIC": r["aic"], "BIC": r["bic"]})
ab = pd.DataFrame(ab).sort_values("AIC")
ab["dAIC"] = ab["AIC"] - ab["AIC"].min()
ab["dBIC"] = ab["BIC"] - ab["BIC"].min()
ab.round(3).to_csv(bao.RESULTS / "E_aic_bic.csv", index=False)

# ---- bootstrap confidence interval of the gamma shape (is alpha = 1, exponential, plausible?)
rng = np.random.default_rng(7)
boot = np.array([gof.fit("Gamma", rng.choice(x, N))[1]["alpha"] for _ in range(2000)])
ci = np.quantile(boot, [0.025, 0.975])
pd.DataFrame([{"alpha_hat": fits["Gamma"][1]["alpha"], "ci_low_95": ci[0], "ci_high_95": ci[1],
               "contains_1": bool(ci[0] <= 1 <= ci[1])}]).round(4).to_csv(
    bao.RESULTS / "E_gamma_shape_bootstrap.csv", index=False)

# ---- tail probabilities: fitted vs observed
ex = []
for thr in (10, 20, 30):
    row = {"x_threshold_ppm": thr, "Delta_threshold_ppm": thr + bao.THRESHOLD,
           "observed_P(X>x)": float((x > thr).mean())}
    for k, (d, _) in fits.items():
        row[f"{k}_P(X>x)"] = float(d.sf(thr))
    row["observed_daily_P(Delta>x+5)"] = float((daily["delta"] > thr + bao.THRESHOLD).mean())
    row["Gamma_daily_P(Delta>x+5)"] = p_event * float(fits["Gamma"][0].sf(thr))
    ex.append(row)
pd.DataFrame(ex).round(5).to_csv(bao.RESULTS / "E_exceedance.csv", index=False)

# ---- seasonal dependence
month = np.array([int(s[5:7]) for s in dates])
seasons = {"DJF": (12, 1, 2), "MAM": (3, 4, 5), "JJA": (6, 7, 8), "SON": (9, 10, 11)}
srows = []
for s, ms in seasons.items():
    xs = x[np.isin(month, ms)]
    d, p = gof.fit("Gamma", xs)
    srows.append({"season": s, "N": len(xs), "mean": xs.mean(), "median": np.median(xs),
                  "gamma_alpha": p["alpha"], "gamma_beta": p["beta"]})
seas = pd.DataFrame(srows)
w, su = x[np.isin(month, seasons["DJF"])], x[np.isin(month, seasons["JJA"])]
ks2 = stats.ks_2samp(w, su)
seas.round(4).to_csv(bao.RESULTS / "E_seasonal.csv", index=False)
pd.DataFrame([{"test": "two-sample K-S, DJF vs JJA", "statistic": ks2.statistic, "p_value": ks2.pvalue}]
             ).round(5).to_csv(bao.RESULTS / "E_seasonal_test.csv", index=False)

# ---- figure 7: survival function on a log axis (tail view) + seasonal survival
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4))
xs = np.sort(x)
emp = 1 - (np.arange(1, N + 1) - 0.4) / (N + 0.2)
a1.scatter(xs, emp, s=8, color=ps.INK2, label="Observed", zorder=3)
xx = np.linspace(0.01, 60, 600)
for k in ("Lognormal", "Gamma", "Weibull"):
    a1.plot(xx, fits[k][0].sf(xx), color=ps.COLOR[k], ls=ps.LINE[k], label=k)
a1.set_yscale("log"); a1.set_ylim(1e-3, 1.05); a1.set_xlim(0, 60)
a1.set_xlabel("x (ppm)"); a1.set_ylabel("Pr{X > x}")
a1.set_title("Upper tail: exceedance probability"); a1.legend()
for s, col in (("DJF", ps.COLOR["Normal"]), ("JJA", ps.COLOR["Lognormal"])):
    v = np.sort(x[np.isin(month, seasons[s])]); n = len(v)
    a2.scatter(v, 1 - (np.arange(1, n + 1) - 0.4) / (n + 0.2), s=8, color=col, label=f"{s} observed (N={n})")
    d, _ = gof.fit("Gamma", v)
    a2.plot(xx, d.sf(xx), color=col, lw=1.4, label=f"{s} gamma fit")
a2.set_yscale("log"); a2.set_ylim(1e-3, 1.05); a2.set_xlim(0, 60)
a2.set_xlabel("x (ppm)"); a2.set_title(f"Winter vs summer (K-S p = {ks2.pvalue:.1e})"); a2.legend()
fig.tight_layout(); fig.savefig(bao.FIGURES / "fig7_tail_and_season.png"); plt.close(fig)

print(ab.round(2).to_string(index=False))
print(pd.read_csv(bao.RESULTS / "E_gamma_shape_bootstrap.csv").to_string(index=False))
print(pd.read_csv(bao.RESULTS / "E_exceedance.csv").to_string(index=False))
print(seas.round(3).to_string(index=False), ks2)
