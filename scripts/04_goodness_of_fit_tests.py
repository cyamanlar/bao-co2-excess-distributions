"""Part D - Chi-square and Kolmogorov-Smirnov tests (alpha = 0.05)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import bao  # noqa: E402
import gof  # noqa: E402
import plotstyle as ps  # noqa: E402

ALPHA, M = 0.05, 10
ps.apply()
_, x = bao.excess_sample()
N = len(x)
fits = gof.fit_all(x)

# ---- chi-square, m = 10 equal-probability classes
chi = {k: gof.chi_square(x, d, m=M, alpha=ALPHA) for k, (d, _) in fits.items()}
pd.DataFrame([{"distribution": k, "chi2": v["chi2"], "df": v["df"], "chi2_crit": v["crit"],
               "p_value": v["p_value"], "decision": "reject H0" if v["reject"] else "accept H0"}
              for k, v in chi.items()]).round(4).to_csv(bao.RESULTS / "D_chi_square.csv", index=False)
cls = {"class": np.arange(1, M + 1)}
for k, v in chi.items():
    cls[f"{k}_upper_limit"] = np.r_[v["edges"], np.inf]
    cls[f"{k}_N_i"] = v["Ni"]
pd.DataFrame(cls).round(3).to_csv(bao.RESULTS / "D_chi_square_classes.csv", index=False)

# ---- sensitivity to the number of classes
sens = []
for m in (6, 8, 10, 12, 15, 20):
    row = {"m": m}
    for k, (d, _) in fits.items():
        r = gof.chi_square(x, d, m=m, alpha=ALPHA)
        row[k] = r["chi2"]; row["chi2_crit"] = r["crit"]
    sens.append(row)
pd.DataFrame(sens).round(2).to_csv(bao.RESULTS / "D_chi_square_sensitivity.csv", index=False)

# ---- Kolmogorov-Smirnov with table and Monte Carlo (estimated-parameter) critical values
ks_rows = []
for k, (d, p) in fits.items():
    D, x_at = gof.ks_statistic(x, d)
    c_tab = gof.ks_critical_table(N, ALPHA)
    c_mc = gof.ks_critical_mc(k, N, shape=p.get("alpha"), alpha=ALPHA)
    ks_rows.append({"distribution": k, "Delta": D, "x_at_max": x_at,
                    "crit_table": c_tab, "decision_table": "reject H0" if D >= c_tab else "accept H0",
                    "crit_estimated_params": c_mc,
                    "decision_estimated_params": "reject H0" if D >= c_mc else "accept H0"})
ks = pd.DataFrame(ks_rows)
ks.round(4).to_csv(bao.RESULTS / "D_kolmogorov_smirnov.csv", index=False)

# ---- figure 5: observed counts per class
fig, ax = plt.subplots(figsize=(9, 3.8))
wd = 0.16
for j, (k, v) in enumerate(chi.items()):
    ax.bar(np.arange(1, M + 1) + (j - 2) * wd, v["Ni"], width=wd * 0.9, color=ps.COLOR[k], label=k)
ax.axhline(N / M, color=ps.INK, lw=1, ls="--")
ax.text(M + 0.45, N / M + 2, f"expected N/m = {N / M:.1f}", ha="right", fontsize=8, color=ps.INK2)
ax.set_xticks(range(1, M + 1))
ax.set_xlabel("Equal-probability class i"); ax.set_ylabel("Observed count $N_i$")
ax.set_title("Chi-square test: observed counts in 10 equal-probability classes")
ax.legend(ncol=5, loc="upper left", fontsize=8)
fig.tight_layout(); fig.savefig(bao.FIGURES / "fig5_chi_square_classes.png"); plt.close(fig)

# ---- figure 6: empirical vs fitted CDF
xs = np.sort(x)
fig, ax = plt.subplots(figsize=(9, 4.3))
ax.step(xs, np.arange(1, N + 1) / N, where="post", color=ps.INK, lw=1.2, label="Empirical $F^*(x)$")
xx = np.logspace(-3, np.log10(50), 800)
for k, (d, _) in fits.items():
    ax.plot(xx, d.cdf(xx), color=ps.COLOR[k], ls=ps.LINE[k], lw=1.6, label=k)
ax.set_xscale("log"); ax.set_xlim(0.01, 50); ax.set_ylim(0, 1.01)
ax.set_xlabel("x (ppm, log axis)"); ax.set_ylabel("Cumulative probability")
ax.set_title("Kolmogorov-Smirnov test: empirical and fitted distribution functions")
ax.legend(loc="upper left")
fig.tight_layout(); fig.savefig(bao.FIGURES / "fig6_ks_cdf.png"); plt.close(fig)

print(pd.read_csv(bao.RESULTS / "D_chi_square.csv").to_string(index=False))
print(pd.read_csv(bao.RESULTS / "D_chi_square_sensitivity.csv").to_string(index=False))
print(ks.round(4).to_string(index=False))
