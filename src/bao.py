"""Data handling for the BAO flask CO2 record: valid samples, daily means,
background curve and polluted-day (event) definition."""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "bao_co2_pfp_flask.csv"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

THRESHOLD = 5.0      # ppm above the background curve -> polluted day
CLIP_K = 2.0         # residuals > CLIP_K * sigma are excluded when fitting the curve
N_HARMONICS = 2      # annual + semi-annual harmonics
T_REF = 2000.0       # reference year for the linear trend (numerical conditioning only)


def load_samples(path=DATA):
    """All flask samples; a sample is valid if the first QC character is '.'."""
    df = pd.read_csv(path, dtype={"qc_flag": str, "flask_id": str, "event_number": str})
    df["valid"] = df["qc_flag"].str[0] == "."
    df["date"] = df["datetime_utc"].str[:10]
    return df


def daily_means(samples):
    """Average the valid samples of each day (time = mean decimal year of the samples)."""
    v = samples[samples["valid"]]
    d = v.groupby("date").agg(n_samples=("co2_ppm", "size"),
                              co2=("co2_ppm", "mean"),
                              t=("decimal_year", "mean")).reset_index()
    return d


def design(t, n_harm=N_HARMONICS):
    cols = [np.ones_like(t), t - T_REF]
    for k in range(1, n_harm + 1):
        cols += [np.cos(2 * np.pi * k * t), np.sin(2 * np.pi * k * t)]
    return np.column_stack(cols)


def background_curve(t, c, k=CLIP_K, max_iter=50):
    """Linear trend + 2 harmonics fitted by least squares. Days whose residual is
    more than k standard deviations ABOVE the curve (polluted air) are excluded
    and the fit is repeated until the set of excluded days no longer changes."""
    A = design(t)
    keep = np.ones(len(c), bool)
    for _ in range(max_iter):
        coef = np.linalg.lstsq(A[keep], c[keep], rcond=None)[0]
        r = c - A @ coef
        new = r <= k * r.std(ddof=1)
        if (new == keep).all():
            break
        keep = new
    coef = np.linalg.lstsq(A[keep], c[keep], rcond=None)[0]
    return A @ coef, coef, keep


def build_daily(path=DATA):
    """Daily table with background b(t), anomaly Delta = C - b and event flag."""
    d = daily_means(load_samples(path))
    b, coef, keep = background_curve(d["t"].to_numpy(), d["co2"].to_numpy())
    d["background"] = b
    d["delta"] = d["co2"] - d["background"]
    d["used_in_fit"] = keep
    d["event"] = d["delta"] > THRESHOLD
    d["excess"] = np.where(d["event"], d["delta"] - THRESHOLD, np.nan)
    return d, coef


def excess_sample(path=DATA):
    """x = Delta - 5 ppm on polluted days (the variable analysed in this project)."""
    d, _ = build_daily(path)
    e = d[d["event"]]
    return e["date"].to_numpy(), e["excess"].to_numpy()
