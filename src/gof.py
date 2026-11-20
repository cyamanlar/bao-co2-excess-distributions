"""Fitting and goodness-of-fit tools for two-parameter continuous distributions.

Estimators (Wilks, Statistical Methods in the Atmospheric Sciences):
  Normal     mu = mean, sigma = s
  Lognormal  mu_y, sigma_y = mean and s of ln x
  Gamma      Thom (1958) approximation to the maximum-likelihood estimator
  Weibull    Justus et al. (1978): alpha = (s/mean)^-1.086
  Gumbel     method of moments
"""
import numpy as np
from scipy import special, stats

DISTRIBUTIONS = ["Normal", "Lognormal", "Gamma", "Weibull", "Gumbel"]
EULER = 0.5772156649


def fit(name, x):
    """Return (frozen scipy distribution, parameter dict)."""
    x = np.asarray(x, float)
    m, s = x.mean(), x.std(ddof=1)
    if name == "Normal":
        return stats.norm(m, s), {"mu": m, "sigma": s}
    if name == "Lognormal":
        y = np.log(x)
        my, sy = y.mean(), y.std(ddof=1)
        return stats.lognorm(sy, scale=np.exp(my)), {"mu_y": my, "sigma_y": sy}
    if name == "Gamma":
        D = np.log(m) - np.log(x).mean()
        a = (1 + np.sqrt(1 + 4 * D / 3)) / (4 * D)
        return stats.gamma(a, scale=m / a), {"alpha": a, "beta": m / a, "D": D}
    if name == "Weibull":
        a = (s / m) ** (-1.086)
        b = m / special.gamma(1 + 1 / a)
        return stats.weibull_min(a, scale=b), {"alpha": a, "beta": b}
    if name == "Gumbel":
        b = s * np.sqrt(6) / np.pi
        return stats.gumbel_r(m - EULER * b, b), {"xi": m - EULER * b, "beta": b}
    raise ValueError(name)


def fit_all(x):
    return {k: fit(k, x) for k in DISTRIBUTIONS}


def plotting_position(n):
    """Cunnane plotting position p_i = (i - 0.4) / (n + 0.2)."""
    i = np.arange(1, n + 1)
    return (i - 0.4) / (n + 0.2)


def probability_plot(x, dist):
    xs = np.sort(x)
    p = plotting_position(len(xs))
    q = dist.ppf(p)
    F = dist.cdf(xs)
    return dict(x=xs, p=p, q=q, F=F,
                r_qq=float(np.corrcoef(xs, q)[0, 1]),
                r_pp=float(np.corrcoef(p, F)[0, 1]))


def chi_square(x, dist, n_par=2, m=10, alpha=0.05):
    """Chi-square test with m equal-probability classes (p_i = 1/m).
    chi2 = sum (N_i - N p_i)^2 / (N p_i), d.f. = m - n_par - 1."""
    n = len(x)
    edges = dist.ppf(np.arange(1, m) / m)
    Ni = np.bincount(np.searchsorted(edges, x, side="right"), minlength=m)
    E = n / m
    chi2 = float(((Ni - E) ** 2 / E).sum())
    df = m - n_par - 1
    crit = float(stats.chi2.ppf(1 - alpha, df))
    return dict(edges=edges, Ni=Ni, expected=E, chi2=chi2, df=df, crit=crit,
                p_value=float(stats.chi2.sf(chi2, df)), reject=chi2 >= crit)


def ks_statistic(x, dist):
    """Delta = max |F(x_i) - F*(x_i)|, checked on both sides of every step of F*."""
    xs = np.sort(x)
    n = len(xs)
    F = dist.cdf(xs)
    i = np.arange(1, n + 1)
    d = np.maximum(np.abs(F - i / n), np.abs(F - (i - 1) / n))
    k = int(np.argmax(d))
    return float(d[k]), float(xs[k])


def ks_critical_table(n, alpha=0.05):
    """Large-sample critical value for a fully specified distribution."""
    c = {0.10: 1.22, 0.05: 1.36, 0.01: 1.63}[alpha]
    return c / np.sqrt(n)


def ks_critical_mc(name, n, shape=None, n_sim=2000, alpha=0.05, seed=1):
    """Critical Delta when the parameters are estimated from the same sample
    (Lilliefors idea): simulate samples, re-fit with the same estimator, take
    the (1 - alpha) quantile of Delta."""
    rng = np.random.default_rng(seed)
    gen = {
        "Normal": lambda: rng.normal(size=n),
        "Lognormal": lambda: rng.lognormal(size=n),
        "Gamma": lambda: rng.gamma(shape or 1.0, size=n),
        "Weibull": lambda: rng.weibull(shape or 1.0, size=n),
        "Gumbel": lambda: rng.gumbel(size=n),
    }[name]
    out = np.empty(n_sim)
    for j in range(n_sim):
        z = gen()
        out[j] = ks_statistic(z, fit(name, z)[0])[0]
    return float(np.quantile(out, 1 - alpha))


def mle_aic_bic(name, x):
    """Maximum-likelihood fit (location fixed at 0 for positive distributions)
    with AIC = -2 lnL + 2k and BIC = -2 lnL + k ln n."""
    x = np.asarray(x, float)
    if name == "Normal":
        d = stats.norm(*stats.norm.fit(x))
    elif name == "Lognormal":
        s, _, sc = stats.lognorm.fit(x, floc=0); d = stats.lognorm(s, scale=sc)
    elif name == "Gamma":
        a, _, sc = stats.gamma.fit(x, floc=0); d = stats.gamma(a, scale=sc)
    elif name == "Weibull":
        a, _, sc = stats.weibull_min.fit(x, floc=0); d = stats.weibull_min(a, scale=sc)
    elif name == "Gumbel":
        d = stats.gumbel_r(*stats.gumbel_r.fit(x))
    elif name == "Exponential":
        d = stats.expon(0, x.mean())
    k = 1 if name == "Exponential" else 2
    ll = float(d.logpdf(x).sum())
    n = len(x)
    return dict(dist=d, loglik=ll, k=k, aic=-2 * ll + 2 * k, bic=-2 * ll + k * np.log(n))
