"""Shared definitions for the figure scripts: data loading, labels, colours, and the Damkohler number."""
import os
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get("DATA", os.path.join(HERE, "..", "data")); OUT = os.environ.get("OUT", os.path.join(HERE, "output")); os.makedirs(OUT, exist_ok=True)   # env overrides, e.g. DATA=../data_v2 OUT=output_v2
plt.rcParams.update({"font.size": 11, "pdf.fonttype": 42})
BLUE, GRAY, RED = "#2166ac", "#525252", "#b2182b"
BEST, HOMOGENEOUS, WORST = (1, -2), (0, 0), (-1, 2)          # (beta_pyr, beta_cal) at |beta_cal| = 2

NETWORKS = pd.read_csv(os.path.join(HERE, "..", "1_inputs", "networks.csv")).set_index("dfn")
SIGMA = NETWORKS["sigma_log_q"]


def load(name):
    return pd.read_csv(os.path.join(DATA, name), keep_default_na=False, na_values=[""])


def save(fig, name):
    for ext, dpi in (("png", 300), ("pdf", 600)):
        fig.savefig(os.path.join(OUT, "%s.%s" % (name, ext)), dpi=dpi, bbox_inches="tight")
    print("saved output/%s.png and .pdf" % name)


def escape_at(ts, t_yr, series="factorial"):
    """U escape [%] at time t_yr, indexed by (dfn, beta_pyr, beta_cal)."""
    d = ts[(ts.series == series) & (ts.t_yr == t_yr)]
    return d.set_index(["dfn", "beta_pyr", "beta_cal"])["U_esc_percent"]


# ---- network-scale Damkohler number (Equations 6-9) and the escape relationship (Equation 10)
N_URANYL, K_U, F_ACETATE = 4.0, 1e-7, 0.17      # uranyl per acetate turnover; U(VI) half-saturation [M]; acetate Monod factor
SECONDS_PER_YEAR = 365.25 * 86400


def damkohler(terms, column="I_over_alpha"):
    """Da = lambda_eff * tau, lambda_eff = N k_b f_Ac <I_O2/alpha>_c / K_U, tau = network transit time / flow-rate factor."""
    tau = terms["dfn"].map(NETWORKS["tau_yr"]).values * SECONDS_PER_YEAR / terms["q_factor"].values
    return N_URANYL * terms["kb"].values * F_ACETATE * np.clip(terms[column].values, 1e-12, None) / K_U * tau


def escape_law(log10_da, log10_m, n):
    return 100.0 * np.exp(-(10 ** (log10_da - log10_m)) ** n)


def fit_escape_law(da, escape):
    from scipy.optimize import curve_fit
    p, _ = curve_fit(escape_law, np.log10(da), escape, p0=[0.0, 0.7])
    return p                                                  # (log10 m, n)


def r_squared(obs, pred):
    return 1 - ((obs - pred) ** 2).sum() / ((obs - obs.mean()) ** 2).sum()


def is_fitting_set(terms):
    """The 180 strong-correlation simulations (|beta_cal| = 2 factorial incl. beta_cal = 0) used to fit Equation 10."""
    return (terms.series == "factorial") & (terms.beta_cal.abs() != 1)
