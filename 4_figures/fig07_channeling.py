#!/usr/bin/env python3
"""Figure 7. (a) Pyrite and calcite arrangement effects at 1,500 yr versus flow channeling, with the effects predicted
from Equation 10; (b) U escape relative to the network maximum versus the initial conduit calcite volume fraction."""
from paper import *
from scipy.optimize import curve_fit

ts = load("escape_timeseries.csv"); E = escape_at(ts, 1500); dfns = NETWORKS.index.values; sig = SIGMA.values
dU_pyr = np.array([E[(d, -1, 2)] - E[(d, 1, 2)] for d in dfns])            # anti-aligned - aligned pyrite, calcite aligned
dU_cal = np.array([E[(d, 1, 2)] - E[(d, 1, -2)] for d in dfns])            # aligned - anti-aligned calcite, pyrite aligned

# effects predicted from the Da of the two arrangements being compared (Equation 10, constants from the ensemble fit)
terms = load("conduit_terms.csv"); fit = terms[is_fitting_set(terms) & (terms.t_yr >= 100)]
p = fit_escape_law(damkohler(fit), fit.U_esc_percent.values)
t15 = terms[(terms.series == "factorial") & (terms.t_yr == 1500)].copy(); t15["pred"] = escape_law(np.log10(damkohler(t15)), *p)
P = t15.set_index(["dfn", "beta_pyr", "beta_cal"])["pred"]
pred_pyr = np.array([P[(d, -1, 2)] - P[(d, 1, 2)] for d in dfns]); pred_cal = np.array([P[(d, 1, 2)] - P[(d, 1, -2)] for d in dfns])
rmse = lambda a, b: np.sqrt(((a - b) ** 2).mean())
print("predicted vs simulated effects: calcite RMSE %.1f pp, pyrite RMSE %.1f pp" % (rmse(pred_cal, dU_cal), rmse(pred_pyr, dU_pyr)))

logistic = lambda x, A, x50, w: A / (1 + np.exp(-(x - x50) / w))
L, _ = curve_fit(logistic, sig, dU_cal, p0=[60, 0.5, 0.05], maxfev=5000); s50, w = L[1], L[2]
rng = np.random.default_rng(0); boot = []
for _ in range(5000):
    i = rng.integers(0, 20, 20)
    try: boot.append(curve_fit(logistic, sig[i], dU_cal[i], p0=[60, 0.5, 0.05], maxfev=4000)[0][1])
    except Exception: pass
print("logistic midpoint %.2f (95%% CI %.2f-%.2f), 10-90%% width %.2f-%.2f; midpoint of predicted effects %.2f"
      % (s50, *np.percentile(boot, [2.5, 97.5]), s50 - w * np.log(9), s50 + w * np.log(9), curve_fit(logistic, sig, pred_cal, p0=[60, 0.5, 0.05], maxfev=5000)[0][1]))

fig, (a, b) = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw=dict(width_ratios=[1.25, 1], wspace=0.22))
x = sig.copy(); order = np.argsort(sig); groups = [[order[0]]]          # spread symbols of networks with nearly equal sigma
for k in order[1:]: (groups[-1].append(k) if sig[k] - sig[groups[-1][-1]] < 0.025 else groups.append([k]))
for g in groups:
    if len(g) > 1: x[g] = sig[g].mean() + 0.025 * (np.arange(len(g)) - (len(g) - 1) / 2)
a.axvline(s50, color="0.35", ls="--", lw=1); a.axhline(0, color="k", lw=0.6)
xx = np.linspace(0.05, 1.03, 300); a.plot(xx, logistic(xx, *L), color=BLUE, lw=1, alpha=0.55)
a.scatter(x, dU_cal, c=BLUE, marker="D", s=44, edgecolor="k", lw=0.4, zorder=3, label="$\\Delta U_{cal}^{pyrA}$")
a.scatter(x, dU_pyr, c=RED, marker="o", s=46, edgecolor="k", lw=0.4, zorder=4, label="$\\Delta U_{pyr}^{calA}$")
a.scatter(x, pred_cal, facecolor="none", edgecolor=BLUE, marker="D", s=44, zorder=3, label="$\\Delta U_{cal}^{pyrA}$ (Eq. 10)")
a.scatter(x, pred_pyr, facecolor="none", edgecolor=RED, marker="o", s=46, zorder=4, label="$\\Delta U_{pyr}^{calA}$ (Eq. 10)")
a.set_xlim(0.05, 1.03); a.set_ylim(-5, 80); a.grid(alpha=0.25); a.legend(fontsize=9, loc="upper left")
a.set_xlabel("$\\sigma_{\\log|\\mathbf{q}|}$"); a.set_ylabel("$\\Delta U$ (pp)"); a.set_title("(a)")

cc = load("conduit_calcite.csv"); cc["U"] = [E[(d, 1, bc)] for d, bc in zip(cc.dfn, cc.beta_cal)]
cc["rel"] = cc.U - cc.groupby("dfn").U.transform("max")
for bc, mk in ((-2, "o"), (-1, "D"), (0, "X"), (1, "s"), (2, "^")):
    d = cc[cc.beta_cal == bc]
    sc = b.scatter(d.conduit_calcite_VF0, d.rel, c=d.dfn.map(SIGMA), cmap="viridis", vmin=0.1, vmax=1.0, marker=mk, s=46, edgecolor="k", lw=0.4, zorder=3, label="$\\beta_{cal}$ = %d" % bc)
b.axvspan(1e-4, 2e-4, color="#1a7a3c", alpha=0.12, lw=0); b.axvline(2.5e-3, color="0.5", ls="--", lw=1)
b.axhline(0, color="k", lw=0.6); b.axhline(-5, color="0.6", ls=":", lw=0.7)
from matplotlib.lines import Line2D
b.legend(handles=[Line2D([], [], marker=mk, ls="", color="0.45", mec="k", label="$\\beta_{cal}$ = %d" % bc) for bc, mk in ((-2, "o"), (-1, "D"), (0, "X"), (1, "s"), (2, "^"))], fontsize=8.5, loc="lower left")
b.set_xscale("log"); b.set_xlim(3e-10, 6e-3); b.set_ylim(-72, 4); b.grid(alpha=0.25, which="both")
b.set_xlabel("conduit calcite VF$_0$"); b.set_ylabel("$U_{esc} - U_{esc}^{max}$ (pp)"); b.set_title("(b)")
fig.colorbar(sc, ax=b, pad=0.02).set_label("$\\sigma_{\\log|\\mathbf{q}|}$")
for lab, m in (("VF0 < 1e-4", cc.conduit_calcite_VF0 < 1e-4), ("1e-4 to 2e-4", cc.conduit_calcite_VF0.between(1e-4, 2e-4)), ("VF0 > 2e-4", cc.conduit_calcite_VF0 > 2e-4)):
    print("%-13s n = %2d, U_esc - U_esc_max from %.1f to %.1f pp" % (lab, m.sum(), cc.rel[m].min(), cc.rel[m].max()))
save(fig, "fig07_channeling")
