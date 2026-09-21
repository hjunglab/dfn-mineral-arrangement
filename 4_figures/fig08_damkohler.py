#!/usr/bin/env python3
"""Figure 8. (a) U escape of the 180 strong-correlation simulations versus Da with the fitted Equation 10;
(b) simulations not used in the fit, evaluated with the same constants."""
from paper import *

terms = load("conduit_terms.csv"); terms = terms[terms.t_yr >= 100].reset_index(drop=True)
Da, E = damkohler(terms), terms.U_esc_percent.values
fitset = is_fitting_set(terms).values
p = fit_escape_law(Da[fitset], E[fitset]); m, n = 10 ** p[0], p[1]
res = E - escape_law(np.log10(Da), *p)
print("Equation 10 fitted to %d points: m = %.2f, n = %.2f, R2 = %.3f" % (fitset.sum(), m, n, r_squared(E[fitset], E[fitset] - res[fitset])))
print("  U_esc = 50%% at Da = %.2f, 10%% at Da = %.2f" % (m * np.log(2) ** (1 / n), m * np.log(10) ** (1 / n)))
VALIDATION = [("weak correlation, $|\\beta_{cal}|$ = 1", (terms.series == "factorial") & (terms.beta_cal.abs() == 1), "#ff7f0e", "s"),
              ("A1: random mineral fields", terms.series == "random", "#2ca02c", "*"),
              ("A2: flow rate", terms.series == "flowrate", "#d62728", "o"),
              ("A3: parameter sensitivity", terms.series == "sensitivity", "#9467bd", "D")]
for lab, mask, _, _ in VALIDATION:
    k = mask.values; print("  %-40s n = %5d  R2 = %.3f  residual %+.1f +/- %.1f pp" % (lab.replace("$", ""), k.sum(), r_squared(E[k], E[k] - res[k]), res[k].mean(), res[k].std()))
pa = fit_escape_law(Da, E); print("refit to all %d points: m = %.2f, n = %.2f, R2 = %.3f" % (len(E), 10 ** pa[0], pa[1], r_squared(E, escape_law(np.log10(Da), *pa))))
big = np.abs(res) >= 10; print("residuals >= 10 pp: %d of %d, %d positive, all in calcite anti-aligned runs: %s" % (big.sum(), len(E), (res[big] > 0).sum(), bool((terms.beta_cal[big] < 0).all())))

fig, (a, b) = plt.subplots(1, 2, figsize=(12.5, 4.8), sharey=True, gridspec_kw=dict(wspace=0.08))
xx = np.logspace(np.log10(Da.min()) - 0.15, np.log10(Da.max()) + 0.15, 400)
for bp, mk in ((1, "^"), (0, "o"), (-1, "v")):
    for bc, col in ((2, BLUE), (0, GRAY), (-2, RED)):
        k = fitset & (terms.beta_pyr == bp).values & (terms.beta_cal == bc).values
        a.scatter(Da[k], E[k], s=13, marker=mk, color=col, alpha=0.55, lw=0, label="pyrite %+d / calcite %+d" % (bp, bc))
a.plot(xx, escape_law(np.log10(xx), *p), "k--", lw=1.4, label="Equation 10: m = %.2f, n = %.2f" % (m, n))
b.scatter(Da[fitset], E[fitset], s=5, color="0.8", alpha=0.5, lw=0)
for lab, mask, col, mk in VALIDATION:
    k = mask.values; b.scatter(Da[k], E[k], s=16 if mk == "*" else 11, marker=mk, color=col, alpha=0.6, lw=0, label=lab)
b.plot(xx, escape_law(np.log10(xx), *p), "k--", lw=1.4, label="Equation 10 (not refitted)")
for ax_, t in ((a, "(a)"), (b, "(b)")):
    ax_.set_xscale("log"); ax_.set_xlim(xx[0], xx[-1]); ax_.set_ylim(-2, 102); ax_.grid(alpha=0.25, which="both"); ax_.set_xlabel("$Da$"); ax_.set_title(t); ax_.legend(fontsize=7.6, loc="upper right")
a.set_ylabel("$U_{esc}$ (%)")
save(fig, "fig08_damkohler")
