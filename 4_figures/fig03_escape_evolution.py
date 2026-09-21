#!/usr/bin/env python3
"""Figure 3. U escape over 3,000 yr for the best, homogeneous and worst arrangements: (a) |beta_cal| = 1,
(b) |beta_cal| = 2 (mean +/- 1 SD over twenty networks, thin lines = networks), (c) best case in four networks."""
from paper import *

ts = load("escape_timeseries.csv"); ts = ts[ts.series == "factorial"]
fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.4), sharey=True)
for a, mag, letter in zip(ax, (1, 2), "ab"):
    for label, (bp, bc), col in (("worst", (-1, mag), RED), ("homogeneous", (0, 0), GRAY), ("best", (1, -mag), BLUE)):
        d = ts[(ts.beta_pyr == bp) & (ts.beta_cal == bc)].pivot(index="t_yr", columns="dfn", values="U_esc_percent")
        a.plot(d.index, d.values, color=col, lw=0.45, alpha=0.22)
        m, sd = d.mean(axis=1), d.std(axis=1, ddof=1)
        a.fill_between(d.index, m - sd, m + sd, color=col, alpha=0.17, lw=0); a.plot(d.index, m, color=col, lw=2.2, label=label)
        print("|beta_cal|=%d %-12s <U_esc> at 100/1500/3000 yr: %.1f / %.1f / %.1f" % (mag, label, m[100], m[1500], m[3000]))
    a.axvline(1500, color="0.6", lw=0.8, ls=":"); a.set_title("(%s) $|\\beta_{cal}|$ = %d" % (letter, mag)); a.set_xlabel("Time (yr)")
cmap = plt.cm.viridis
for dfn in (6, 9, 15, 19):
    d = ts[(ts.dfn == dfn) & (ts.beta_pyr == 1) & (ts.beta_cal == -2)]
    ax[2].plot(d.t_yr, d.U_esc_percent, lw=2.2, color=cmap((SIGMA[dfn] - 0.1) / 0.9), label="DFN%d" % dfn)
ax[2].set_title("(c) best, $|\\beta_{cal}|$ = 2"); ax[2].set_xlabel("Time (yr)"); ax[2].legend(fontsize=9)
ax[0].set_ylabel("$U_{esc}$ (%)"); ax[0].legend(fontsize=9, loc="lower right")
for a in ax: a.set_xlim(0, 3000); a.set_ylim(0, 100); a.grid(alpha=0.25)
fig.tight_layout(); save(fig, "fig03_escape_evolution")
