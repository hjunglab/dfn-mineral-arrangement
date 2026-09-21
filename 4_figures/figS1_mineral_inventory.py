#!/usr/bin/env python3
"""Figure S1. Network-wide fraction of the initial pyrite and calcite remaining (mean +/- 1 SD over twenty networks)."""
from paper import *

inv = load("mineral_inventory.csv")
fig, ax = plt.subplots(2, 2, figsize=(11, 7.4), sharex=True)
for j, mag in enumerate((1, 2)):
    for label, (bp, bc), col in (("worst", (-1, mag), RED), ("homogeneous", (0, 0), GRAY), ("best", (1, -mag), BLUE)):
        d = inv[(inv.beta_pyr == bp) & (inv.beta_cal == bc)]
        for i, col_name in enumerate(("pyrite_remaining", "calcite_remaining")):
            P = d.pivot(index="t_yr", columns="dfn", values=col_name); m, sd = P.mean(axis=1), P.std(axis=1, ddof=1)
            ax[i][j].fill_between(P.index, m - sd, m + sd, color=col, alpha=0.18, lw=0); ax[i][j].plot(P.index, m, color=col, lw=2.2, label=label)
        if mag == 2: print("%-12s remaining at 3,000 yr: pyrite %.2f, calcite %.3f" % (label, d[d.t_yr == 3000].pyrite_remaining.mean(), d[d.t_yr == 3000].calcite_remaining.mean()))
    ax[0][j].set_title("$|\\beta_{cal}|$ = %d" % mag); ax[1][j].set_xlabel("Time (yr)")
    for i in range(2): ax[i][j].grid(alpha=0.25); ax[i][j].set_xlim(0, 3000); ax[i][j].text(0.97, 0.92, "(%s)" % "abcd"[2 * i + j], transform=ax[i][j].transAxes, ha="right", fontsize=13)
ax[0][0].set_ylabel("pyrite remaining (fraction of initial)"); ax[1][0].set_ylabel("calcite remaining (fraction of initial)")
for j in range(2): ax[0][j].set_ylim(0, 1.02); ax[1][j].set_ylim(0.975, 1.001)
ax[0][0].legend(fontsize=9, loc="lower left"); fig.tight_layout(); save(fig, "figS1_mineral_inventory")
