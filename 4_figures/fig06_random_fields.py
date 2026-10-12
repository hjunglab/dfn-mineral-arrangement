#!/usr/bin/env python3
"""Figure 6. U escape at 1,500 yr for the flow-uncorrelated random mineral fields (Series A1) compared with the
homogeneous, best and worst arrangements. A short black bar marks the mean of the five random fields."""
from paper import *

ts = load("escape_timeseries.csv"); E = escape_at(ts, 1500)
rnd = ts[(ts.series == "random") & (ts.t_yr == 1500)]
nets = sorted(rnd.dfn.unique()); rng = np.random.default_rng(1)
fig, a = plt.subplots(figsize=(6.2, 4.5))
for k, dfn in enumerate(nets):
    v = rnd[rnd.dfn == dfn].U_esc_percent.values; first = k == 0
    a.scatter(k + rng.uniform(-0.09, 0.09, len(v)), v, s=32, color=plt.cm.viridis((SIGMA[dfn] - 0.1) / 0.9), edgecolor="k", lw=0.4, zorder=4, label="random fields" if first else None)
    a.hlines(v.mean(), k - 0.22, k + 0.22, color="k", lw=1.8, zorder=5, label="random fields, mean" if first else None)   # mean of the five random fields
    a.plot(k, E[(dfn,) + HOMOGENEOUS], "s", ms=8, mfc="none", mew=1.6, color=GRAY, label="homogeneous" if first else None)
    a.plot(k, E[(dfn,) + BEST], "v", ms=8, color=BLUE, label="best" if first else None)
    a.plot(k, E[(dfn,) + WORST], "^", ms=8, color=RED, label="worst" if first else None)
    print("DFN%-2d random: range %.0f pp, mean - homogeneous %+.1f pp | best-worst range %.0f pp" % (dfn, np.ptp(v), v.mean() - E[(dfn,) + HOMOGENEOUS], E[(dfn,) + WORST] - E[(dfn,) + BEST]))
a.set_xticks(range(len(nets))); a.set_xticklabels(["DFN%d\n$\\sigma$ = %.2f" % (d, SIGMA[d]) for d in nets])
a.set_ylim(-2, 100); a.grid(alpha=0.25, axis="y"); a.set_ylabel("$U_{esc}$ at 1,500 yr (%)"); a.legend(fontsize=8.5, loc="lower left")
save(fig, "fig06_random_fields")
