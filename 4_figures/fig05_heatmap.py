#!/usr/bin/env python3
"""Figure 5. Ensemble-mean U escape at 1,500 yr for the 3 x 3 pyrite-calcite arrangements at (a) |beta_cal| = 1, (b) |beta_cal| = 2."""
from paper import *
from matplotlib.patches import Rectangle

E = escape_at(load("escape_timeseries.csv"), 1500)
LAB = ["aligned", "homogeneous", "anti-aligned"]
mats = {mag: np.array([[[E[(d, bp, bc)] for d in NETWORKS.index] for bc in (-mag, 0, mag)] for bp in (1, 0, -1)]) for mag in (1, 2)}
lo = min(m.mean(2).min() for m in mats.values()); hi = max(m.mean(2).max() for m in mats.values())
norm, cmap = plt.Normalize(lo, hi), plt.get_cmap("RdBu_r")
fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.0), gridspec_kw=dict(wspace=0.32))
for a, mag, letter in zip(ax, (1, 2), "ab"):
    M, S = mats[mag].mean(2), mats[mag].std(2, ddof=1)
    for i in range(3):
        for j in range(3):                                   # vector rectangles (an image would be rasterized in the PDF)
            a.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor=cmap(norm(M[i, j])), edgecolor="none"))
            tc = "white" if not 0.28 < norm(M[i, j]) < 0.72 else "black"
            a.text(j, i - 0.12, "%.1f" % M[i, j], ha="center", va="center", fontsize=12, fontweight="bold", color=tc)
            a.text(j, i + 0.18, "(± %.1f)" % S[i, j], ha="center", va="center", fontsize=9, color=tc)
    a.set_xlim(-0.5, 2.5); a.set_ylim(2.5, -0.5); a.set_xticks(range(3)); a.set_yticks(range(3))
    a.set_xticklabels(["anti-aligned", "homogeneous", "aligned"]); a.set_yticklabels(LAB if mag == 1 else [""] * 3)
    a.set_xlabel("Calcite arrangement"); a.set_title("(%s) $|\\beta_{cal}|$ = %d" % (letter, mag))
    print("|beta_cal| = %d, rows pyrite aligned/homogeneous/anti, columns calcite anti/homogeneous/aligned:\n" % mag, np.round(M, 1))
ax[0].set_ylabel("Pyrite arrangement")
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
cb = fig.colorbar(sm, ax=ax, shrink=0.85, pad=0.02); cb.set_label("$\\langle U_{esc} \\rangle$ at 1,500 yr (%)"); cb.solids.set_edgecolor("face")
save(fig, "fig05_heatmap")
