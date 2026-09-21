#!/usr/bin/env python3
"""Table S4. Sensitivity of Equation 10 to the flux cutoff that defines the conduit cells."""
from paper import *

scan = load("cutoff_scan.csv"); scan = scan[scan.t_yr >= 100].reset_index(drop=True); E = scan.U_esc_percent.values
fitset, a1, a2 = (scan.series == "factorial").values, (scan.series == "random").values, (scan.series == "flowrate").values
rows = []
for c in (50, 60, 75, 80, 85, 90, 95, 98, 100):
    Da = damkohler(scan, "I_over_alpha_%d" % c); p = fit_escape_law(Da[fitset], E[fitset]); pred = escape_law(np.log10(Da), *p)
    rows.append([c, 10 ** p[0], p[1], r_squared(E[fitset], pred[fitset]), np.sqrt(((E - pred)[fitset] ** 2).mean()),
                 r_squared(E[a1], pred[a1]), r_squared(E[a2], pred[a2]), scan.loc[fitset, "Vc_over_V_%d" % c].mean()])
t = pd.DataFrame(rows, columns=["cutoff_percent", "m", "n", "R2_ensemble", "RMSE_pp", "R2_A1_random", "R2_A2_flowrate", "Vc_over_V"])
t.to_csv(os.path.join(OUT, "tableS4_cutoff_scan.csv"), index=False, float_format="%.3f"); print(t.round(3).to_string(index=False))
