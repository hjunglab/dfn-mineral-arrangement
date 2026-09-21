#!/usr/bin/env python3
"""Conduit-averaged terms of the Damkohler number at every 100-yr snapshot of every run.

For the conduit cells (common.conduit_mask) and each snapshot:
    I_O2         = < K_O2 / (K_O2 + [O2]) >            volume-weighted mean oxygen-inhibition factor
    inv_alpha    = < [UO2++] / [U(VI)]_total >         volume-weighted mean free-uranyl fraction (1/alpha)
    I_over_alpha = < I_O2 / alpha >                    the combined term entering lambda_eff
Output: data/conduit_terms.csv (90 % cutoff, used in the paper) and data/cutoff_scan.csv (I_over_alpha for
cutoffs 50-100 %, Table S4).        Usage: python extract_conduit_terms.py [n_processes=4]
"""
import os, sys, csv, glob
import numpy as np
from multiprocessing import Pool
from common import runs, read_escape, read_vtk_field, cell_volumes, darcy_speed, conduit_mask, deck_parameters, DATA, KEYS

CUTOFFS = [0.50, 0.60, 0.75, 0.80, 0.85, 0.90, 0.95, 0.98, 1.00]


def one(r):
    n, vol = cell_volumes(r["seed"]); speed = darcy_speed(r["seed"], n, r["flow_dir"])
    masks = {c: conduit_mask(speed, vol, c) for c in CUTOFFS}
    kb, ko2 = deck_parameters(os.path.join(r["run_dir"], "deck.in"))
    t, e = read_escape(os.path.join(r["run_dir"], "deck-mas.dat"))
    rows = []
    for k, vtk in enumerate(sorted(glob.glob(os.path.join(r["run_dir"], "deck-0*.vtk")))):
        o2 = np.clip(read_vtk_field(vtk, "Total_O2(aq)", n), 0, None)
        tot = np.clip(read_vtk_field(vtk, "Total_UO2++", n), 1e-300, None); free = np.clip(read_vtk_field(vtk, "Free_UO2++", n), 1e-300, None)
        I, f = ko2 / (ko2 + o2), free / tot
        avg = lambda x, m: float((vol[m] * x[m]).sum() / vol[m].sum())
        m90 = masks[0.90]
        rows.append([r[x] for x in KEYS] + [r["q_factor"], 100 * k, "%.4f" % np.interp(100 * k, t, e), kb, ko2,
                    avg(I, m90), avg(f, m90), avg(I * f, m90)] + [avg(I * f, masks[c]) for c in CUTOFFS]
                    + [float(vol[masks[c]].sum() / vol.sum()) for c in CUTOFFS])
    return rows


if __name__ == "__main__":
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    with Pool(nproc) as p: res = p.map(one, runs(), chunksize=1)
    head = KEYS + ["q_factor", "t_yr", "U_esc_percent", "kb", "K_O2"]
    with open(os.path.join(DATA, "conduit_terms.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(head + ["I_O2", "inv_alpha", "I_over_alpha"])
        for rows in res:
            for x in rows: w.writerow(x[:len(head) + 3])
    with open(os.path.join(DATA, "cutoff_scan.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(head + ["I_over_alpha_%d" % round(100 * c) for c in CUTOFFS] + ["Vc_over_V_%d" % round(100 * c) for c in CUTOFFS])
        for rows in res:
            for x in rows:
                if x[0] in ("factorial", "random", "flowrate") and (x[0] != "factorial" or abs(int(x[3])) != 1):   # series used in Table S4
                    w.writerow(x[:len(head)] + x[len(head) + 3:])
    print("wrote conduit_terms.csv and cutoff_scan.csv (%d runs)" % len(res))
