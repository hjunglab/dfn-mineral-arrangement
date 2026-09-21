#!/usr/bin/env python3
"""Initial conduit calcite volume fraction of the five calcite arrangements under aligned pyrite (Figure 7b)
-> data/conduit_calcite.csv"""
import os, csv
import numpy as np
from common import runs, read_vtk_field, cell_volumes, darcy_speed, conduit_mask, DATA

with open(os.path.join(DATA, "conduit_calcite.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["dfn", "beta_cal", "conduit_calcite_VF0", "network_calcite_VF0"])
    for r in runs():
        if r["series"] == "factorial" and r["beta_pyr"] == 1:
            n, vol = cell_volumes(r["seed"]); m = conduit_mask(darcy_speed(r["seed"], n), vol)
            c0 = np.clip(read_vtk_field(os.path.join(r["run_dir"], "deck-000.vtk"), "Calcite_VF", n), 0, None)
            w.writerow([r["dfn"], r["beta_cal"], "%.6e" % ((vol[m] * c0[m]).sum() / vol[m].sum()), "%.6e" % ((vol * c0).sum() / vol.sum())])
print("wrote conduit_calcite.csv")
