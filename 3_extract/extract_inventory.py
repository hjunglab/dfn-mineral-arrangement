#!/usr/bin/env python3
"""Network-wide fraction of the initial pyrite and calcite remaining (Figure S1) -> data/mineral_inventory.csv
Best (pyr+1, cal-|b|), homogeneous and worst (pyr-1, cal+|b|) arrangements at |beta_cal| = 1 and 2."""
import os, csv, glob
import numpy as np
from common import runs, read_vtk_field, cell_volumes, DATA

WANTED = {(1, -2), (1, -1), (0, 0), (-1, 2), (-1, 1)}
with open(os.path.join(DATA, "mineral_inventory.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["dfn", "beta_pyr", "beta_cal", "t_yr", "pyrite_remaining", "calcite_remaining"])
    for r in runs():
        if r["series"] == "factorial" and (r["beta_pyr"], r["beta_cal"]) in WANTED:
            n, vol = cell_volumes(r["seed"]); p0 = c0 = None
            for k, vtk in enumerate(sorted(glob.glob(os.path.join(r["run_dir"], "deck-0*.vtk")))):
                p = (read_vtk_field(vtk, "Pyrite_VF", n) * vol).sum(); c = (read_vtk_field(vtk, "Calcite_VF", n) * vol).sum()
                if k == 0: p0, c0 = p, c
                w.writerow([r["dfn"], r["beta_pyr"], r["beta_cal"], 100 * k, "%.6f" % (p / p0), "%.6f" % (c / c0)])
print("wrote mineral_inventory.csv")
