#!/usr/bin/env python3
"""U escape time series of every run (10-yr resolution) -> data/escape_timeseries.csv"""
import os, csv
import numpy as np
from common import runs, read_escape, DATA, KEYS

grid = np.arange(10, 3001, 10.0)
with open(os.path.join(DATA, "escape_timeseries.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(KEYS + ["t_yr", "U_esc_percent"])
    for r in runs():
        t, e = read_escape(os.path.join(r["run_dir"], "deck-mas.dat"))
        for tt, ee in zip(grid, np.interp(grid, t, e)):
            w.writerow([r[k] for k in KEYS] + [int(tt), "%.4f" % ee])
print("wrote escape_timeseries.csv")
