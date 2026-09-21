#!/usr/bin/env python3
"""Channeling index, transit time and simulated P10 of each network (Tables 1 and 2).

sigma_log_q : standard deviation over fractures of log10 of the volume-averaged Darcy flux magnitude q_f
tau         : mean hydraulic transit time = porosity x fracture volume / Q
P10         : fractures intersected per metre along 2,100 vertical scanlines, by family (compare with measured P10)
"""
import os, csv
import numpy as np
from common import WORK, NETWORKS, cell_volumes, darcy_speed

PHI, Q = 0.25, 5e-7 * 365.25            # porosity, flow rate [m3 yr-1]
LX, LY, LZ = 6.0, 3.5, 2.0


def simulated_p10(seed):
    g = os.path.join(WORK, "networks", "s%d" % seed, "dfnGen_output")
    rad = np.array([[float(a[0]), int(a[2])] for a in (l.split() for l in open(g + "/radii_Final.dat").read().splitlines()[2:] if l.strip())])
    nv = np.array([[float(x) for x in l.split()] for l in open(g + "/normal_vectors.dat").read().splitlines() if l.strip()])
    tr = np.array([[float(x) for x in p[:3]] for p in (l.split() for l in open(g + "/translations.dat").read().splitlines()[1:]) if len(p) >= 3 and (len(p) == 3 or p[3] != "R")])
    hits = np.zeros(4); nlines = 0
    for x0 in np.linspace(-LX / 2 * 0.98, LX / 2 * 0.98, 60):
        for y0 in np.linspace(-LY / 2 * 0.98, LY / 2 * 0.98, 35):
            nlines += 1
            for (r, fam), nrm, c in zip(rad, nv, tr):
                if abs(nrm[2]) < 1e-9: continue
                z = (np.dot(nrm, c) - nrm[0] * x0 - nrm[1] * y0) / nrm[2]
                if abs(z) <= LZ / 2 + 1e-9 and np.linalg.norm(np.array([x0, y0, z]) - c) <= r: hits[int(fam)] += 1
    return hits[1:] / nlines / LZ


print("dfn,seed,n_fractures,sigma_log_q,tau_yr,P10_A,P10_B,P10_C")
for row in NETWORKS:
    seed = int(row["seed"]); n, vol = cell_volumes(seed); speed = darcy_speed(seed, n)
    mat = np.array([int(l.split()[0]) for l in open(os.path.join(WORK, "networks", "s%d" % seed, "materialid.dat")).read().splitlines()[3:]][:n])
    qf = np.array([(speed[mat == f] * vol[mat == f]).sum() / vol[mat == f].sum() for f in sorted(set(mat))])
    p10 = simulated_p10(seed)
    print("%s,%d,%d,%.3f,%.2f,%.3f,%.3f,%.3f" % (row["dfn"], seed, len(qf), np.log10(qf[qf > 0]).std(), PHI * vol.sum() / Q, *p10))
