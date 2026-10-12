#!/usr/bin/env python3
"""Write one reactive transport deck with a prescribed pyrite and calcite arrangement.

    python make_decks.py <network_dir> <velocity_vtk> <deck_out> --beta-pyr B --beta-cal B
    python make_decks.py <network_dir> <velocity_vtk> <deck_out> --random-seed K [--sigma 1.0]

Flow-correlated arrangement (Equations 4-5 of the paper): fracture f receives the raw loading
w_f = (q_f / median(q))**beta, where q_f is the volume-averaged Darcy flux magnitude of the fracture.
beta > 0 concentrates the mineral in fast-flow fractures (aligned), beta < 0 in slow-flow fractures
(anti-aligned), beta = 0 is homogeneous. Random arrangement: w_f = exp(sigma * xi_f), xi_f ~ N(0, 1),
drawn independently for pyrite and calcite. In every case the loadings are normalized so that the
network-mean volume fraction of each mineral is unchanged.

The per-fracture volume fractions are written into the chemistry template as one
REGION / TRANSPORT_CONDITION / CONSTRAINT / INITIAL_CONDITION per fracture (PFLOTRAN reads the cell
ids of each fracture from <deck_dir>/regions.h5, created here).
"""
import os, re, sys, argparse
import numpy as np
import h5py

HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("network_dir"); ap.add_argument("velocity_vtk"); ap.add_argument("deck_out")
ap.add_argument("--beta-pyr", type=float, default=0.0); ap.add_argument("--beta-cal", type=float, default=0.0)
ap.add_argument("--random-seed", type=int, default=None, help="flow-uncorrelated random field instead of beta")
ap.add_argument("--sigma", type=float, default=1.0, help="ln-standard deviation of the random weights")
ap.add_argument("--template", default=os.path.join(HERE, "..", "1_inputs", "chemistry_template.in"))
ap.add_argument("--restart", default="flow-restart.h5", help="steady flow checkpoint, relative to the deck")
a = ap.parse_args()

MEAN_VF = {"Pyrite": 0.5e-5, "Calcite": 2.5e-3}                 # network-mean volume fractions (bulk basis)
# reactive surface area [m2 m-3 bulk] = BET [m2 g-1] x density [g cm-3] x 1e6 x volume fraction (Marty et al., 2015)
BET_RHO = {"Quartz": (0.05, 2.65), "K-Feldspar": (0.11, 2.56), "Albite": (0.08, 2.62), "Annite": (1.7, 3.0),
           "Ripidolite-14A": (0.0027, 2.9), "Calcite": (0.7, 2.71), "Pyrite": (0.05, 5.0)}
SILICATE_VF = {"Quartz": 0.29117, "K-Feldspar": 0.2768, "Albite": 0.1782, "Annite": 9.218e-5, "Ripidolite-14A": 3.746e-6}


def area(mineral, vf):
    bet, rho = BET_RHO[mineral]; return bet * (rho * 1e6) * vf


# ---- mesh: cell volumes and the fracture each cell belongs to
mat = [int(l.split()[0]) for l in open(os.path.join(a.network_dir, "materialid.dat")).read().splitlines()[3:]]
n = len(mat); vol = np.zeros(n)
with open(os.path.join(a.network_dir, "full_mesh_vol_area.uge")) as f:
    for _ in range(int(f.readline().split()[1])):
        t = f.readline().split(); vol[int(t[0]) - 1] = float(t[4])
cells = {}
for cid, fr in enumerate(mat, 1): cells.setdefault(fr, []).append(cid)
fracs = sorted(cells)
frac_vol = {fr: vol[np.array(cells[fr]) - 1].sum() for fr in fracs}
V = sum(frac_vol.values())


# ---- fracture-averaged Darcy flux magnitude q_f and its median
def vtk_scalar(name):
    txt = open(a.velocity_vtk).read(); i = txt.find("SCALARS " + name); k = txt.find("\n", txt.find("LOOKUP_TABLE", i))
    v = []
    for tok in txt[k + 1:].split():
        try: v.append(float(tok))
        except ValueError: break
    return np.array(v[:n])


speed = np.sqrt(vtk_scalar("Vlx") ** 2 + vtk_scalar("Vly") ** 2 + vtk_scalar("Vlz") ** 2)
q = {fr: (speed[np.array(cells[fr]) - 1] * vol[np.array(cells[fr]) - 1]).sum() / frac_vol[fr] for fr in fracs}
q_med = np.median([q[fr] for fr in fracs])


# ---- raw loadings and mass-conserving normalization
def loadings(mineral, beta, offset):
    if a.random_seed is not None:
        rng = np.random.default_rng(a.random_seed + offset)
        w = {fr: float(np.exp(a.sigma * rng.standard_normal())) for fr in fracs}
    else:
        w = {fr: (1.0 if beta == 0 else (q[fr] / max(q_med, 1e-30)) ** beta if beta > 0
                  else (1.0 / max(q[fr] / max(q_med, 1e-30), 1e-30)) ** (-beta)) for fr in fracs}
    s = sum(w[fr] * frac_vol[fr] for fr in fracs)
    return {fr: w[fr] * MEAN_VF[mineral] * V / s for fr in fracs}


pyr = loadings("Pyrite", a.beta_pyr, 0)
cal = loadings("Calcite", a.beta_cal, 1000)

# ---- deck blocks
SIL = "".join("    %-15s %s  %.4e\n" % (m, "%.5g" % SILICATE_VF[m], area(m, SILICATE_VF[m])) for m in SILICATE_VF)
WATER = """    H+        8.1       P
    Na+       9.583d-4  Z
    K+        8.696d-6  T
    Ca++      3.730d-4  T
    Mg++      5.143d-5  T
    Cl-       4.456d-5  T
    HCO3-     1.400d-3  T
    SO4--     7.006d-5  T
    F-        2.474d-4  T
    Al+++     1.0d-8    M  K-Feldspar
    Fe++      1.495d-7  T
    UO2++     1.798d-8  T
    O2(aq)    1.0d-12   T
    SiO2(aq)  3.011d-4  T
    U++++     1.0d-12   T
    Acetate-  1.0d-10   T"""                                  # initial reducing groundwater (DB1-7)


def minerals(fr):
    return (SIL + "    Pyrite          %.4e  %.4e\n" % (pyr[fr], area("Pyrite", pyr[fr]))
            + "    Calcite         %.4e  %.4e\n" % (cal[fr], area("Calcite", cal[fr]))
            + "    Uraninite       1.d-10    1.d4 m^2/kg\n    Ferrihydrite    0.0d0     1.0d0")


deck_dir = os.path.dirname(os.path.abspath(a.deck_out)); os.makedirs(deck_dir, exist_ok=True)
h5 = os.path.join(deck_dir, "regions.h5")
with h5py.File(h5, "w") as h:
    g = h.create_group("Regions")
    for fr in fracs: g.create_group("frac%d" % fr).create_dataset("Cell Ids", data=np.array(cells[fr], dtype=np.int32))

blocks = "\n".join(
    "REGION frac%d\n  FILE %s\nEND\nTRANSPORT_CONDITION tc_frac%d\n  TYPE dirichlet_zero_gradient\n"
    "  CONSTRAINT_LIST\n    0.d0 gw_frac%d\n  /\nEND\nCONSTRAINT gw_frac%d\n  CONCENTRATIONS\n%s\n  /\n  MINERALS\n%s\n  /\nEND"
    % (fr, h5, fr, fr, fr, WATER, minerals(fr)) for fr in fracs)
ics = "\n".join("INITIAL_CONDITION\n  FLOW_CONDITION pinned_p\n  TRANSPORT_CONDITION tc_frac%d\n  REGION frac%d\nEND" % (fr, fr) for fr in fracs)

t = open(a.template).read()
assert t.count("#@FRACTURE_BLOCKS@") == 1 and t.count("#@FRACTURE_INITIAL_CONDITIONS@") == 1
t = t.replace("#@FRACTURE_BLOCKS@", blocks).replace("#@FRACTURE_INITIAL_CONDITIONS@", ics)
t = re.sub(r"FILENAME flow_v05c-restart\.h5", "FILENAME " + a.restart, t)
open(a.deck_out, "w").write(t)

totP = sum(pyr[fr] * frac_vol[fr] for fr in fracs) / V; totC = sum(cal[fr] * frac_vol[fr] for fr in fracs) / V
print("%s: %d fractures | mean VF pyrite %.3e calcite %.3e | VF range pyrite %.2e-%.2e calcite %.2e-%.2e"
      % (os.path.basename(a.deck_out), len(fracs), totP, totC, min(pyr.values()), max(pyr.values()), min(cal.values()), max(cal.values())))
