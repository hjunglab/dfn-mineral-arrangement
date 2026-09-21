"""Shared helpers for the extraction scripts: directory layout, PFLOTRAN readers, conduit definition.

The scripts read the raw simulation output under $WORK (default: ../work), laid out by 2_simulation/:
    networks/s<seed>/                                   mesh and grid files
    flow/s<seed>[_Q<factor>]/flow-vel-001.vtk           steady Darcy velocity
    rxn/factorial/s<seed>/pyr<+b>_cal<+b>/              deck.in, deck-mas.dat, deck-000.vtk ... deck-030.vtk
    rxn/random/s<seed>/field<k>/
    rxn/flowrate/s<seed>/Q<factor>_<O2fixed|O2scaled>/pyr<+b>_cal<+b>/
    rxn/sensitivity/s<seed>/<name>_<value>/pyr<+b>_cal<+b>/
"""
import os, re, csv, glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("WORK", os.path.join(HERE, "..", "work"))
DATA = os.path.join(HERE, "..", "data")
NETWORKS = list(csv.DictReader(open(os.path.join(HERE, "..", "1_inputs", "networks.csv"))))
DFN_OF = {int(r["seed"]): int(r["dfn"]) for r in NETWORKS}
K_O2_DEFAULT, KB_DEFAULT = 1e-6, 3e-7

_FIX = re.compile(r"(\d)([-+]\d{3})$")          # Fortran drops the 'E' of three-digit exponents: 2.5-132


def _float(tok):
    try: return float(tok)
    except ValueError: return float(_FIX.sub(r"\1e\2", tok))


def read_vtk_field(path, name, n):
    """Cell values of scalar <name> from an ASCII PFLOTRAN VTK file."""
    txt = open(path).read()
    m = re.search(r"SCALARS %s\s+\w+\s+\d+\s*\nLOOKUP_TABLE \w+\s*\n" % re.escape(name), txt)
    return np.array([_float(t) for t in txt[m.end():].split()[:n]])


def read_escape(mas_file):
    """Time [yr] and instantaneous U escape [%] = total U flux out / total U flux in (U(VI) + U(IV))."""
    names = re.findall(r'"([^"]*)"', open(mas_file).readline()); i = {k: j for j, k in enumerate(names)}
    A = np.array([l.split() for l in open(mas_file).read().strip().split("\n")[1:]], float)
    out = A[:, i["outflow UO2++ [mol/y]"]] + A[:, i["outflow U++++ [mol/y]"]]
    inn = A[:, i["inflow UO2++ [mol/y]"]] + A[:, i["inflow U++++ [mol/y]"]]
    with np.errstate(divide="ignore", invalid="ignore"):
        return A[:, 0], np.where(inn > 0, -out / inn * 100, np.nan)


def cell_volumes(seed):
    with open(os.path.join(WORK, "networks", "s%d" % seed, "full_mesh_vol_area.uge")) as f:
        n = int(f.readline().split()[1]); vol = np.zeros(n)
        for _ in range(n):
            t = f.readline().split(); vol[int(t[0]) - 1] = float(t[4])
    return n, vol


def darcy_speed(seed, n, flow_dir=None):
    vtk = os.path.join(flow_dir or os.path.join(WORK, "flow", "s%d" % seed), "flow-vel-001.vtk")
    return np.sqrt(sum(read_vtk_field(vtk, c, n) ** 2 for c in ("Vlx", "Vly", "Vlz")))


def conduit_mask(speed, vol, cutoff=0.90):
    """Conduit = smallest set of cells that carries <cutoff> of the network-wide sum of |q|_i V_i."""
    w = speed * vol; order = np.argsort(-w); cum = np.cumsum(w[order]) / w.sum()
    mask = np.zeros(len(w), bool); mask[order[:np.searchsorted(cum, cutoff) + 1]] = True
    return mask


def deck_parameters(deck):
    """Microbial rate constant k_b and O2 inhibition constant K_O2 actually used in a deck."""
    txt = open(deck).read(); blk = txt[re.search(r"^\s*MICROBIAL_REACTION", txt, re.M).start():]
    num = lambda s: float(s.replace("d", "e").replace("D", "e"))
    kb = num(re.search(r"RATE_CONSTANT\s+([0-9.]+[dDeE][-+]?\d+)", blk).group(1))
    ko2 = num(re.search(r"SPECIES_NAME O2\(aq\)[\s\S]*?THRESHOLD_CONCENTRATION\s+([0-9.]+[dDeE][-+]?\d+)", blk).group(1))
    return kb, ko2


def _betas(label):
    m = re.match(r"pyr([-+]\d+)_cal([-+]\d+)$", label); return int(m.group(1)), int(m.group(2))


def completed(run_dir, last="deck-030.vtk"):
    return os.path.exists(os.path.join(run_dir, "deck-mas.dat")) and os.path.exists(os.path.join(run_dir, last))


def runs():
    """All completed runs as dicts: series, dfn, seed, beta_pyr, beta_cal, field, condition, q_factor, run_dir, flow_dir."""
    out = []
    def add(series, seed, run_dir, bp="", bc="", field="", condition="", qf=1.0, flow_dir=None):
        if completed(run_dir):
            out.append(dict(series=series, dfn=DFN_OF[seed], seed=seed, beta_pyr=bp, beta_cal=bc, field=field, condition=condition,
                            q_factor=qf, run_dir=run_dir, flow_dir=flow_dir or os.path.join(WORK, "flow", "s%d" % seed)))
    for seed in sorted(DFN_OF, key=DFN_OF.get):
        base = os.path.join(WORK, "rxn")
        for d in sorted(glob.glob(os.path.join(base, "factorial", "s%d" % seed, "pyr*_cal*"))):
            bp, bc = _betas(os.path.basename(d)); add("factorial", seed, d, bp, bc)
        for d in sorted(glob.glob(os.path.join(base, "random", "s%d" % seed, "field*"))):
            add("random", seed, d, field=int(os.path.basename(d)[5:]))
        for d in sorted(glob.glob(os.path.join(base, "flowrate", "s%d" % seed, "Q*", "pyr*_cal*"))):
            cond = os.path.basename(os.path.dirname(d)); qf = float(cond.split("_")[0][1:]); bp, bc = _betas(os.path.basename(d))
            add("flowrate", seed, d, bp, bc, condition=cond, qf=qf, flow_dir=os.path.join(WORK, "flow", "s%d_Q%s" % (seed, cond.split("_")[0][1:])))
        for d in sorted(glob.glob(os.path.join(base, "sensitivity", "s%d" % seed, "*", "pyr*_cal*"))):
            bp, bc = _betas(os.path.basename(d)); add("sensitivity", seed, d, bp, bc, condition=os.path.basename(os.path.dirname(d)))
    return out


KEYS = ["series", "dfn", "beta_pyr", "beta_cal", "field", "condition"]
