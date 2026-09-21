# dfn-mineral-arrangement

Reactive transport simulations of uranium escape in discrete fracture networks with flow-correlated pyrite and calcite arrangements.

This repository accompanies

> Jung, H., Namgung, S., Kim, H.-r., Lee, S., & Kwon, J.-S. (2026). Flow-correlated mineral arrangement controls uranium escape in fractured media under oxidizing intrusion. Water Resources Research (submitted).

It contains (i) the inputs and scripts that generate the simulations, (ii) the quantities extracted from them, and
(iii) the scripts that turn those quantities into the figures and numbers of the paper.

## Reproduce the figures (no simulation needed)

```bash
conda env create -f environment.yml && conda activate dfn-mineral-arrangement
bash 4_figures/make_all.sh          # a few minutes; writes 4_figures/output/
```

Each script prints the numbers quoted in the paper (for example the constants of Equation 10, the validation
statistics, and the logistic midpoint of the calcite transition) next to the figure it draws.

| Paper item | Script | Input |
|---|---|---|
| Figure 3 (U escape over time) | `4_figures/fig03_escape_evolution.py` | `data/escape_timeseries.csv` |
| Figure 5 (arrangement heat map) | `4_figures/fig05_heatmap.py` | `data/escape_timeseries.csv` |
| Figure 6 (random mineral fields) | `4_figures/fig06_random_fields.py` | `data/escape_timeseries.csv` |
| Figure 7 (channeling and calcite threshold) | `4_figures/fig07_channeling.py` | `data/escape_timeseries.csv`, `data/conduit_terms.csv`, `data/conduit_calcite.csv` |
| Figure 8 (Damköhler relationship) | `4_figures/fig08_damkohler.py` | `data/conduit_terms.csv` |
| Figure S1 (mineral inventories) | `4_figures/figS1_mineral_inventory.py` | `data/mineral_inventory.csv` |
| Table S4 (conduit cutoff) | `4_figures/tableS4_cutoff_scan.py` | `data/cutoff_scan.csv` |
| Table 2 (network properties) | `3_extract/network_properties.py` | raw output |
| Figures 1, 2, 4 (3-D views) | `4_figures/fig01_02_04_3d_views.py` | raw output, `pyvista` |

Figures 1, 2 and 4 and Table 2 need the raw simulation output (meshes and VTK snapshots), which is too large to
distribute here; regenerate it with the scripts in `2_simulation/`.

## Repository layout (the folders follow the workflow)

```
1_inputs/       fracture statistics, PFLOTRAN flow deck and chemistry template, thermodynamic database, list of networks
2_simulation/   network generation, deck generation with a prescribed mineral arrangement, run scripts for every series
3_extract/      raw PFLOTRAN output -> the CSV files in data/
data/           extracted results used by the figure scripts (see data/README.md for the column definitions)
4_figures/      figures and tables of the paper
```

## Re-run the simulations

The simulations were run with dfnWorks v2.9.6 and the PFLOTRAN build shipped in its Docker image
(`ees16/dfnworks:v2.9.6`). Inside the container, with this repository mounted:

```bash
cd 2_simulation
python3 build_network.py 500          # one network per seed (seeds: column 2 of 1_inputs/networks.csv)
bash run_flow.sh 500                  # steady flow field, Q = 5e-7 m3/d
bash run_factorial.sh 500             # 3 pyrite x 5 calcite flow-correlated arrangements, 3,000 yr each
bash run_random.sh                    # Series A1: flow-uncorrelated random mineral fields
bash run_flowrate.sh                  # Series A2: flow rate x3 and x0.3, with and without O2 compensation
bash run_sensitivity.sh               # Series A3: one-at-a-time parameter sensitivity
cd ../3_extract && for s in extract_escape extract_conduit_terms extract_conduit_calcite extract_inventory; do python3 $s.py; done
```

Output goes to `work/` (override with `WORK=/path`). One reactive transport run takes roughly 15-90 minutes on four
cores; the full set is 457 runs. `env.sh` holds the PFLOTRAN path, MPI ranks per run (`NP`) and the number of
concurrent runs (`NPAR`). An arrangement is written as `pyr<beta_pyr>_cal<beta_cal>`, e.g. `pyr+1_cal-2` is the best
case (pyrite aligned, calcite anti-aligned at |beta_cal| = 2), `pyr+0_cal+0` the homogeneous reference and
`pyr-1_cal+2` the worst case.

The cleaned scripts were checked against the production runs of the paper: `build_network.py` reproduces the
production network for a given seed, `run_flow.sh` reproduces the velocity field exactly, `make_decks.py`
reproduces the production decks line by line, and the extraction scripts reproduce the archived values.

## Notes on the inputs

* `1_inputs/hanford_neaUOH4.dat` is the `hanford.dat` thermodynamic database distributed with PFLOTRAN with two
  entries changed: log K (25 °C) of U(OH)4(aq) from 4.57 to 10.0 (NEA value; Grenthe et al., 2021) and the molar
  volume of UO2(am) (not used in these simulations).
* Mineral rate constants, reaction orders and BET surface areas follow Marty et al. (2015); groundwater
  compositions are from Kwon et al. (2026); fracture statistics are from the KURT borehole survey (Namgung et al., 2025).

## License

Code: MIT (`LICENSE`). Data in `data/` and `1_inputs/networks.csv`: CC BY 4.0 (`data/LICENSE`). The thermodynamic
database remains under the terms of the PFLOTRAN distribution.
