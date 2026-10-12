# Data dictionary

All files were produced by the scripts in `3_extract/` from the raw PFLOTRAN output of the 457 simulations (300 factorial +
20 random-field + 48 flow-rate + 89 sensitivity runs; one sensitivity configuration was numerically intractable and is absent).
The simulations use the near-equilibrium uraninite treatment described in the top-level README.

**Keys shared by the files**

| Column | Meaning |
|---|---|
| `series` | `factorial` (flow-correlated arrangements, 300 runs), `random` (Series A1, 20), `flowrate` (Series A2, 48), `sensitivity` (Series A3, 89) |
| `dfn` | network number, 1-20 in order of increasing flow channeling (`1_inputs/networks.csv`) |
| `beta_pyr`, `beta_cal` | flow-alignment exponents of pyrite (+1, 0, -1) and calcite (+2, +1, 0, -1, -2); empty for random fields. Best = (+1, -2), homogeneous = (0, 0), worst = (-1, +2) |
| `field` | index of the random mineral field (1-5); empty otherwise |
| `condition` | flow-rate series: `Q3_O2fixed`, `Q3_O2scaled`, `Q0.3_O2fixed`, `Q0.3_O2scaled`; sensitivity series: `<parameter>_<value>` with parameters `kb`, `KO2`, `O2`, `kcal`, `kpyr` |
| `t_yr` | time since the onset of the oxidizing intrusion (yr) |
| `U_esc_percent` | instantaneous U escape = total U flux out / total U flux in x 100 |

**Files**

| File | Content |
|---|---|
| `escape_timeseries.csv` | U escape of every run at 10-yr resolution |
| `conduit_terms.csv` | per 100-yr snapshot: `q_factor` (flow-rate multiplier), `kb`, `K_O2` (values used in the run), and the conduit volume averages `I_O2` = <K_O2/(K_O2+[O2])>, `inv_alpha` = <[UO2++]/[U(VI)]total>, `I_over_alpha` = <I_O2/alpha>. Conduit = smallest set of cells carrying 90 % of the network-wide sum of \|q\|V |
| `cutoff_scan.csv` | `I_over_alpha_<c>` and conduit volume fraction `Vc_over_V_<c>` for conduit cutoffs c = 50-100 % (Table S4) |
| `conduit_calcite.csv` | initial conduit and network-mean calcite volume fraction for the five calcite arrangements under aligned pyrite (Figure 7b) |
| `mineral_inventory.csv` | network-wide fraction of initial pyrite and calcite remaining (Figure S1) |
