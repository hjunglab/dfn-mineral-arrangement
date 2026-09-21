#!/bin/bash
# Regenerate every figure and table that depends only on data/ (a few minutes on a laptop).
set -e; cd "$(dirname "$0")"
for s in fig03_escape_evolution fig05_heatmap fig06_random_fields fig07_channeling fig08_damkohler figS1_mineral_inventory tableS4_cutoff_scan; do
  echo "=== $s"; python "$s.py"
done
