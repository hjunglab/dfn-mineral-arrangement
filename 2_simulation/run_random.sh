#!/bin/bash
# Series A1: flow-uncorrelated random mineral fields (sigma = 1), five fields on each of four networks.
#   bash run_random.sh [seeds="507 500 535 513"] [fields="1 2 3 4 5"]      (DFN6, DFN9, DFN15, DFN19)
set -u; source "$(dirname "$0")/env.sh"
SEEDS=${1:-"507 500 535 513"}; FIELDS=${2:-"1 2 3 4 5"}
for K in $FIELDS; do for SEED in $SEEDS; do
  FLOW="$WORK/flow/s$SEED"; RD="$WORK/rxn/random/s$SEED/field$K"; mkdir -p "$RD"
  [ -f "$RD/deck.in" ] || python3 "$HERE/make_decks.py" "$WORK/networks/s$SEED" "$FLOW/flow-vel-001.vtk" "$RD/deck.in" --random-seed $K --sigma 1.0
  throttle; run_reactive "$RD" "$SEED" "$FLOW" &
done; done; wait
