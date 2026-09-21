#!/bin/bash
# Flow-correlated arrangement factorial (Section 2.3): 3 pyrite x 5 calcite arrangements per network, 3,000 yr.
#   bash run_factorial.sh [seeds="all"]        e.g.  bash run_factorial.sh "500 513"
# Requires: build_network.py and run_flow.sh done for each seed.
set -u; source "$(dirname "$0")/env.sh"
SEEDS=${1:-all}; [ "$SEEDS" = all ] && SEEDS=$(all_seeds)
for SEED in $SEEDS; do
  FLOW="$WORK/flow/s$SEED"
  for BC in 2 1 0 -1 -2; do for BP in 1 0 -1; do          # calcite outer, pyrite inner: spreads the stiff pyrite anti-aligned runs
    RD="$WORK/rxn/factorial/s$SEED/$(label $BP $BC)"; mkdir -p "$RD"
    [ -f "$RD/deck.in" ] || python3 "$HERE/make_decks.py" "$WORK/networks/s$SEED" "$FLOW/flow-vel-001.vtk" "$RD/deck.in" --beta-pyr $BP --beta-cal $BC
    throttle; run_reactive "$RD" "$SEED" "$FLOW" &
  done; done
done; wait
