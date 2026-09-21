#!/bin/bash
# Series A3: one-at-a-time parameter sensitivity on three networks (DFN6, DFN9, DFN19) and three arrangements.
#   bash run_sensitivity.sh [seeds="507 500 513"]
# The configuration O2 = 5.6d-5 / seed 513 / worst arrangement is numerically intractable over 3,000 yr and is skipped.
set -u; source "$(dirname "$0")/env.sh"
SEEDS=${1:-"507 500 513"}
for P in "kb:1.d-7:1.d-6" "O2:1.4d-5:5.6d-5" "kcal:0.1:10" "kpyr:0.1:10" "KO2:6.d-7:6.d-5"; do
  IFS=: read -r NAME V1 V2 <<< "$P"
  for VAL in $V1 $V2; do for SEED in $SEEDS; do for ARR in "0 0" "1 -2" "-1 2"; do set -- $ARR
    [ "$NAME $VAL $SEED $1" = "O2 5.6d-5 513 -1" ] && continue
    FLOW="$WORK/flow/s$SEED"; RD="$WORK/rxn/sensitivity/s$SEED/${NAME}_$VAL/$(label $1 $2)"; mkdir -p "$RD"
    if [ ! -f "$RD/deck.in" ]; then
      python3 "$HERE/make_decks.py" "$WORK/networks/s$SEED" "$FLOW/flow-vel-001.vtk" "$RD/deck.in" --beta-pyr $1 --beta-cal $2
      python3 "$HERE/set_parameter.py" "$RD/deck.in" "$NAME" "$VAL"
    fi
    throttle; run_reactive "$RD" "$SEED" "$FLOW" &
  done; done; done
done; wait
