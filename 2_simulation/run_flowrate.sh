#!/bin/bash
# Series A2: flow rate x3 and x0.3, each with the inlet O2 concentration unchanged ("O2fixed") or scaled by 1/Q
# so that the O2 loading rate stays at its baseline value ("O2scaled"); homogeneous, best and worst arrangements.
#   bash run_flowrate.sh [seeds="507 500 535 513"]
set -u; source "$(dirname "$0")/env.sh"
SEEDS=${1:-"507 500 535 513"}
for SEED in $SEEDS; do for QF in 3 0.3; do
  [ -f "$WORK/flow/s${SEED}_Q$QF/flow-restart.h5" ] || bash "$HERE/run_flow.sh" "$SEED" "$QF"
  FLOW="$WORK/flow/s${SEED}_Q$QF"
  RATE=$(python3 -c "print(('%.4e' % (5e-7*$QF)).replace('e','d'))")
  for O2MODE in O2fixed O2scaled; do
    O2=$(python3 -c "print(('%.4e' % (2.813e-5/($QF if '$O2MODE'=='O2scaled' else 1))).replace('e','d'))")
    for ARR in "0 0" "1 -2" "-1 2"; do set -- $ARR
      RD="$WORK/rxn/flowrate/s$SEED/Q${QF}_$O2MODE/$(label $1 $2)"; mkdir -p "$RD"
      if [ ! -f "$RD/deck.in" ]; then
        python3 "$HERE/make_decks.py" "$WORK/networks/s$SEED" "$WORK/flow/s$SEED/flow-vel-001.vtk" "$RD/deck.in" --beta-pyr $1 --beta-cal $2   # arrangement is defined on the baseline flow field
        sed -i.bak -e "s/RATE 5\.d-7 m^3\/d/RATE $RATE m^3\/d/" -e "s/RATE -5\.d-7 m^3\/d/RATE -$RATE m^3\/d/" -e "s/O2(aq)    2\.813d-5  T/O2(aq)    $O2  T/" "$RD/deck.in" && rm -f "$RD/deck.in.bak"
      fi
      throttle; run_reactive "$RD" "$SEED" "$FLOW" &
    done
  done
done; done; wait
