#!/bin/bash
# Steady flow field of one network:  bash run_flow.sh <seed> [Q_factor=1]
# Output: $WORK/flow/s<seed>[_Q<factor>]/flow-restart.h5 (checkpoint) and flow-vel-001.vtk (cell-centred Darcy velocity)
set -eu; source "$(dirname "$0")/env.sh"
SEED=$1; QF=${2:-1}
FD="$WORK/flow/s$SEED"; [ "$QF" != 1 ] && FD="${FD}_Q$QF"
mkdir -p "$FD"; cd "$FD"
for f in $MESH_FILES; do ln -sf "$WORK/networks/s$SEED/$f" "$f"; done
RATE=$(python3 -c "print(('%.4e' % (5e-7*$QF)).replace('e','d'))")            # baseline Q = 5e-7 m3/d
sed -e "s/RATE 5\.d-7 m^3\/d/RATE $RATE m^3\/d/" -e "s/RATE -5\.d-7 m^3\/d/RATE -$RATE m^3\/d/" "$INPUTS/flow.in" > flow.in
mpirun -np 1 "$PFLOTRAN" -pflotranin flow.in > flow.log 2>&1
ls -la flow-restart.h5 flow-vel-001.vtk
