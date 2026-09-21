# Common settings, sourced by every run script. Override from the environment if needed.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export INPUTS="$HERE/../1_inputs"
export WORK="${WORK:-$HERE/../work}"                 # all simulation output goes here
export PFLOTRAN="${PFLOTRAN:-/dfnWorks/bin/pflotran}"  # path inside the dfnWorks v2.9 Docker image
export NP="${NP:-4}"                                 # MPI ranks per reactive transport run
export NPAR="${NPAR:-3}"                             # runs executed concurrently
[ -d /dfnWorks/lib/petsc/arch-linux2-c-debug/bin ] && export PATH=/dfnWorks/lib/petsc/arch-linux2-c-debug/bin:$PATH  # mpirun
MESH_FILES="full_mesh_vol_area.uge full_mesh.vtk dfn_properties.h5 boundary_left_w.ex boundary_right_e.ex boundary_top.ex boundary_bottom.ex boundary_front_n.ex boundary_back_s.ex pinned.ex"

# seeds of the twenty networks of the paper (column 2 of networks.csv)
all_seeds() { tail -n +2 "$INPUTS/networks.csv" | cut -d, -f2; }
# arrangement label used in directory names, e.g. "pyr+1_cal-2"
label() { printf "pyr%+d_cal%+d" "$1" "$2"; }
# wait until fewer than NPAR background jobs are running
throttle() { while [ "$(jobs -rp | wc -l)" -ge "$NPAR" ]; do sleep 5; done; }

# run_reactive <run_dir> <seed> <flow_dir>: link mesh, flow restart and database next to deck.in, then run PFLOTRAN
run_reactive() {
  local rd=$1 seed=$2 flow=$3
  if grep -q "^rc=0" "$rd/run.log" 2>/dev/null; then echo "  skip (done): $rd"; return; fi
  ( cd "$rd" || exit 1
    for f in $MESH_FILES; do ln -sf "$WORK/networks/s$seed/$f" "$f"; done
    ln -sf "$flow/flow-restart.h5" flow-restart.h5
    ln -sf "$INPUTS/hanford_neaUOH4.dat" hanford_neaUOH4.dat
    mpirun -np "$NP" "$PFLOTRAN" -pflotranin deck.in > run.log 2>&1; echo "rc=$?" >> run.log )
  echo "  finished: $rd ($(tail -1 "$rd/run.log"))"
}
