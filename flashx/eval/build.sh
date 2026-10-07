#!/usr/bin/env bash
# Build the four executables of the evaluation (run prepare.sh once before):
#   $FLASHX/eval_<problem>_ni   NewImpl on the UniformMesh (--with-unit=.../UniformMeshNewImpl)
#   $FLASHX/eval_<problem>_ref  plain Spark (the reference)
# for <problem> in vortex (IsentropicVortex) and sod (Sod), all 2D, +spark +nofbs, with
# HDF5 IO. The sub-unit is (re)installed from this repo first, so the build always uses
# the working copy's mesh and driver code. Logs: $FLASHX/eval_<problem>_<kind>/{setup,make}.log
# usage: flashx/eval/build.sh        (settings: flashx/eval/config.sh)
set -euo pipefail
source "$(dirname "$0")/config.sh"

[ -f "$FLASHX/sites/$SITE/Makefile.h" ] || die "no site file for $SITE, run prepare.sh first"

log "installing UniformMeshNewImpl (dim=$DIM, halo=$HALO)"
"$UM_REPO/flashx/install_newimpl.sh" "$FLASHX" "$DIM" "$HALO"

build_one() {  # $1 = simulation, $2 = objdir, $3.. = extra setup arguments
    local sim=$1 obj=$2; shift 2
    log "setup + make $obj ($sim)"
    rm -rf "${FLASHX:?}/$obj"
    # ./setup (not ./bin/setup.py): it changes into bin/, where setup_shortcuts.txt is;
    # otherwise every +shortcut is silently ignored
    (cd "$FLASHX" && ./setup "$sim" -auto -"$DIM" +spark +nofbs -site="$SITE" -objdir="$obj" "$@") \
        > "$FLASHX/setup_$obj.log" 2>&1 || die "setup of $obj failed, see $FLASHX/setup_$obj.log"
    mv "$FLASHX/setup_$obj.log" "$FLASHX/$obj/setup.log"
    grep -q "IO/IOMain/hdf5/parallel/NoFbs" "$FLASHX/$obj/Units" || die "$obj: HDF5 IO (parallel/NoFbs) was not selected"
    (cd "$FLASHX/$obj" && make -j "$JOBS" > make.log 2>&1) || die "make of $obj failed, see $FLASHX/$obj/make.log"
    [ -x "$FLASHX/$obj/flashx" ] || die "$obj: no flashx executable"
}

built=" "
for p in "${PROBLEMS[@]}"; do
    read -r _ sim _ build <<< "$p"
    [[ "$built" == *" $build "* ]] && continue   # problems sharing a build differ only in parameters
    build_one "$sim" "$(objdir "$build" ni)" --with-unit=physics/Hydro/HydroMain/Spark/UniformMeshNewImpl
    grep -q "hy_prepareAdvance" "$FLASHX/$(objdir "$build" ni)/Hydro.F90" || \
        die "$(objdir "$build" ni): Hydro.F90 is not the UniformMeshNewImpl driver"
    build_one "$sim" "$(objdir "$build" ref)"
    built+="$build "
done
log "build done"
