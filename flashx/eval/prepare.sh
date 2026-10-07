#!/usr/bin/env bash
# One-time environment for the UniformMeshNewImpl evaluation (idempotent: finished
# steps are detected and skipped).
#   1. check the Flash-X checkout (RIKEN fork, branch MOL, pinned commit)
#   2. build HDF5 $HDF5_VERSION: parallel into $HDF5_PAR (Flash-X IO), serial into $HDF5_SER (sfocu)
#   3. install the site file sites/$SITE/Makefile.h (MPI_PATH, HDF5_PATH substituted)
#   4. remove the two objects MOL's Spark/Makefile lists without sources (plain Spark does not build otherwise)
#   5. build tools/sfocu against the serial HDF5
# usage: flashx/eval/prepare.sh        (settings: flashx/eval/config.sh)
set -euo pipefail
source "$(dirname "$0")/config.sh"

# 1. Flash-X checkout
[ -d "$FLASHX/source/physics/Hydro/HydroMain/Spark/NewImpl" ] || \
    die "$FLASHX is not a checkout of the RIKEN fork's MOL branch (no Spark/NewImpl)"
head=$(git -C "$FLASHX" rev-parse HEAD)
if [ "$head" != "$FLASHX_COMMIT" ]; then
    [ "${ALLOW_OTHER_COMMIT:-0}" = 1 ] || \
        die "$FLASHX is at $head, expected $FLASHX_COMMIT (set ALLOW_OTHER_COMMIT=1 to continue anyway)"
    log "warning: $FLASHX is at $head, not the pinned $FLASHX_COMMIT"
fi
log "Flash-X: $FLASHX @ $head"

# 2. HDF5
build_hdf5() {  # $1 = prefix, $2 = parallel (yes/no)
    local prefix=$1 parallel=$2 settings="$1/lib/libhdf5.settings"
    if [ -f "$settings" ] && grep -q "HDF5 Version: $HDF5_VERSION" "$settings" && \
       grep -q "Parallel HDF5: $parallel" "$settings"; then
        log "HDF5 (parallel=$parallel) already in $prefix"
        return
    fi
    [ -f "$HDF5_TARBALL" ] || { log "downloading $HDF5_URL"; curl -fL -o "$HDF5_TARBALL" "$HDF5_URL"; }
    local src; src=$(mktemp -d "${TMPDIR:-/tmp}/hdf5-build.XXXXXX")
    tar xzf "$HDF5_TARBALL" -C "$src"
    # the release tarball's entries start with ./hdf5-<version>/, so locate configure
    local top; top=$(dirname "$(find "$src" -maxdepth 3 -name configure -type f | head -1)")
    local cc=gcc extra=""
    if [ "$parallel" = yes ]; then cc=mpicc; extra=--enable-parallel; fi
    log "building HDF5 $HDF5_VERSION (parallel=$parallel) into $prefix (log: $src/build.log)"
    (cd "$top" && CC=$cc ./configure --prefix="$prefix" $extra --disable-fortran --disable-cxx \
         --enable-build-mode=production --disable-tests --disable-tools-tests && \
     make -j "$JOBS" && make install) > "$src/build.log" 2>&1 || die "HDF5 build failed, see $src/build.log"
    grep -q "Parallel HDF5: $parallel" "$settings" || die "HDF5 in $prefix is not parallel=$parallel"
    rm -rf "$src"
}
build_hdf5 "$HDF5_PAR" yes
build_hdf5 "$HDF5_SER" no

# 3. site file
sitefile="$FLASHX/sites/$SITE/Makefile.h"
mkdir -p "$(dirname "$sitefile")"
# paths under $HOME are written as $(HOME)/..., as in the committed site file, so an
# unchanged setting reproduces it exactly (and a rerun neither changes nor backs it up)
make_path() { case "$1" in "$HOME"/*) echo "\$(HOME)${1#"$HOME"}" ;; *) echo "$1" ;; esac; }
new=$(mktemp)
sed -e "s|^MPI_PATH .*|MPI_PATH   = $(make_path "$MPI_PATH")|" \
    -e "s|^HDF5_PATH .*|HDF5_PATH  = $(make_path "$HDF5_PAR")|" \
    "$EVAL_DIR/site/Makefile.h" > "$new"
if [ -f "$sitefile" ] && ! cmp -s "$new" "$sitefile"; then
    cp "$sitefile" "$sitefile.bak.$(date +%Y%m%d%H%M%S)"
    log "backed up the previous $sitefile"
fi
mv "$new" "$sitefile"
log "site file: $sitefile (HDF5_PATH = $HDF5_PAR)"

# 4. MOL's Spark/Makefile lists hy_getFaceFlux.o and hy_updateSolution.o, whose sources
#    exist only as NewImpl templates; remove them so plain Spark builds (the
#    UniformMeshNewImpl sub-unit lists both objects in its own Makefile).
spark_mk="$FLASHX/source/physics/Hydro/HydroMain/Spark/Makefile"
if grep -qE 'hy_getFaceFlux\.o|hy_updateSolution\.o' "$spark_mk"; then
    sed -i.orig -e 's/ hy_getFaceFlux\.o//' -e 's/ hy_updateSolution\.o//' "$spark_mk"
    log "removed hy_getFaceFlux.o/hy_updateSolution.o from $spark_mk (original: $spark_mk.orig)"
else
    log "Spark/Makefile already builds plain Spark"
fi

# 5. sfocu (serial HDF5, no MPI, no PnetCDF)
if [ ! -x "$SFOCU" ] || [ "${REBUILD_SFOCU:-0}" = 1 ]; then
    log "building sfocu"
    (cd "$FLASHX/tools/sfocu" && make clean > /dev/null 2>&1 || true
     make SITE="$SITE" NO_MPI=True NO_NCDF=True CCOMP=gcc \
          CFLAGS_HDF5="-I$HDF5_SER/include -DH5_USE_18_API" \
          LIB_HDF5="-L$HDF5_SER/lib -lhdf5 -lz -ldl -lm -Wl,-rpath,$HDF5_SER/lib") \
        > "$FLASHX/tools/sfocu/build.log" 2>&1 || die "sfocu build failed, see $FLASHX/tools/sfocu/build.log"
fi
log "sfocu: $SFOCU"
log "prepare done"
