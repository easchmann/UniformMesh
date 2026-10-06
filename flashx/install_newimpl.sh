#!/usr/bin/env bash
# Install the NewImpl Spark sub-unit physics/Hydro/HydroMain/Spark/UniformMeshNewImpl
# into a Flash-X checkout on the RIKEN fork's MOL branch.
#
# This sub-unit drives the standalone NewImpl whole-step driver (hy_prepareAdvance)
# from a Flash-X Hydro() entry point (replacing Spark's Hydro_prepBlock/Hydro_advance
# on the mesh), while keeping um_init / um_copyIn / um_copyOut.
#
# usage: flashx/install_newimpl.sh <flash-x dir> [dim=2d] [halo=4]
#   dim  : 2d or 3d, must match the -2d / -3d of the Flash-X setup
#   halo : the UniformMesh halo count; NewImpl needs MESH_NHALO >= 3 (the driver hard-
#          asserts guardLo <= lo-3 and guardHi >= hi+3). Default 4 (macros/halo_4.ini).
# The layout is always AoS, which matches Flash-X's unk(var,i,j,k).
#
# The NewImpl .F90-mc templates use plain @macro(args)@ syntax, which Flash-X's
# bin/macroProcessor.py cannot expand (it only understands @M macro@). They are
# therefore pre-expanded here with macro_expand.py from ~/newimpl-work and the
# generated .F90 placed in the sub-unit, so the checkout is self-contained.
#
# Sources: ~/UniformMesh (mesh core + this script) and ~/newimpl-work (NewImpl).
# Neither is modified. Does NOT touch the wrapper-path Spark/UniformMesh sub-unit.
set -euo pipefail

FLASHX=${1:?usage: $0 <flash-x dir> [dim=2d] [halo=4]}
DIM=${2:-2d}
HALO=${3:-4}
PYTHON=${PYTHON:-python3}
NEWIMPL=${NEWIMPL:-$HOME/newimpl-work}

REPO=$(cd "$(dirname "$0")/.." && pwd)
SPARK="$FLASHX/source/physics/Hydro/HydroMain/Spark"
DEST="$SPARK/UniformMeshNewImpl"

[ -d "$SPARK" ] || { echo "error: $FLASHX is not a Flash-X checkout (no $SPARK)" >&2; exit 1; }
[ -f "$REPO/macros/dim_$DIM.ini" ] || { echo "error: unknown dim '$DIM' (2d or 3d)" >&2; exit 1; }
HALO_INI="$REPO/macros/halo_$HALO.ini"
[ -f "$HALO_INI" ] || { echo "error: no macros/halo_$HALO.ini (NewImpl needs >= 3)" >&2; exit 1; }
[ -f "$NEWIMPL/macro_expand.py" ] || { echo "error: $NEWIMPL is not newimpl-work (no macro_expand.py)" >&2; exit 1; }

mkdir -p "$DEST"

# 1. generated header: general + dim + layout (AoS) + halo
"$PYTHON" "$REPO/ini_to_h.py" "$REPO/macros/general.ini" "$REPO/macros/dim_$DIM.ini" \
    "$REPO/macros/layout_aos.ini" "$HALO_INI" -o "$DEST/mesh_config.h"

# 2. copy the mesh core
cp "$REPO/mesh.c" "$REPO/mesh.h" "$REPO/mesh_bind.c" "$REPO/mesh_f.F90" "$DEST/"

# 3. the NewImpl macro library + templates (plain @macro@ dialect) live in newimpl-work.
#    Copy the .ini macros into the sub-unit so the directory is self-contained.
#    Expansion order matters: the base storage macros (hydro_helpers etc.) must be
#    loaded before hydro_layout_uniformmesh.ini, which overrides cell_ref/face_ref and
#    the allocate/release helpers for the mesh-backed storage.
# Stage the NewImpl .ini macro files in a temp dir (NOT in $DEST): Flash-X's setup
# scans every *.ini in a unit directory as its own macro defs ([section] headers),
# so the NewImpl macro library must not be installed into the sub-unit. Only the
# pre-expanded .F90 sources go into $DEST.
STAGE=$(mktemp -d)
trap 'rm -rf "$STAGE"' EXIT
for ini in hydro_helpers update_helpers diagnostic_helpers driver_helpers \
           hydro_layout_uniformmesh; do
    cp "$NEWIMPL/$ini.ini" "$STAGE/$ini.ini"
done

# NOTE: macro_expand.py's -d is nargs='+' with default argparse action (store), so
# repeated -d flags OVERWRITE instead of accumulate. Pass every definition file in a
# single -d group:  -d <f1> <f2> ... <fN>
DEFS="$STAGE/hydro_helpers.ini $STAGE/update_helpers.ini $STAGE/diagnostic_helpers.ini \
      $STAGE/driver_helpers.ini $STAGE/hydro_layout_uniformmesh.ini"
DEFS_ARGS="-d $DEFS"

expand_newimpl() {
    # macro_expand.py: -d <defs...> -i <input .F90-mc> -o <output .F90>
    "$PYTHON" "$NEWIMPL/macro_expand.py" $DEFS_ARGS -i "$NEWIMPL/$1.F90-mc" -o "$DEST/$2"
}

for src in hy_getFaceFlux hy_updateSolution hy_shockDetect hy_computeDt; do
    expand_newimpl "$src" "$src.F90"
done
# uniformmesh_grid.F90-mc carries its own @macro_imports/@macro_declarations
expand_newimpl uniformmesh_grid uniformmesh_grid.F90
cp "$NEWIMPL/hydro_grid_contract.F90" "$DEST/hydro_grid_contract.F90"
# the whole-step driver (generated from pre-expanded storage allocs in Hydro's body)
expand_newimpl hy_prepareAdvance hy_prepareAdvance.F90

# 4. the sub-unit source files (Config, Makefile, Hydro.F90-mc driver, um_flashx).
#    Hydro.F90-mc carries NOVARIANTS and no @M macros; Flash-X's setup copies/macro-
#    processes it to the top-level Hydro.F90 in the object tree (replacing Spark's).
cp "$REPO/flashx/UniformMeshNewImpl/Config"     "$DEST/Config"
cp "$REPO/flashx/UniformMeshNewImpl/Makefile"   "$DEST/Makefile"
cp "$REPO/flashx/UniformMeshNewImpl/um_flashx.F90" "$DEST/um_flashx.F90"
cp "$REPO/flashx/UniformMeshNewImpl/Hydro.F90-mc" "$DEST/Hydro.F90-mc"

echo "installed UniformMeshNewImpl into $DEST (dim=$DIM, halo=$HALO, layout=aos)"
