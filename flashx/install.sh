#!/usr/bin/env bash
# Install the UniformMesh into a Flash-X checkout as the Spark sub-unit
# physics/Hydro/HydroMain/Spark/UniformMesh (see flashx/README.md).
#
# usage: flashx/install.sh <flash-x dir> [dim=2d] [halo=8]
#   dim  : 2d or 3d, must match the -2d / -3d of the Flash-X setup
#   halo : must equal Spark's GUARDCELLS (8 for RK2 with the default WENO/MP5, 4 for TVD/lim03/fog)
# The layout is always AoS, which matches Flash-X's unk(var,i,j,k).
set -euo pipefail

FLASHX=${1:?usage: $0 <flash-x dir> [dim=2d] [halo=8]}
DIM=${2:-2d}
HALO=${3:-8}
PYTHON=${PYTHON:-python3}

REPO=$(cd "$(dirname "$0")/.." && pwd)
SPARK="$FLASHX/source/physics/Hydro/HydroMain/Spark"
DEST="$SPARK/UniformMesh"

[ -d "$SPARK" ] || { echo "error: $FLASHX is not a Flash-X checkout (no $SPARK)" >&2; exit 1; }
[ -f "$REPO/macros/dim_$DIM.ini" ] || { echo "error: unknown dim '$DIM' (2d or 3d)" >&2; exit 1; }
HALO_INI=""
if [ "$HALO" != 1 ]; then
    HALO_INI="$REPO/macros/halo_$HALO.ini"
    [ -f "$HALO_INI" ] || { echo "error: no macros/halo_$HALO.ini" >&2; exit 1; }
fi

mkdir -p "$DEST"
# only the generated header goes into Flash-X: its setup treats *.ini files as its own macros
"$PYTHON" "$REPO/ini_to_h.py" "$REPO/macros/general.ini" "$REPO/macros/dim_$DIM.ini" \
    "$REPO/macros/layout_aos.ini" $HALO_INI -o "$DEST/mesh_config.h"
cp "$REPO/mesh.c" "$REPO/mesh.h" "$REPO/mesh_bind.c" "$REPO/mesh_f.F90" "$DEST/"
cp "$REPO/flashx/UniformMesh/"* "$DEST/"

echo "installed UniformMesh into $DEST (dim=$DIM, halo=$HALO, layout=aos)"
