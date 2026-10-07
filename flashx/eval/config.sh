# Shared settings for flashx/eval/{prepare,build,run}.sh. Every value can be
# overridden from the environment, e.g.  FLASHX=~/other/Flash-X ./build.sh
# Sourced, not executed.

EVAL_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
UM_REPO=$(cd "$EVAL_DIR/../.." && pwd)

# Flash-X checkout of the RIKEN fork, branch MOL, at the commit the evaluation was made with
FLASHX=${FLASHX:-$HOME/Flash-X-MOL}
FLASHX_COMMIT=${FLASHX_COMMIT:-44188e65bc44fac1982b803aa38cc6d5acfa244e}
SITE=${SITE:-supercomp}

# toolchain: system Open MPI + GNU compilers
MPI_PATH=${MPI_PATH:-/usr/lib64/openmpi}
export PATH="$MPI_PATH/bin:$PATH"
JOBS=${JOBS:-16}

# HDF5 1.14.5: parallel build for Flash-X's IO (UG + nofbs selects hdf5/parallel/NoFbs),
# serial build for sfocu
HDF5_VERSION=${HDF5_VERSION:-1.14.5}
HDF5_TARBALL=${HDF5_TARBALL:-$HOME/hdf5-$HDF5_VERSION.tar.gz}
HDF5_URL=${HDF5_URL:-https://support.hdfgroup.org/releases/hdf5/v1_14/v1_14_5/downloads/hdf5-$HDF5_VERSION.tar.gz}
HDF5_PAR=${HDF5_PAR:-$HOME/opt/hdf5-par}
HDF5_SER=${HDF5_SER:-$HOME/opt/hdf5}

# mesh build of the sub-unit (install_newimpl.sh): dimension and halo width
DIM=${DIM:-2d}
HALO=${HALO:-4}

SFOCU="$FLASHX/tools/sfocu/sfocu"

# problems: name, Flash-X simulation, parameter file in this directory
PROBLEMS=("vortex IsentropicVortex vortex.par" "sod Sod sod.par")
# block splits (nbx x nby) for the split-invariance test; 1x1 is the baseline
SPLITS_vortex=(1x1 2x2 4x4 8x8 16x16 4x1 1x8)
SPLITS_sod=(1x1 4x1 8x2 16x4 32x1)
# split used for the NewImpl-vs-Spark comparison
COMPARE_vortex=4x4
COMPARE_sod=4x1

# object directories in $FLASHX: eval_<problem>_ni (NewImpl on the UniformMesh), eval_<problem>_ref (plain Spark)
objdir() { echo "eval_$1_$2"; }

log() { echo "[$(date +%H:%M:%S)] $*"; }
die() { echo "error: $*" >&2; exit 1; }
