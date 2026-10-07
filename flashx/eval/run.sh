#!/usr/bin/env bash
# Run the UniformMeshNewImpl evaluation with the executables from build.sh:
#   A. block-split invariance: every split of a problem must give a bitwise identical
#      final checkpoint to the 1x1 run (sfocu SUCCESS)                      -> PASS/FAIL
#   B. NewImpl vs plain Spark, same parameters:
#      vortex (smooth): max sfocu mag error <= VORTEX_TOL                   -> PASS/FAIL
#      sod (shocks): the schemes differ near shocks by design               -> INFO only
#      (error norms, and the sum of velx, which is 0 for a solution that keeps the
#      problem's mirror symmetry about x = 0.25)
#   C. both kinds take the same number of steps                            -> PASS/FAIL
# Results: $OUT/summary.txt (+ versions.txt, every run directory and sfocu report).
# Exit status 1 if any check fails.
# usage: flashx/eval/run.sh        (settings: flashx/eval/config.sh, OUT=<dir> to choose the output)
set -euo pipefail
source "$(dirname "$0")/config.sh"

VORTEX_TOL=${VORTEX_TOL:-1e-12}
OUT=${OUT:-$FLASHX/eval_runs/$(date +%Y%m%d-%H%M%S)}
mkdir -p "$OUT"
SUMMARY="$OUT/summary.txt"
: > "$SUMMARY"
failures=0

[ -x "$SFOCU" ] || die "no sfocu, run prepare.sh first"

report() {  # $1 = PASS|FAIL|INFO, rest = text
    local status=$1; shift
    printf '%-4s  %s\n' "$status" "$*" | tee -a "$SUMMARY"
    [ "$status" != FAIL ] || failures=$((failures + 1))
}

# what was run: commits (and uncommitted changes), libraries, compilers
{
    echo "date:      $(date -Is)"
    echo "host:      $(hostname)"
    echo "Flash-X:   $FLASHX @ $(git -C "$FLASHX" rev-parse HEAD)"
    git -C "$FLASHX" status --short --untracked-files=no | sed 's/^/           /'
    echo "UniformMesh: $UM_REPO @ $(git -C "$UM_REPO" rev-parse HEAD)"
    git -C "$UM_REPO" status --short | sed 's/^/           /'
    echo "HDF5 (IO): $(grep 'HDF5 Version' "$HDF5_PAR/lib/libhdf5.settings" | sed 's/.*: //'), $HDF5_PAR"
    echo "mpif90:    $(mpif90 --version | head -1)"
    echo "mesh:      dim=$DIM halo=$HALO"
} > "$OUT/versions.txt"

lastchk() { ls "$1"/*_hdf5_chk_* 2>/dev/null | sort | tail -1; }
nsteps() { grep -cE '^ +[0-9]+ [0-9]\.[0-9]+E[-+][0-9]+ ' "$1/run.out" || true; }
# max over variables of one column of sfocu's second table (2 = L1-ErrNorm, 3 = Mag Error)
sfocu_max() {
    awk -F'|' -v c="$2" '/^Var *[|] *L1-ErrNorm/ {t = 1; next}
        t && $1 ~ /^[a-z]+ *$/ {v = $c + 0; if (v > m) m = v} END {printf "%.3e", m}' "$1"
}
sfocu_col() {  # $1 = report, $2 = variable, $3 = column of the second table
    awk -F'|' -v var="$2" -v c="$3" '/^Var *[|] *L1-ErrNorm/ {t = 1; next}
        t {n = $1; gsub(/ /, "", n); if (n == var) {gsub(/ /, "", $c); print $c; exit}}' "$1"
}

run_case() {  # $1 = executable dir, $2 = run dir, $3 = par file, $4 = nbx, $5 = nby (empty: plain Spark)
    local exe=$1 dir=$2 par=$3
    mkdir -p "$dir"
    if [ -n "${4:-}" ]; then
        sed -e "s/^um_nblockx .*/um_nblockx = $4/" -e "s/^um_nblocky .*/um_nblocky = $5/" \
            -e "s/^um_dumpInterval .*/um_dumpInterval = 0/" "$par" > "$dir/flash.par"
    else
        grep -v '^um_' "$par" > "$dir/flash.par"   # plain Spark does not know the um_ parameters
    fi
    (cd "$dir" && mpirun -np 1 "$exe/flashx" > run.out 2>&1 < /dev/null) || {
        report FAIL "$(basename "$dir"): flashx exited with an error, see $dir/run.out"; return; }
    [ -n "$(lastchk "$dir")" ] || report FAIL "$(basename "$dir"): no checkpoint written"
}

for p in "${PROBLEMS[@]}"; do
    read -r name sim par <<< "$p"
    ni="$FLASHX/$(objdir "$name" ni)"; ref="$FLASHX/$(objdir "$name" ref)"
    [ -x "$ni/flashx" ] && [ -x "$ref/flashx" ] || die "missing executables for $name, run build.sh first"
    splits_var="SPLITS_$name[@]"; compare_var="COMPARE_$name"
    splits=("${!splits_var}"); compare=${!compare_var}

    log "$name: running ${#splits[@]} block splits and the Spark reference"
    for s in "${splits[@]}"; do
        run_case "$ni" "$OUT/$name/ni_$s" "$EVAL_DIR/$par" "${s%x*}" "${s#*x}"
    done
    run_case "$ref" "$OUT/$name/ref" "$EVAL_DIR/$par"

    echo "== $name ($sim): block-split invariance (vs 1x1)" | tee -a "$SUMMARY"
    base=$(lastchk "$OUT/$name/ni_1x1")
    for s in "${splits[@]:1}"; do
        rep="$OUT/$name/sfocu_ni_${s}_vs_1x1.txt"
        "$SFOCU" "$base" "$(lastchk "$OUT/$name/ni_$s")" > "$rep" 2>&1 || true
        if [ "$(tail -1 "$rep")" = SUCCESS ]; then
            report PASS "$name split $s: bitwise identical to 1x1"
        else
            report FAIL "$name split $s: differs from 1x1 (max mag error $(sfocu_max "$rep" 3)), see $rep"
        fi
    done

    echo "== $name ($sim): NewImpl ($compare) vs plain Spark" | tee -a "$SUMMARY"
    rep="$OUT/$name/sfocu_ni_${compare}_vs_ref.txt"
    "$SFOCU" "$(lastchk "$OUT/$name/ni_$compare")" "$(lastchk "$OUT/$name/ref")" > "$rep" 2>&1 || true
    mag=$(sfocu_max "$rep" 3); l1=$(sfocu_max "$rep" 2)
    if [ "$name" = vortex ]; then
        if awk -v m="$mag" -v t="$VORTEX_TOL" 'BEGIN {exit !(m <= t)}'; then
            report PASS "$name: max mag error $mag <= $VORTEX_TOL (max L1 $l1)"
        else
            report FAIL "$name: max mag error $mag > $VORTEX_TOL (max L1 $l1), see $rep"
        fi
    else
        report INFO "$name: max mag error $mag, max L1 $l1, sum(velx) NewImpl $(sfocu_col "$rep" velx 5) / Spark $(sfocu_col "$rep" velx 9), see $rep"
    fi
    a=$(nsteps "$OUT/$name/ni_$compare"); b=$(nsteps "$OUT/$name/ref")
    if [ "$a" = "$b" ] && [ "$a" -gt 0 ]; then
        report PASS "$name: same number of steps ($a)"
    else
        report FAIL "$name: steps NewImpl $a vs Spark $b"
    fi
done

# D. exact-solution errors, symmetry, plots (information only; needs numpy, h5py, matplotlib)
echo "== plots and errors vs exact solution" | tee -a "$SUMMARY"
if python3 -c "import numpy, h5py, matplotlib" 2> /dev/null; then
    python3 "$EVAL_DIR/plot.py" "$OUT" --compare "vortex=$COMPARE_vortex" "sod=$COMPARE_sod" | tee -a "$SUMMARY" ||
        echo "INFO  plot.py reported an error (see above); the checks are unaffected" | tee -a "$SUMMARY"
else
    echo "INFO  skipped: python3 lacks numpy/h5py/matplotlib (python3 -m pip install --user h5py); run flashx/eval/plot.py $OUT later" | tee -a "$SUMMARY"
fi

echo "== $failures failed check(s); results in $OUT" | tee -a "$SUMMARY"
[ "$failures" -eq 0 ]
