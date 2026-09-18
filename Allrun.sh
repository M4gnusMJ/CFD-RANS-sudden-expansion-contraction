#!/bin/sh
# Group 22 - run all three mesh levels from the same base_case.
#
#   ./Allrun.sh            run coarse, medium and fine
#   ./Allrun.sh medium     run one level only
#
# Each level runs in its own copy (run_<level>/) so the three results coexist,
# which is what the GCI study needs. CSVs land in postprocess/data/.

cd "$(dirname "$0")"

LEVELS=${*:-"coarse medium fine"}

die() {
    echo "  FAILED: $1"
    echo "  ---- last 25 lines of $2 ----"
    tail -25 "$2" | sed 's/^/  /'
    exit 1
}

for LEVEL in $LEVELS; do
    echo "=================== $LEVEL ==================="
    RUN="run_$LEVEL"

    rm -rf "$RUN"
    cp -r base_case "$RUN" || exit 1
    cp "case_design/blockMeshDict.$LEVEL" "$RUN/system/blockMeshDict" || exit 1
    rm -rf "$RUN/constant/polyMesh"

    echo "  blockMesh..."
    blockMesh -case "$RUN" > "$RUN/log.blockMesh" 2>&1 \
        || die blockMesh "$RUN/log.blockMesh"

    echo "  checkMesh..."
    checkMesh -case "$RUN" > "$RUN/log.checkMesh" 2>&1 \
        || die checkMesh "$RUN/log.checkMesh"
    grep -E "cells:|Max aspect ratio|non-orthogonality|Max skewness|Mesh OK|\*\*\*" \
        "$RUN/log.checkMesh" | sed 's/^/    /'

    echo "  simpleFoam..."
    simpleFoam -case "$RUN" > "$RUN/log.simpleFoam" 2>&1 \
        || die simpleFoam "$RUN/log.simpleFoam"
    grep -E "SIMPLE solution converged|^End" "$RUN/log.simpleFoam" | sed 's/^/    /' \
        || echo "    WARNING: no convergence message - check $RUN/log.simpleFoam"

    mkdir -p postprocess/data

    # --- pressure along the axis -------------------------------------------
    # OpenFOAM names the file after every sampled field, sorted alphabetically
    # (line_epsilon_k_p_U.csv), so glob it and find the p column by its header
    # rather than assuming a position.
    AX=$(ls "$RUN"/postProcessing/axisLine/*/line_*.csv 2>/dev/null | tail -1)
    if [ -n "$AX" ]; then
        awk -F',' 'NR==1 {for (i=1;i<=NF;i++) if ($i=="p") c=i; next}
                   c    {print $1","$c}' "$AX" \
            > "postprocess/data/${LEVEL}_axis.csv"
        echo "    axis: $(wc -l < "postprocess/data/${LEVEL}_axis.csv") rows from $(basename "$AX")"
    else
        echo "    NO axisLine output - is foam_functions included in controlDict?"
    fi

    # --- y+ along the wall --------------------------------------------------
    # Face values from the patch surface: x y z yPlus, space separated,
    # comment lines start with #.
    YP=$(ls "$RUN"/postProcessing/yPlusSurface/*/yPlus_wall.raw 2>/dev/null | tail -1)
    if [ -n "$YP" ]; then
        awk '!/^#/ && NF>=4 {print $1","$4}' "$YP" \
            > "postprocess/data/${LEVEL}_yplus.csv"
        echo "    y+  : $(wc -l < "postprocess/data/${LEVEL}_yplus.csv") rows from $(basename "$YP")"
    else
        echo "    NO yPlusSurface output"
    fi

    # --- y+ summary straight from the solver, as a cross-check --------------
    YD=$(ls "$RUN"/postProcessing/yPlusWall/*/yPlus.dat 2>/dev/null | tail -1)
    [ -n "$YD" ] && tail -1 "$YD" | awk '{printf "    y+  : min %.3f  max %.2f  avg %.3f (yPlus.dat)\n", $3, $4, $5}'
done

echo
echo "Done. Outside the container:  python3 postprocess/plot_pressure.py"
