#!/bin/sh
# Group 22 - run all three mesh levels from the same base_case.
#
#   ./Allrun.sh            run coarse, medium and fine
#   ./Allrun.sh medium     run one level only
#
# Each level runs in its own copy (run_<level>/) so the three results coexist,
# which is what the GCI study needs. CSVs land in postprocess/data/.
# Medium/fine start from the latest saved coarse/medium solution when available
# and geometrically compatible; otherwise they use the base-case initial fields.

cd "$(dirname "$0")"

LEVELS=${*:-"coarse medium fine"}
python3 postprocess/sample_wall_profiles.py --generate || exit 1

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

    echo "  renumberMesh..."
    renumberMesh -case "$RUN" -overwrite > "$RUN/log.renumberMesh" 2>&1 \
        || die renumberMesh "$RUN/log.renumberMesh"

    echo "  checkMesh..."
    checkMesh -case "$RUN" > "$RUN/log.checkMesh" 2>&1 \
        || die checkMesh "$RUN/log.checkMesh"
    grep -E "cells:|Max aspect ratio|non-orthogonality|Max skewness|Mesh OK|\*\*\*" \
        "$RUN/log.checkMesh" | sed 's/^/    /'

    SOURCE=""
    case "$LEVEL" in
        medium) SOURCE=run_coarse ;;
        fine) SOURCE=run_medium ;;
    esac
    if [ -n "$SOURCE" ] && [ -d "$SOURCE/constant/polyMesh" ]; then
        SOURCE_TIME=$(foamListTimes -case "$SOURCE" -latestTime -noZero)
        # Avoid mapping an old run made before a geometry change.
        SOURCE_VERTICES=$(foamDictionary "$SOURCE/system/blockMeshDict" -entry vertices -value) || exit 1
        TARGET_VERTICES=$(foamDictionary "$RUN/system/blockMeshDict" -entry vertices -value) || exit 1
        if [ -n "$SOURCE_TIME" ] && [ "$SOURCE_VERTICES" = "$TARGET_VERTICES" ]; then
            echo "  mapFields from $SOURCE (time $SOURCE_TIME)..."
            mapFields "$SOURCE" -case "$RUN" -sourceTime "$SOURCE_TIME" \
                -consistent -mapMethod interpolate > "$RUN/log.mapFields" 2>&1 \
                || die mapFields "$RUN/log.mapFields"
        else
            echo "  Using base initial fields: $SOURCE has no saved solution or different geometry."
        fi
    elif [ -n "$SOURCE" ]; then
        echo "  Using base initial fields: $SOURCE is unavailable."
    fi

    echo "  decomposePar (8 processes)..."
    decomposePar -case "$RUN" > "$RUN/log.decomposePar" 2>&1 \
        || die decomposePar "$RUN/log.decomposePar"

    echo "  simpleFoam (8 processes)..."
    mpirun -np 8 simpleFoam -case "$RUN" -parallel > "$RUN/log.simpleFoam" 2>&1 \
        || die simpleFoam "$RUN/log.simpleFoam"
    if grep -q 'SIMPLE solution converged' "$RUN/log.simpleFoam"; then
        grep 'SIMPLE solution converged' "$RUN/log.simpleFoam" | sed 's/^/    /'
    elif grep -q 'equationInitialResidual: convergedUz condition satisfied' "$RUN/log.simpleFoam"; then
        echo "    Component-wise convergence criteria satisfied."
    else
        die "convergence criteria not satisfied (iteration cap or interrupted run)" "$RUN/log.simpleFoam"
    fi

    echo "  reconstructPar (latest time)..."
    reconstructPar -case "$RUN" -latestTime > "$RUN/log.reconstructPar" 2>&1 \
        || die reconstructPar "$RUN/log.reconstructPar"

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

    # --- upstream velocity and wall diagnostics ----------------------------
    python3 postprocess/sample_wall_profiles.py "$LEVEL" || exit 1

    # --- y+ summary straight from the solver, as a cross-check --------------
    YD=$(ls "$RUN"/postProcessing/yPlusWall/*/yPlus.dat 2>/dev/null | tail -1)
    [ -n "$YD" ] && tail -1 "$YD" | awk '{printf "    y+  : min %.3f  max %.2f  avg %.3f (yPlus.dat)\n", $3, $4, $5}'
done

echo
echo "=================== postprocessing ==================="
mkdir -p postprocess/logs
for SCRIPT in plot_pressure minor_losses velocity_profiles gci; do
    echo "  $SCRIPT..."
    python3 "postprocess/$SCRIPT.py" > "postprocess/logs/$SCRIPT.log" 2>&1 \
        || die "$SCRIPT" "postprocess/logs/$SCRIPT.log"
done
echo "Done. Postprocessing logs: postprocess/logs/"
