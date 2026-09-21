#!/bin/sh
# Rebuild every analysis and plot from the existing postprocess/data exports.
cd "$(dirname "$0")" || exit 1
mkdir -p postprocess/logs || exit 1
STATUS=0
for SCRIPT in plot_pressure minor_losses velocity_profiles gci; do
    echo "  $SCRIPT..."
    if python3 "postprocess/$SCRIPT.py" > "postprocess/logs/$SCRIPT.log" 2>&1; then
        echo "    OK (postprocess/logs/$SCRIPT.log)"
    else
        echo "    FAILED (postprocess/logs/$SCRIPT.log)"
        tail -n 25 "postprocess/logs/$SCRIPT.log"
        STATUS=1
    fi
done
echo "Postprocessing finished. Figures: postprocess/figs/; tables: postprocess/tables/"
echo "Combined pressure/y+ plot: postprocess/pressure_comparison.png"
exit "$STATUS"
