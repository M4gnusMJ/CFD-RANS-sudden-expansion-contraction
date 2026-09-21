# Current state

M1–M3 complete: revised mesh family targets nominal developed-wall y+ 1 / 0.667 / 0.444 with 34,720 / 78,120 / 175,770 cells. User superseded the former coarser family. Generator adds two axial planes 20h from the steps, relaxes internal radial spacing, retains outer-wall spacing, and matches sizes across interfaces and grading joins. README explains helpers and rerun requirements.

Actual isolated meshes pass standard checkMesh and renumberMesh -overwrite. Max aspect ~333, nonorthogonality <0.81 degrees, skewness .332. Coarse max radial growth15.8%, axial2.8%. Optional allGeometry determinant check still fails, also reproduced on previous coarse; see VERIFICATION.md. Production solver runs, postprocess/data and plots were not modified this turn. Improved minor-loss predictions not established.

Relevant code: case_design/case_design.py plus generated dictionaries, case_params.json and design_summary.txt. Prior uncommitted Allpostprocess/baseline plotting changes retained. OpenFOAM v2412 environment /usr/lib/openfoam/openfoam2412/etc/bashrc. Mesh verification artifacts /tmp/physics-mesh-46zuikzd, spacing script /tmp/check_mesh_spacing.py, comparison plot inspected.

Next user action: run all three levels with ./Allrun.sh before interpreting new GCI, since metadata now describes the revised family. Eight-process run, convergence checks, mapping, renumber -overwrite, full postprocessing suite retained. No outstanding implementation work.
