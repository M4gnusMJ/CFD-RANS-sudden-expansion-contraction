# Current state

Repository: assignment_3_expansion_contraction. Existing pressure GCI work committed by user before this task; initial working tree clean.
OpenFOAM v2412 available after sourcing /usr/lib/openfoam/openfoam2412/etc/bashrc. Saved final fields and meshes present for all levels.
Existing wallProfiles includes x=2.7 but not 8.4 and is not exported. Wedge wall is at R*cos(2 degrees) at z=0, so wall distance must use actual sampled face coordinates.
V1-V3 complete; see VERIFICATION.md. Generated 18 profiles (240 points each), two plots, three analysis tables, per-mesh raw exports and metadata. Allrun generates the include and exports future results automatically.

Use `python3 postprocess/sample_wall_profiles.py --resample` in an OpenFOAM shell to refresh saved-run samples, then `python3 postprocess/velocity_profiles.py` to analyse only data-folder exports.

Sampling correction: cellPoint gave nonzero velocity offsets below the first cell. cellPointFace improved but did not eliminate interpolation artefacts there. Retain all raw samples with below_first_cell flag, but restrict wall-scaled plots/sublayer metrics to y+ at or above the local wall-cell y+. Do not interpret denser interpolation points as independent resolved cells.

Findings: fine-grid mean-profile changes between 3D and 4D/2D below 0.2% area-weighted RMS of nominal bulk speed. Fine wall y+ approximately .44 and .61. Viscous-sublayer RMS departures .31%/.67%; contraction log-region departure and mesh sensitivity remain substantial. No experimental validation claimed. No remaining implementation tasks.
