# Upstream velocity diagnostics

User requests velocity profiles 3 upstream diameters before expansion and contraction to investigate developed flow and turbulence-model wall behaviour. Keep numerical inputs in postprocess/data; output figures and tables in existing directories.

Acceptance: extract saved fields without solving again; sample target and adjacent 4D/2D stations; use measured wall shear for inner scaling; report mesh comparisons, streamwise profile changes, wall y+, and limitations without asserting experimental validation; integrate export into future runs and document reproducible commands.

Assumption: diameter means the local upstream diameter (x=2.7 and 8.4 m for current geometry). Adjacent stations are diagnostic comparisons.

## Run and mesh revision (2026-09-21)
Explain mapped time vs convergence; reduce high aspect ratios based on cell locations; current medium resolution becomes new fine, retain systematic ~1.5 refinement; run on eight MPI processes, preserve renumber -overwrite and mapping. Validate actual meshes and short parallel runs with reconstruction/export in temporary cases. Preserve production solutions.

Confirmed: old medium becomes new fine. Added request: raise iteration cap and accept only convergence-based termination; permit larger normalized Uz residual while retaining Ux/Uy thresholds.

## M: physics-focused mesh revision (2026-09-21)
User supersedes R2 resolution targets: coarse design y+ about 1, successive directional refinement about 1.5. Minimal effective generator change: remove persistent wall-scale clustering at internal old-radius interface, retain lip and outer-wall resolution, broad downstream refinement and smooth spacing. Preserve geometry, solver and production data. Verify all three actual meshes and spacing/refinement; improved loss predictions require new converged simulations and are not a mesh-only acceptance claim.
