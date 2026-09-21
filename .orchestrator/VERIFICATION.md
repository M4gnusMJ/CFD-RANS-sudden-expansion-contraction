# Upstream diagnostics verification

Scope: V1-V3 in the current workspace; new sampling/analysis scripts, generated dictionaries, Allrun and README changes, exported data/figures/tables. Existing pressure GCI behavior retained. OpenFOAM v2412; Python with NumPy, pandas, SciPy and Matplotlib.

| Criterion | Check and observed result | Verdict |
|---|---|---|
| Saved solution extraction | `source /usr/lib/openfoam/openfoam2412/etc/bashrc` then `python3 postprocess/sample_wall_profiles.py --resample`; 1,440 profile points per mesh; times 1861/3480/6274. Logs under run directories; function objects only. | passed |
| Required stations and adjacent comparisons | Six stations per mesh including x=2.7/8.4; 18 profiles and 12 comparisons in final analysis. | passed |
| Wall scaling and development metrics | `python3 -m unittest discover -s postprocess -p 'test_*.py'`: six tests pass, including existing pressure GCI tests, exact linear viscous scaling, parabolic area-weighted change, numeric time selection and invalid wall data. | passed |
| Data-folder analysis | `MPLCONFIGDIR=/tmp/velocity-matplotlib python3 postprocess/velocity_profiles.py`; outputs generated solely from exported CSV/JSON. All 18 bulk-integral deviations below 0.2%; finite summary metrics. | passed |
| Visual output | Both PNGs opened and inspected; readable mesh comparisons and wall-law guides. Final wall-scaled curves omit unsupported below-first-cell interval. | passed |
| Sampling integration | `foamDictionary base_case/system/controlDict -entry functions -value` parses generated include; `sh -n Allrun.sh` passes. Actual diagnostic dictionary run with OpenFOAM on all three meshes. | passed |
| Missing data | Run analysis from /tmp against an empty temporary data directory: exit 1 with extraction instructions. | passed |
| Patch hygiene | `git diff --check` passes. | passed |

No new time-marching solve was run; future Allrun solver execution not exercised. Sampling functions were executed on all existing final fields, and base-case dictionary integration was parsed. Wall-law consistency does not establish experimental accuracy or equilibrium turbulence; contraction-side disagreement remains visible.

## 2026-09-21 solver startup repair

Scope: duplicate yPlus and wallShearStress registrations introduced by the diagnostics include. Normal runs now reuse yPlusWall/wallShear; standalone resampling creates its own producers. Existing geometry edits preserved. Tested dirty workspace files:
- `postprocess/sample_wall_profiles.py` SHA256 `cdd174bf06656bff914609a309fd12b0e4d7937339573789cf3503719ec372ad`
- `base_case/system/upstream_diagnostics` SHA256 `640b89987defe36a1109026373878c28bc33bb777803defa8db51a2c0520f96f`

- Passed: `simpleFoam -case /tmp/gci-yplus-check-7q3fjs3e/coarse` exited 0 after two iterations, with wall outputs. Temporary copy of current coarse case; endTime and writeInterval set to 2 only in that copy.
- Passed: `simpleFoam -case /tmp/gci-yplus-check-7q3fjs3e/coarse -postProcess -latestTime -dict /tmp/gci-yplus-check-7q3fjs3e/coarse/system/upstreamDiagnosticsDict` exited 0; six profile files and two wall surface files present and nonempty. Logs: `/tmp/gci-yplus-check-7q3fjs3e/coarse/log.verify-solver` and `log.verify-resample`.
- Initial check exposed the additional wallShearStress duplicate; fixed before passing checks. Initial resampling check lacked saved p because writeInterval was 2000; repeated with writeInterval 2.
- Full convergence and mesh quality for the longer geometry: not verified. No production solve or data export performed.

## 2026-09-21 Allrun initialization speedup

- Tested `Allrun.sh`, SHA256 `bdbb98ee60f015869aa0211f01ef93fa3abecf3e9462774196cf922982c3c644`.
- Tested `base_case/system/controlDict`, SHA256 `72e6dd2bae52ea0953a2568c8ffe30d848b483219b6ddcbf047dda616524c4c5`.
- Passed: `sh -n Allrun.sh`; scoped `git diff --check`.
- Passed: temporary full Allrun with endTime/writeInterval 2, all three levels, renumberMesh -overwrite, coarse-to-medium and medium-to-fine mapFields interpolate, solvers and exports. Logs under `/tmp/gci-allrun-check-cojc5kww` (`log.test`, each run/log.renumberMesh and log.mapFields).
- Passed: standalone medium with temporarily absent coarse source used base initial fields and completed (`log.test-no-source`).
- writePrecision 12 removed wedge-rounding warnings encountered with precision 8. Final mapping logs have 51 lines each. Coarse bandwidth 4343 -> 54.
- Full convergence/runtime speedup and different-geometry fallback execution not verified. No production run or existing data changed.

## R1–R4: mesh revision, parallel workflow and convergence acceptance

- `Allrun.sh` SHA256 `5212fb1069e4217f0b9190718d814b2055425b9c9eef91480333840f0ba6cfff`.
- `case_design/case_design.py` SHA256 `73ec6a9b513dec1f093f3445c08017558c249876f2ad938abd7639558fe35360`.
- `base_case/system/controlDict` SHA256 `335f2424a555c784bf341d0d5775759fff0fa9c0096b95f06b7aeab770a4faf9`.
- `base_case/system/fvSolution` SHA256 `dabc6b875f90a9dbe31bb6eb2c40e797a9ec009307acec0c346b2248472fd73c`.
- `base_case/system/convergence_controls` SHA256 `6a88657d4be0fcc1b44daebdd7abd5b6bae22651fbb7ae1e6e17cf1178bcda4b`.
- `base_case/system/decomposeParDict` SHA256 `362facc7b6033f0e9fb92b5c6203e1066bd3a4f93fbc79f640d4a5cac51d2680`.
- `base_case/system/foam_functions` SHA256 `ddd0c0620e4dec376e4c3126c6f038b1c2a491c31d2e34dadce4ef3181157fe0`.

- R1 passed: existing coarse/medium logs ended at 10000 without SIMPLE convergence; mapFields logs recorded sourceTime 10000 (not 1000). Coarse max cell-centre |Uz| 1.317e-13 m/s despite normalized residual .01416. High-aspect cells occupy x=6.254–8.733 m, near r=.05 m, axial sizes .179–.236 m.
- R2 passed: generated actual meshes with 7200/16033/36000 cells. blockMesh, renumberMesh -overwrite, checkMesh all pass; maximum aspect ratios 466.516/468.210/472.277, zero nonorthogonality, max skewness .33171. Outputs under `/tmp/gci-parallel-check-lz4r2hux`.
- R3 passed: `sh /tmp/gci-parallel-check-lz4r2hux/Allrun.sh` in OpenFOAM environment using 8 MPI processes, relaxed test-only residual thresholds and endTime20. All three stop via grouped component convergence at iteration3, reconstruct that iteration, map coarse->medium->fine, export 1440 velocity samples plus wall data per level. log.parallel-test records the full workflow. Test-only thresholds do not establish physical convergence. MPI sandbox sockets were unavailable; the same test passed with approved escalated execution.
- R4 passed: test-only Ux threshold -1 while all other thresholds2, endTime4. Allrun exits1 at the cap, without accepting the other five satisfied conditions, reconstructing or exporting. Existing test CSV hashes unchanged. log.rejection-test records failure as expected.
- Passed: shell syntax, scoped diff hygiene, generated cell counts vs metadata, eight partitions per level, six exported profiles and synchronized solution time3.
- Production residual thresholds, same group structure: p1e-3, Ux/Uy2e-4, Uz.02, k/epsilon5e-5; strict standard SIMPLE criteria retained; iteration cap100000. Full production convergence/runtime speedup not verified. Existing production runs/data preserved; all meshes must be rerun before GCI using the updated metadata.

## M1–M3: physics-focused mesh revision

Tested in repository cwd, OpenFOAM v2412, dirty workspace; pre-existing plotting/workflow edits preserved.
Generator SHA256 883412ebcb0379191d2fb0a93e6d2ddab68c78a587d6673c07290a2b258ad59e.
Dictionary hashes coarse/medium/fine:
2f3a74413f2324205108a60e88f47921a712aba42f6cde78b68e2ee5bc0ce52e,
ce066f1facd113c00b73e5f32fb4906df2790675d8bce7e426ec14e9e256a0c2,
2d7e23807ad498c8810761b1d6c3e17a412b01df0e532a78fa2895c3f579ff8d.

| Criterion | Executed check / actual evidence | Verdict |
|---|---|---|
| Family and wall targets | Ran generator; metadata counts34720/78120/175770, targets1/.667/.444. Asserted sqrt count ratios1.5 and target ratios1.5. Generated dictionaries byte-identical to tested copies after final helper validation edits. | passed |
| Actual mesh connectivity and standard quality | In temporary copies: blockMesh, renumberMesh -overwrite, checkMesh for all3. All report Mesh OK. Max aspects332.423/332.855/333.142, nonorthogonality.802/.799/.797deg, skewness.332. | passed |
| Smooth redistribution and interface matching | /tmp/check_mesh_spacing.py reads actual generated points: coarse maximum adjacent radial factor1.15824, axial1.02827; refined levels lower. At x3/12 interface widths match .1762655mm, at x4/11 both1mm. Large-pipe wall height .6437581mm constant; all prescribed spacings scale by1.5. | passed |
| Visual layout | Opened /tmp/physics-mesh-46zuikzd/mesh_comparison.png from actual meshes. Checked lip clustering relaxes, outer wall stays refined, and broad step regions have axial resolution. | passed |
| Optional allGeometry determinant check | checkMesh -allGeometry -allTopology: small determinant warning on all3 (986/2264/5071cells). Previous coarse from HEAD also fails (999cells). New coarse min1.4233e-5 vs old1.0726e-5. No hidden threshold relaxation. | failed |
| Improved losses and measured y+ | No production solve; requires converged reruns of entire family. | not verified |

Logs and spacing report: /tmp/physics-mesh-46zuikzd/{coarse,medium,fine,previous-coarse}/log.*, spacing.txt and design.log.
Scoped git diff --check passes. Geometry, turbulence model, solver settings, production runs and datasets preserved.
The requested mesh revision is delivered with standard checks passing; optional determinant warning remains explicit, and no solution-accuracy claim is made.
