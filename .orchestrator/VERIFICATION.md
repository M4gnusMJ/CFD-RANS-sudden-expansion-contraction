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
