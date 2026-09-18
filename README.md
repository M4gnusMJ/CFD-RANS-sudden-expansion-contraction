# CFD Assignment 3 — Group 22

```
base_case/      OpenFOAM case
case_design/    mesh generator + analytical baseline
postprocess/    GCI and plotting scripts
```

## 1. Generate the meshes

```bash
python3 case_design/case_design.py
```

Writes `blockMeshDict.coarse/.medium/.fine` (16k / 36k / 81k cells) and
`design_summary.txt` with all parameters.

## 2. Run a case

```bash
cp case_design/blockMeshDict.medium base_case/system/blockMeshDict
blockMesh -case base_case
checkMesh -case base_case
simpleFoam -case base_case
```

Repeat for `.coarse` and `.fine`.

## 3. Post-process

Drop two CSVs per run into `postprocess/data/`:

```
<level>_axis.csv     columns: x, p
<level>_yplus.csv    columns: x, yPlus
```

`<level>` = coarse | medium | fine. Then:

```bash
python3 postprocess/plot_pressure.py
python3 postprocess/gci.py
```

`postprocess/foam_functions` has the OpenFOAM function objects that write these
files.
