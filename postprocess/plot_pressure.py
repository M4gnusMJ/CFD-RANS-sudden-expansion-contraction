#!/usr/bin/env python3
"""
Group 22 - post-processing: OpenFOAM vs. Bernoulli-with-losses.

Reads whatever mesh levels exist in data/ and produces:
  - pressure along the symmetry axis, all levels + the analytical curve
  - the extrapolated linear friction lines and the dp read off at each step
    (exactly the construction in the exercise sheet's figure)
  - y+ along the wall for every level
  - the extracted dp values printed so they can go straight into gci.py

Data contract - each run drops two CSVs so nobody is blocked on anybody:
  data/<level>_axis.csv    columns: x, p          (p is kinematic, m2/s2)
  data/<level>_yplus.csv   columns: x, yPlus
where <level> is coarse | medium | fine.

Run with no data at all and it still draws the analytical curve, so the script
can be reviewed and tested before the first simulation finishes.

Usage:  python3 plot_pressure.py
"""

import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

# --- geometry / flow ---------------------------------------------------------
# Read from case_design/case_params.json so there is ONE place to change a
# parameter. The fallbacks below are only used if that file is missing; if you
# see the warning, re-run case_design.py rather than editing numbers here.
_DEFAULTS = {"x_exp": 3.0, "x_con": 9.0, "Ltot": 12.0, "d1": 0.10, "d2": 0.20,
             "dp_exp": 0.75, "dp_con": -2.625, "f1": 0.02575, "U1": 2.0}
_PARAMS_PATH = os.path.join(HERE, "..", "case_design", "case_params.json")
try:
    with open(_PARAMS_PATH) as _fh:
        P = {**_DEFAULTS, **json.load(_fh)}
    _PARAMS_OK = True
except (OSError, ValueError):
    P = dict(_DEFAULTS)
    _PARAMS_OK = False

X_EXP, X_CON, L_TOT = P["x_exp"], P["x_con"], P["Ltot"]
D1, D2 = P["d1"], P["d2"]
H_STEP = (D2 - D1) / 2.0
LEVELS = ["coarse", "medium", "fine"]
COLORS = {"coarse": "#d6604d", "medium": "#4393c3", "fine": "#2166ac"}


def read_csv(path, cols=2):
    """Tolerant CSV reader: skips headers and comment lines."""
    if not os.path.exists(path):
        return None
    out = []
    with open(path) as fh:
        for row in csv.reader(fh):
            if not row or row[0].lstrip().startswith(("#", "x", "X")):
                continue
            try:
                out.append(tuple(float(row[i]) for i in range(cols)))
            except (ValueError, IndexError):
                continue
    return sorted(out) if out else None


def linfit(pts):
    """Least-squares slope/intercept. Pure python, no numpy needed."""
    n = len(pts)
    sx = sum(p[0] for p in pts)
    sy = sum(p[1] for p in pts)
    sxx = sum(p[0] ** 2 for p in pts)
    sxy = sum(p[0] * p[1] for p in pts)
    den = n * sxx - sx * sx
    if abs(den) < 1e-30:
        return 0.0, sy / n
    a = (n * sxy - sx * sy) / den
    b = (sy - a * sx) / n
    return a, b


def window(pts, x0, x1):
    return [p for p in pts if x0 <= p[0] <= x1]


def extract_dp(axis):
    """Extend the linear friction lines from each pipe section to the step and
    take the difference there - the red dashed-line construction in the handout.

    Exclusion zones: the flow already feels a step upstream of it, and downstream
    there is recirculation plus redevelopment. UP_EXCL / DN_EXCL below set how
    much is thrown away; widen them if the printed fit looks contaminated on the
    plot. The windows actually used are printed so they can be checked.
    """
    res = {}
    UP_EXCL_1, UP_EXCL_2 = 3 * D1, 3 * D2      # upstream of a step
    DN_EXCL = 20 * H_STEP                       # downstream of a step
    res["windows"] = {}

    up = window(axis, X_EXP - 12 * D1, X_EXP - UP_EXCL_1)
    dn = window(axis, X_EXP + DN_EXCL, X_CON - UP_EXCL_2)
    if len(up) >= 5 and len(dn) >= 5:
        a1, b1 = linfit(up)
        a2, b2 = linfit(dn)
        res["expansion"] = (a2 * X_EXP + b2) - (a1 * X_EXP + b1)
        res["fit_exp"] = ((a1, b1), (a2, b2))
        res["windows"]["expansion"] = (up[0][0], up[-1][0], dn[0][0], dn[-1][0])

    up = window(axis, X_EXP + DN_EXCL, X_CON - UP_EXCL_2)
    dn = window(axis, X_CON + DN_EXCL, L_TOT - UP_EXCL_1)
    if len(up) >= 5 and len(dn) >= 5:
        a1, b1 = linfit(up)
        a2, b2 = linfit(dn)
        res["contraction"] = (a2 * X_CON + b2) - (a1 * X_CON + b1)
        res["fit_con"] = ((a1, b1), (a2, b2))
        res["windows"]["contraction"] = (up[0][0], up[-1][0], dn[0][0], dn[-1][0])

    # friction factor from the measured slope in pipe 1:  f = -2 D (dp/dx) / U^2
    sec1 = window(axis, 0.35 * X_EXP, X_EXP - 6 * D1)
    if len(sec1) >= 5:
        a, _ = linfit(sec1)
        res["f1_measured"] = -2.0 * D1 * a / (P["U1"] ** 2)
    return res


def main():
    analytical = read_csv(os.path.join(HERE, "..", "case_design",
                                       "analytical_pressure.csv"))
    if analytical is None:
        analytical = read_csv(os.path.join(HERE, "analytical_pressure.csv"))

    found = {}
    for lv in LEVELS:
        ax = read_csv(os.path.join(DATA, f"{lv}_axis.csv"))
        yp = read_csv(os.path.join(DATA, f"{lv}_yplus.csv"))
        if ax or yp:
            found[lv] = {"axis": ax, "yplus": yp}

    print("Group 22 - post-processing")
    print("=" * 52)
    if _PARAMS_OK:
        print(f"parameters: d1={D1} d2={D2} Re1={P.get('Re1')} "
              f"steps at x={X_EXP} and x={X_CON}  (case_params.json)")
    else:
        print("WARNING: case_design/case_params.json not found - using fallback")
        print("         geometry. Run case_design.py so these stay in sync.")
    if analytical:
        print(f"analytical curve: {len(analytical)} points loaded")
    else:
        print("analytical curve: NOT FOUND (run case_design.py first)")
    if not found:
        print("\nNo simulation data in data/ yet.")
        print("Drop <level>_axis.csv and <level>_yplus.csv there as runs finish;")
        print("this script picks up whatever exists. Drawing analytical only.\n")

    # --- extracted quantities -------------------------------------------------
    table = {}
    for lv, d in found.items():
        if not d["axis"]:
            continue
        r = extract_dp(d["axis"])
        table[lv] = r
        print(f"\n{lv}:")
        if "expansion" in r:
            print(f"  dp/rho across expansion   = {r['expansion']:+.5f} m2/s2"
                  f"   (analytical {P['dp_exp']:+.5f})")
        if "contraction" in r:
            print(f"  dp/rho across contraction = {r['contraction']:+.5f} m2/s2"
                  f"   (analytical {P['dp_con']:+.5f})")
        for nm, w in r.get("windows", {}).items():
            print(f"  fit windows, {nm:<12} {w[0]:.2f}-{w[1]:.2f} m  and  "
                  f"{w[2]:.2f}-{w[3]:.2f} m")
        if "f1_measured" in r:
            print(f"  f1 from measured slope    = {r['f1_measured']:.5f}"
                  f"        (Haaland {P['f1']:.5f})")
        if d["yplus"]:
            ys = [p[1] for p in d["yplus"]]
            print(f"  y+ range                  = {min(ys):.3f} .. {max(ys):.3f}"
                  f"   (must stay <~ 1)")

    if len(table) == 3 and all("contraction" in table[l] for l in LEVELS):
        print("\nAll three levels present -> feed these into gci.py:")
        print("  contraction:", ", ".join(
            f"{l}={table[l]['contraction']:.5f}" for l in ["fine", "medium", "coarse"]))
        print("  expansion:  ", ", ".join(
            f"{l}={table[l].get('expansion', float('nan')):.5f}"
            for l in ["fine", "medium", "coarse"]))

    # --- plots ----------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nmatplotlib missing - numbers above are still valid.")
        return

    nplots = 2 if any(d.get("yplus") for d in found.values()) else 1
    fig, axes = plt.subplots(nplots, 1, figsize=(10, 5 * nplots), squeeze=False)
    ax = axes[0][0]

    if analytical:
        ax.plot([p[0] for p in analytical], [p[1] for p in analytical],
                "k--", lw=1.6, label="Bernoulli with losses", zorder=5)
    for lv, d in found.items():
        if d["axis"]:
            ax.plot([p[0] for p in d["axis"]], [p[1] for p in d["axis"]],
                    lw=1.8, color=COLORS[lv], label=f"simpleFoam, {lv}")

    # the dashed extrapolation construction, drawn for the finest level available
    for lv in ["fine", "medium", "coarse"]:
        if lv in table and "fit_con" in table[lv]:
            for key, xs in (("fit_exp", (0.3 * X_EXP, X_CON)),
                            ("fit_con", (X_EXP, L_TOT))):
                if key not in table[lv]:
                    continue
                for (a, b) in table[lv][key]:
                    ax.plot(xs, [a * xs[0] + b, a * xs[1] + b],
                            ":", color="#b2182b", lw=0.9, alpha=0.7, zorder=1)
            break

    ax.axvline(X_EXP, color="#4d9221", ls="--", lw=1, alpha=0.8)
    ax.axvline(X_CON, color="#b2182b", ls="--", lw=1, alpha=0.8)
    ax.set_xlabel("x  [m]")
    ax.set_ylabel(r"$p/\rho$  [m$^2$/s$^2$]")
    ax.set_title("Pressure along the symmetry axis - Group 22")
    ax.grid(alpha=0.3)
    ax.legend(frameon=False, fontsize=9)

    if nplots == 2:
        ax2 = axes[1][0]
        for lv, d in found.items():
            if d["yplus"]:
                ax2.plot([p[0] for p in d["yplus"]], [p[1] for p in d["yplus"]],
                         lw=1.5, color=COLORS[lv], label=lv)
        ax2.axhline(1.0, color="k", ls="--", lw=1, label=r"$y^+ = 1$")
        ax2.axvline(X_EXP, color="#4d9221", ls="--", lw=1, alpha=0.6)
        ax2.axvline(X_CON, color="#b2182b", ls="--", lw=1, alpha=0.6)
        ax2.set_xlabel("x  [m]")
        ax2.set_ylabel(r"$y^+$")
        ax2.set_title(r"$y^+$ along the wall - must stay below ~1 everywhere")
        ax2.grid(alpha=0.3)
        ax2.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    out = os.path.join(HERE, "pressure_comparison.png")
    fig.savefig(out, dpi=160)
    print(f"\nWrote {os.path.basename(out)}")


if __name__ == "__main__":
    sys.exit(main())
