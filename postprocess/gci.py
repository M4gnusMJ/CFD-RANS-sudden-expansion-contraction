#!/usr/bin/env python3
"""GCI of the inlet friction factor from generated postprocess/data CSVs.

Run from any directory; optionally pass --data-dir for another CSV directory.
Cell counts and geometry come from case_design/case_params.json.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import root_scalar

from plot_pressure import extract_dp, read_csv


HERE = Path(__file__).resolve().parent
LEVELS = ("coarse", "medium", "fine")


def load_run_data(data_dir):
    """Load measured friction factors and unweighted wall-sample y+ medians."""
    params_path = HERE.parent / "case_design" / "case_params.json"
    with params_path.open() as stream:
        params = json.load(stream)
    cells = np.array([params["cells"][level] for level in LEVELS])
    domain_area = (
        (params["L1"] + params["L3"]) * params["d1"] / 2
        + params["L2"] * params["d2"] / 2
    )
    if not np.all(np.isfinite(cells)) or not np.all(cells > 0) or not np.all(np.diff(cells) > 0):
        raise ValueError("Cell counts must be positive and increase from coarse to fine.")
    if not np.isfinite(domain_area) or domain_area <= 0:
        raise ValueError("The representative domain area must be positive and finite.")

    phi, y_plus = [], []
    for level in LEVELS:
        samples = {}
        for quantity in ("axis", "yplus"):
            path = Path(data_dir) / f"{level}_{quantity}.csv"
            values = read_csv(path)
            if not values or not np.all(np.isfinite(values)):
                raise ValueError(f"Missing, empty, or non-finite run data: {path}")
            samples[quantity] = values
        measured = extract_dp(samples["axis"])
        if "f1_measured" not in measured:
            raise ValueError(f"{level}: insufficient pressure samples in the inlet fit window.")
        phi.append(measured["f1_measured"])
        y_plus.append(np.median([value for _, value in samples["yplus"]]))
    return cells, domain_area, np.array(phi), y_plus


def epsilon_of_monitoring_variable(phi_1, phi_2, phi_3):
    return np.array([phi_2 - phi_1, phi_3 - phi_2])


def q_p(r_21, r_32, s, p):
    return np.log((r_21**p - s) / (r_32**p - s))


def apparent_order_root_eq(p, r_21, r_32, phi_1, phi_2, phi_3):
    """Order equation with grids numbered from coarse (1) to fine (3)."""
    epsilon_21, epsilon_32 = epsilon_of_monitoring_variable(
        phi_1, phi_2, phi_3
    )
    s = np.sign(epsilon_32 / epsilon_21)
    return (
        p * np.log(r_32)
        + q_p(r_21, r_32, s, p)
        - np.log(np.abs(epsilon_21 / epsilon_32))
    )


def main(data_dir=HERE / "data"):
    cells, domain_area, phi, y_plus = load_run_data(data_dir)
    phi_1, phi_2, phi_3 = phi
    differences = np.diff(phi)
    if np.any(phi == 0) or np.any(differences == 0):
        raise ValueError("GCI requires nonzero friction factors and nonzero grid differences.")
    if differences[0] * differences[1] < 0:
        raise ValueError("Friction factors oscillate across grids; monotonic GCI is not applicable.")

    # Representative cell size for a two-dimensional mesh, in metres.
    h = np.sqrt(domain_area / cells)
    r_21, r_32 = h[:-1] / h[1:]
    ratios = np.array([r_21, r_32])

    solution = root_scalar(
        apparent_order_root_eq,
        args=(r_21, r_32, phi_1, phi_2, phi_3),
        x0=2.0,
        x1=2.1,
    )
    if not solution.converged or not np.isfinite(solution.root) or solution.root <= 0:
        raise RuntimeError("Could not determine a positive apparent order of convergence.")
    p = solution.root

    # Each pair uses its finer-grid value as the relative-error reference.
    refinement = ratios**p
    phi_ext = (refinement * phi[1:] - phi[:-1]) / (refinement - 1)
    approximate_relative_error = np.abs((phi[1:] - phi[:-1]) / phi[1:])
    extrapolated_relative_error = np.abs((phi_ext - phi[1:]) / phi_ext)
    gci_fine = 1.25 * approximate_relative_error / (refinement - 1)

    results = pd.DataFrame(
        {
            "mesh": LEVELS,
            "n": cells,
            "h [mm]": 1000 * h,
            "phi = X_n": phi,
            "r": ["", f"r_21 = {r_21:.3f}", f"r_32 = {r_32:.3f}"],
            "yPlus med.": y_plus,
            "phi_ext": ["", *[f"{value:.6f}" for value in phi_ext]],
            "e_a": ["", *[f"{100 * value:.2f}%" for value in approximate_relative_error]],
            "e_ext": ["", *[f"{100 * value:.2f}%" for value in extrapolated_relative_error]],
            "GCI fine": ["", *[f"{100 * value:.3f}%" for value in gci_fine]],
            "p": ["", "", f"{p:.2f}"],
        }
    )

    print(f"Run data: {Path(data_dir).resolve()}")
    print("phi = inlet Darcy friction factor extracted from the axis pressure slope")
    print(f"Grid refinement factors: r_21 = {r_21:.6f}, r_32 = {r_32:.6f}")
    print(f"Apparent order of convergence: {p:.6f}")
    print("GCI fine refers to the finer mesh in each pair (rows 2 and 3).")
    print(results.to_string(index=False, formatters={
        "h [mm]": "{:.5f}".format,
        "phi = X_n": "{:.5f}".format,
        "yPlus med.": "{:.2f}".format,
    }))
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=HERE / "data")
    args = parser.parse_args()
    try:
        results = main(args.data_dir)
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        parser.exit(1, f"GCI error: {exc}\n")
    # Display `styled_results` in a notebook to render the formatted table.
    styled_results = results.style.format(
        {"h [mm]": "{:.5f}", "phi = X_n": "{:.5f}", "yPlus med.": "{:.2f}"}
    ).hide(axis="index")
