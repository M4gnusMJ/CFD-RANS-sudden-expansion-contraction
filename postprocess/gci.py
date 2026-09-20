#!/usr/bin/env python3
"""GCI of inlet friction and contraction loss factors from generated data CSVs.

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
from minor_losses import calculate_losses


HERE = Path(__file__).resolve().parent
LEVELS = ("coarse", "medium", "fine")


def load_run_data(data_dir):
    """Load both indicators and wall-sample y+ medians from the same run CSVs."""
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

    friction, contraction, y_plus = [], [], []
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
        friction.append(measured["f1_measured"])
        losses, _ = calculate_losses(samples["axis"], level)
        contraction.append(next(row["K_sim"] for row in losses if row["step"] == "contraction"))
        y_plus.append(np.median([value for _, value in samples["yplus"]]))
    indicators = {
        "Inlet Darcy friction factor": np.array(friction),
        "Contraction loss coefficient K": np.array(contraction),
    }
    return cells, domain_area, indicators, y_plus


def epsilon_of_monitoring_variable(phi_1, phi_2, phi_3):
    return np.array([phi_2 - phi_1, phi_3 - phi_2])


def q_p(r_21, r_32, s, p):
    return np.log((r_21**p - s) / (r_32**p - s))


def apparent_order_root_eq(p, r_21, r_32, phi_1, phi_2, phi_3):
    """Absolute-value order formula, reindexed from coarse (1) to fine (3).

    A positive result from this formula does not by itself establish convergence.
    """
    epsilon_21, epsilon_32 = epsilon_of_monitoring_variable(
        phi_1, phi_2, phi_3
    )
    s = np.sign(epsilon_32 / epsilon_21)
    return p - np.abs(
        np.log(np.abs(epsilon_21 / epsilon_32)) - q_p(r_21, r_32, s, p)
    ) / np.log(r_32)


def indicator_table(cells, domain_area, phi, y_plus):
    """Report numerical GCI separately from the observed convergence status."""
    # Representative cell size for a two-dimensional mesh, in metres.
    h = np.sqrt(domain_area / cells)
    r_21, r_32 = h[:-1] / h[1:]
    ratios = np.array([r_21, r_32])

    results = pd.DataFrame(
        {
            "mesh": LEVELS,
            "n": cells,
            "h [mm]": 1000 * h,
            "phi = X_n": phi,
            "r": ["", f"r_21 = {r_21:.3f}", f"r_32 = {r_32:.3f}"],
            "yPlus med.": y_plus,
            "phi_ext": ["", "N/A", "N/A"],
            "e_a": ["", *[
                f"{100 * abs((fine - coarse) / fine):.2f}%" if fine != 0 else "N/A"
                for coarse, fine in zip(phi[:-1], phi[1:])
            ]],
            "e_ext": ["", "N/A", "N/A"],
            "GCI fine": ["", "N/A", "N/A"],
            "p": ["", "", "N/A"],
            "grid trend": ["", "", "Undetermined"],
        }
    )

    differences = np.diff(phi)
    if np.any(phi == 0) or np.any(differences == 0):
        return results, "GCI unavailable: zero indicator value or zero grid difference."
    if differences[0] * differences[1] < 0:
        results.loc[2, "grid trend"] = "Oscillatory"
        return results, "GCI unavailable: oscillatory grid sequence; monotonic extrapolation is not applicable."
    # The positive-order model has this lower bound as p approaches zero.
    difference_ratio = abs(differences[0] / differences[1])
    positive_order_supported = difference_ratio > np.log(r_21) / np.log(r_32)
    results.loc[2, "grid trend"] = (
        "Monotonic convergence" if positive_order_supported else "Non-convergent"
    )
    caveat = ""
    if not positive_order_supported:
        caveat = (
            "\nNumerical GCI and extrapolation only: this grid sequence does not support "
            "a positive-order convergence model, so GCI is not a validated uncertainty estimate. "
            f"Absolute changes: coarse→medium={abs(differences[0]):.6g}, "
            f"medium→fine={abs(differences[1]):.6g}."
        )

    solution = root_scalar(
        apparent_order_root_eq,
        args=(r_21, r_32, *phi),
        x0=2.0,
        x1=2.1,
    )
    if not solution.converged or not np.isfinite(solution.root) or solution.root <= 0:
        return results, "GCI unavailable: solver could not determine a positive convergence order."
    p = solution.root
    refinement = ratios**p
    phi_ext = (refinement * phi[1:] - phi[:-1]) / (refinement - 1)
    # Each pair uses its finer-grid value as the relative-error reference.
    approximate_relative_error = np.abs(differences / phi[1:])
    gci_fine = 1.25 * approximate_relative_error / (refinement - 1)
    results.loc[1:, "phi_ext"] = [f"{value:.6f}" for value in phi_ext]
    results.loc[1:, "e_ext"] = [
        f"{100 * abs((ext - fine) / ext):.2f}%" if ext != 0 else "N/A"
        for ext, fine in zip(phi_ext, phi[1:])
    ]
    results.loc[1:, "GCI fine"] = [f"{100 * value:.3f}%" for value in gci_fine]
    results.loc[2, "p"] = f"{p:.2f}"
    return results, f"Apparent order (absolute-value formula): {p:.6f}" + caveat


def main(data_dir=HERE / "data", output_dir=HERE / "tables"):
    cells, domain_area, indicators, y_plus = load_run_data(data_dir)

    print(f"Run data: {Path(data_dir).resolve()}")
    print("Both indicators use axis CSVs; y+ medians use wall CSVs.")
    print("Mesh counts and geometry: case_design/case_params.json")
    print("Contraction K uses extrapolated step pressures, small-pipe bulk velocity and alpha = 1.")
    print("GCI fine refers to the finer mesh in each pair (rows 2 and 3).")
    tables = []
    for name, phi in indicators.items():
        table, status = indicator_table(cells, domain_area, phi, y_plus)
        print("\n" + "─" * 110)
        print(f"Indicator: {name}")
        print(status)
        print(table.to_string(index=False, formatters={
            "h [mm]": "{:.5f}".format,
            "phi = X_n": "{:.6f}".format,
            "yPlus med.": "{:.2f}".format,
        }))
        table.insert(0, "indicator", name)
        table["status"] = status
        tables.append(table)
    results = pd.concat(tables, ignore_index=True)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / "gci.csv"
    results.to_csv(output, index=False)
    print(f"Saved {output}")
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=HERE / "data")
    parser.add_argument("--output-dir", type=Path, default=HERE / "tables", help="Result table directory")
    args = parser.parse_args()
    try:
        results = main(args.data_dir, args.output_dir)
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        parser.exit(1, f"GCI error: {exc}\n")
    # Display `styled_results` in a notebook to render the formatted table.
    styled_results = results.style.format(
        {"h [mm]": "{:.5f}", "phi = X_n": "{:.5f}", "yPlus med.": "{:.2f}"}
    ).hide(axis="index")
