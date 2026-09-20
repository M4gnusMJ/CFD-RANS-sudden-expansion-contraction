#!/usr/bin/env python3
"""GCI of inlet friction, contraction losses and the axial pressure profile.

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


def pressure_profile_table(data_dir, cells, domain_area):
    """Celik profile procedure: mean local order, then pointwise fine-grid GCI.

    All three profiles must be sampled at the same locations and pressure datum.
    No interpolation is used, to avoid adding interpolation error near the steps.
    """
    profiles = []
    for level in LEVELS:
        path = Path(data_dir) / f"{level}_axis.csv"
        samples = read_csv(path)
        if samples is None or len(samples) < 2 or not np.all(np.isfinite(samples)):
            raise ValueError(f"Missing, insufficient or non-finite axis data: {path}")
        profile = np.asarray(samples)
        if np.any(np.diff(profile[:, 0]) <= 0):
            raise ValueError(f"Axis coordinates must be unique and increasing: {path}")
        profiles.append(profile)
    x = profiles[-1][:, 0]
    if any(a.shape != profiles[-1].shape or not np.allclose(
            a[:, 0], x, rtol=0, atol=1e-10) for a in profiles[:-1]):
        raise ValueError("GCI requires all axis CSVs sampled at the same x coordinates; "
                         "regenerate the profiles on a common sampling line.")
    h = np.sqrt(domain_area / np.asarray(cells))
    r_cm, r_mf = h[:-1] / h[1:]
    if not np.all(np.isfinite(h)) or min(r_cm, r_mf) <= 1:
        raise ValueError("GCI requires positive cell sizes decreasing from coarse to fine.")
    phi = np.array([a[:, 1] for a in profiles])
    differences = np.diff(phi, axis=0)
    orders = np.full(len(x), np.nan)
    trends, order_status = [], []
    for i, (d_cm, d_mf) in enumerate(differences.T):
        # Only roundoff-scale differences are treated as indistinguishable.
        tol = 16 * np.finfo(float).eps * max(1.0, np.max(np.abs(phi[:, i])))
        if abs(d_cm) <= tol or abs(d_mf) <= tol:
            trends.append("Undetermined")
            order_status.append("Zero/roundoff grid difference; excluded from mean order")
            continue
        s = np.sign(d_cm) * np.sign(d_mf)
        log_ratio = np.log(abs(d_cm)) - np.log(abs(d_mf))
        trends.append("Oscillatory" if s < 0 else (
            "Monotonic convergence" if log_ratio > np.log(np.log(r_cm) / np.log(r_mf))
            else "Non-convergent"))

        def equation(order):
            # Stable log(r**p - s), avoiding overflow and cancellation at p ~ 0.
            a, b = order * np.log([r_cm, r_mf])
            if s > 0:
                q = a - b + np.log(-np.expm1(-a)) - np.log(-np.expm1(-b))
            else:
                q = np.logaddexp(a, 0) - np.logaddexp(b, 0)
            return order - abs(log_ratio - q) / np.log(r_mf)

        lower, upper = 1e-10, 2.0
        while equation(upper) < 0 and upper < 1024:
            upper *= 2
        if equation(lower) < 0 <= equation(upper):
            solution = root_scalar(equation, bracket=(lower, upper), method="brentq")
            if solution.converged:
                orders[i] = solution.root
        order_status.append("Included" if np.isfinite(orders[i]) else
                            "No resolved positive order; excluded from mean order")
    valid = np.isfinite(orders)
    if not valid.any():
        raise ValueError("No resolved local apparent orders; cannot form pressure p_average.")
    p_average = orders[valid].mean()
    # Absolute GCI avoids dividing by gauge pressure, which may be zero.
    uncertainty = 1.25 * np.abs(differences[1]) / np.expm1(p_average * np.log(r_mf))
    relative = np.full(len(x), np.nan)
    np.divide(100 * uncertainty, np.abs(phi[2]), out=relative, where=phi[2] != 0)
    return pd.DataFrame({
        "x [m]": x,
        "p_coarse/rho [m2/s2]": phi[0],
        "p_medium/rho [m2/s2]": phi[1],
        "p_fine/rho [m2/s2]": phi[2],
        "p_local": orders, "p_average": p_average,
        "grid trend": trends, "order status": order_status,
        "GCI absolute [m2/s2]": uncertainty, "GCI relative [%]": relative,
        "lower/rho [m2/s2]": phi[2] - uncertainty,
        "upper/rho [m2/s2]": phi[2] + uncertainty,
    })


def plot_pressure_gci(table, figures_dir):
    """Plot every profile point, with a readable subset of the GCI error bars."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = table["x [m]"].to_numpy()
    fine = table["p_fine/rho [m2/s2]"].to_numpy()
    uncertainty = table["GCI absolute [m2/s2]"].to_numpy()
    # Include the largest bar even when it falls between displayed samples.
    indices = np.unique(np.r_[np.linspace(0, len(x) - 1, min(65, len(x)), dtype=int),
                              np.argmax(uncertainty)])
    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.plot(x, table["p_coarse/rho [m2/s2]"], color="#999999", lw=1,
            linestyle=":", label="Coarse")
    ax.plot(x, table["p_medium/rho [m2/s2]"], color="#e69f00", lw=1,
            linestyle="--", label="Medium")
    ax.fill_between(x, fine - uncertainty, fine + uncertainty, color="#2166ac", alpha=.12)
    ax.plot(x, fine, color="#2166ac", lw=1.5, label="Fine")
    ax.errorbar(x[indices], fine[indices], yerr=uncertainty[indices], fmt="none",
                ecolor="#2166ac", elinewidth=.9, capsize=2, label="Fine-grid GCI (Fs = 1.25)")
    p_average = table["p_average"].iloc[0]
    oscillatory = 100 * table["grid trend"].eq("Oscillatory").mean()
    nonconvergent = 100 * table["grid trend"].eq("Non-convergent").mean()
    ax.set(xlabel="x [m]", ylabel=r"$p/\rho$ [m$^2$/s$^2$]",
           title=f"Axial pressure with GCI error bars — average order = {p_average:.3f}")
    ax.grid(alpha=.25)
    ax.legend(frameon=False)
    fig.text(.5, .02, f"Oscillatory: {oscillatory:.1f}% | Non-convergent: {nonconvergent:.1f}% | "
             "Bars use the global mean order; local convergence is not guaranteed.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .05, 1, 1))
    figures_dir = Path(figures_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)
    output = figures_dir / "pressure_gci.png"
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(f"Saved {output}")


def main(data_dir=HERE / "data", output_dir=HERE / "tables",
         figures_dir=HERE / "figs", no_plot=False):
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
    profile = pressure_profile_table(data_dir, cells, domain_area)
    profile_output = output_dir / "pressure_gci.csv"
    profile.to_csv(profile_output, index=False)
    print(f"Saved {profile_output}")
    orders = profile["p_local"].dropna()
    print(f"Pressure local order: {orders.min():.6f} to {orders.max():.6f}; "
          f"average = {orders.mean():.6f} ({len(orders)}/{len(profile)} points)")
    print("Pressure grid trends:", profile["grid trend"].value_counts().to_dict())
    print(f"Maximum absolute pressure GCI: {profile['GCI absolute [m2/s2]'].max():.6g} m2/s2")
    print("Profile bars use the mean order; they do not establish local convergence.")
    if not no_plot:
        plot_pressure_gci(profile, figures_dir)
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=HERE / "data")
    parser.add_argument("--output-dir", type=Path, default=HERE / "tables", help="Result table directory")
    parser.add_argument("--figures-dir", type=Path, default=HERE / "figs")
    parser.add_argument("--no-plot", action="store_true", help="Write tables only")
    args = parser.parse_args()
    try:
        results = main(args.data_dir, args.output_dir, args.figures_dir, args.no_plot)
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        parser.exit(1, f"GCI error: {exc}\n")
    # Display `styled_results` in a notebook to render the formatted table.
    styled_results = results.style.format(
        {"h [mm]": "{:.5f}", "phi = X_n": "{:.5f}", "yPlus med.": "{:.2f}"}
    ).hide(axis="index")
