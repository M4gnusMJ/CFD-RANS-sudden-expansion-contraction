#!/usr/bin/env python3
"""Estimate minor-loss K from developed pressure fits evaluated at each step.

Uses generated data/<level>_axis.csv and the fits in plot_pressure.py.
Pressure is kinematic (p/rho). For the horizontal pipe, with kinetic-energy
correction factors of one, as in case_design.py:
    K = [2 * (p_up/rho - p_down/rho) + U_up**2 - U_down**2] / U_ref**2
U_ref is the small-pipe bulk velocity: upstream for expansion, downstream for
contraction. Axis pressure is used as an approximation to section pressure.
"""

import argparse
import csv
import math
from pathlib import Path

import plot_pressure as pressure


HERE = Path(__file__).resolve().parent


def loss_coefficient(p_up, p_down, u_up, u_down, u_ref):
    """Return dimensionless K using kinematic pressures and bulk velocities."""
    if not all(math.isfinite(v) for v in (p_up, p_down, u_up, u_down, u_ref)):
        raise ValueError("Pressures and velocities must be finite.")
    if u_ref <= 0:
        raise ValueError("The reference velocity must be positive.")

    delta_p = p_up - p_down

    # Bernoulli: remove the pressure change associated with bulk acceleration.
    return (2 * delta_p + u_up**2 - u_down**2) / u_ref**2


def calculate_losses(axis, level):
    """Reuse the developed-region fits, evaluating both sides at the same x."""
    if not axis or not all(math.isfinite(v) for point in axis for v in point):
        raise ValueError(f"{level}: empty or non-finite axis pressure data.")
    fits = pressure.extract_dp(axis)
    params = pressure.P
    u_small, u_large = params["U1"], params["U2"]
    rows = []
    for step, fit_key, x, u_up, u_down, analytical in (
        ("expansion", "fit_exp", params["x_exp"], u_small, u_large, params["K_SE"]),
        ("contraction", "fit_con", params["x_con"], u_large, u_small, params["K_SC"]),
    ):
        if fit_key not in fits:
            raise ValueError(f"{level}: insufficient developed-region samples for {step}.")
        (a_up, b_up), (a_down, b_down) = fits[fit_key]
        p_up, p_down = a_up * x + b_up, a_down * x + b_down
        k = loss_coefficient(p_up, p_down, u_up, u_down, u_small)
        w = fits["windows"][step]
        rows.append({
            "mesh": level, "step": step, "x_m": x,
            "p_up_over_rho": p_up, "p_down_over_rho": p_down,
            "dp_down_minus_up_over_rho": p_down - p_up,
            "loss_m2_s2": k * u_small**2 / 2,
            "U_ref_m_s": u_small, "K_sim": k, "K_analytical": analytical,
            "fit_up_start_m": w[0], "fit_up_end_m": w[1],
            "fit_down_start_m": w[2], "fit_down_end_m": w[3],
        })
    return rows, fits


def plot_losses(runs, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(len(runs), 1, figsize=(10, 4 * len(runs)), squeeze=False)
    for ax, (level, axis, rows, fits) in zip(axes[:, 0], runs):
        ax.plot(*zip(*axis), color=pressure.COLORS[level], label="Simulation axis pressure")
        for row, key in zip(rows, ("fit_exp", "fit_con")):
            x = row["x_m"]
            for (a, b), start, end in zip(
                fits[key],
                (row["fit_up_start_m"], x),
                (x, row["fit_down_end_m"]),
            ):
                ax.plot([start, end], [a * start + b, a * end + b], "k--", lw=1)
            p_up, p_down = row["p_up_over_rho"], row["p_down_over_rho"]
            ax.annotate("", xy=(x, p_up), xytext=(x, p_down),
                        arrowprops={"arrowstyle": "<->", "color": "#b2182b"})
            ax.annotate(f"{row['step']}\nK = {row['K_sim']:.4f}",
                        xy=(x, max(p_up, p_down)), xytext=(0, 12),
                        textcoords="offset points", ha="center", va="bottom")
            ax.axvline(x, color="grey", lw=0.7, alpha=0.4)
        ax.set(title=level, xlabel="x [m]", ylabel="p/rho [m²/s²]")
        ax.margins(y=0.25)
        ax.grid(alpha=0.25)
        ax.legend(loc="lower left")
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=HERE / "data")
    parser.add_argument("--output-dir", type=Path, default=HERE / "tables", help="Result table directory")
    parser.add_argument("--figures-dir", type=Path, default=HERE / "figs", help="Plot directory")
    parser.add_argument("--no-plot", action="store_true", help="Only print and save the CSV table")
    args = parser.parse_args(argv)
    try:
        if not pressure._PARAMS_OK:
            raise ValueError("Generate case_design/case_params.json before calculating losses.")
        runs = []
        for level in pressure.LEVELS:
            path = args.data_dir / f"{level}_axis.csv"
            if not path.exists():
                print(f"Skipping {level}: {path} is missing")
                continue
            axis = pressure.read_csv(path)
            rows, fits = calculate_losses(axis, level)
            runs.append((level, axis, rows, fits))
        if not runs:
            raise ValueError(f"No axis CSVs found in {args.data_dir}")
        args.output_dir.mkdir(parents=True, exist_ok=True)
        output = args.output_dir / "minor_losses.csv"
        all_rows = [row for _, _, rows, _ in runs for row in rows]
        with output.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(all_rows[0]))
            writer.writeheader()
            writer.writerows(all_rows)
        print("K uses small-pipe bulk velocity; alpha = 1; pressures extrapolated to the step.")
        print(f"{'mesh':<8} {'step':<12} {'x [m]':>6} {'dp/rho [m²/s²]':>16} {'K simulation':>14} {'K analytical':>14}")
        for row in all_rows:
            print(f"{row['mesh']:<8} {row['step']:<12} {row['x_m']:6.2f} "
                  f"{row['dp_down_minus_up_over_rho']:16.6f} "
                  f"{row['K_sim']:14.6f} {row['K_analytical']:14.6f}")
        print(f"Saved {output}")
        if not args.no_plot:
            args.figures_dir.mkdir(parents=True, exist_ok=True)
            figure = args.figures_dir / "minor_losses.png"
            plot_losses(runs, figure)
            print(f"Saved {figure}")
    except (OSError, ValueError, KeyError, ImportError) as exc:
        parser.exit(1, f"Minor-loss calculation error: {exc}\n")


if __name__ == "__main__":
    main()
