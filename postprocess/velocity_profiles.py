#!/usr/bin/env python3
"""Investigate upstream development and wall behaviour from exported data only."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import trapezoid

HERE = Path(__file__).resolve().parent
LEVELS = ("coarse", "medium", "fine")
COLORS = {"coarse": "#d6604d", "medium": "#e69f00", "fine": "#2166ac"}


def scale_profile(profile, wall, station, nu):
    """Use local wall-face axial shear, not a developed-flow pressure assumption."""
    radius, x = station["wall_radius"], station["x"]
    faces = wall[np.isclose(wall.r, radius, rtol=1e-5, atol=1e-9)]
    faces = faces.groupby("x", as_index=False).mean(numeric_only=True).sort_values("x")
    if len(faces) < 2 or not faces.x.min() <= x <= faces.x.max():
        raise ValueError(f"Wall faces do not bracket x={x} at radius {radius}")
    shear = np.interp(x, faces.x, faces.tau_x_over_rho)
    u_tau = np.sqrt(abs(shear))
    if not np.isfinite(u_tau) or u_tau <= 0 or not np.isfinite(nu) or nu <= 0:
        raise ValueError("Inner scaling requires finite, nonzero wall shear and positive viscosity")
    actual_radius = np.interp(x, faces.x, faces.r)
    result = profile.sort_values("r").copy()
    if len(result) < 3 or np.any(np.diff(result.r) <= 0):
        raise ValueError("Profiles need at least three unique radial coordinates")
    result["wall_distance"] = actual_radius - result.r
    if np.any(result.wall_distance <= 0) or np.any(result.r < 0):
        raise ValueError("Profile points must lie between axis and wall")
    result["r_over_R"] = result.r / actual_radius
    result["U_over_Ubulk"] = result.Ux / station["bulk_velocity"]
    result["yplus"] = result.wall_distance * u_tau / nu
    result["uplus"] = result.Ux / u_tau
    result["nut_over_nu"] = result.nut / nu
    result["u_tau"] = u_tau
    wall_yplus = np.interp(x, faces.x, faces.wall_yplus)
    result["below_first_cell"] = result.yplus < wall_yplus
    return result, wall_yplus, shear


def profile_change(target, neighbour):
    """RMS and maximum shape change, expressed as percent of nominal bulk speed."""
    eta = np.linspace(0, 1, 1001)
    def interpolate(frame):
        return np.interp(eta, np.r_[0., frame.r_over_R, 1.],
                         np.r_[frame.U_over_Ubulk.iloc[0], frame.U_over_Ubulk, 0.])
    change = interpolate(target) - interpolate(neighbour)
    rms = np.sqrt(trapezoid(2 * eta * change ** 2, eta))
    return 100 * rms, 100 * np.max(np.abs(change))


def analyse(data_dir):
    profiles, summary, changes = {}, [], []
    for level in LEVELS:
        meta = json.loads((data_dir / f"{level}_wall_metadata.json").read_text())
        samples = pd.read_csv(data_dir / f"{level}_velocity_profiles.csv")
        walls = pd.read_csv(data_dir / f"{level}_wall_diagnostics.csv")
        for frame in (samples, walls):
            if frame.empty or not np.isfinite(frame.select_dtypes(include="number")).all().all():
                raise ValueError(f"{level}: empty or non-finite exported data")
        for name, station in meta["stations"].items():
            selected = samples[samples.station == name]
            if not np.allclose(selected.x, station["x"]):
                raise ValueError(f"{name}: sampled x differs from metadata")
            scaled, wall_yplus, shear = scale_profile(selected, walls, station, meta["nu"])
            profiles[level, name] = scaled
            u_tau = scaled.u_tau.iloc[0]
            eta = np.r_[0., scaled.r_over_R, 1.]
            axial = np.r_[scaled.Ux.iloc[0], scaled.Ux, 0.]
            bulk = trapezoid(2 * eta * axial, eta)
            # Sampling inside the wall-to-first-cell interval is interpolation,
            # not independently resolved evidence of the viscous profile.
            viscous = scaled[(~scaled.below_first_cell) & (scaled.yplus <= 5)]
            log_layer = scaled[(scaled.yplus >= 30) & (scaled.wall_distance / station["wall_radius"] <= .2)]
            summary.append({
                "mesh": level, "station": name, "x [m]": station["x"],
                "solution_time": meta["solution_time"],
                "u_tau [m/s]": u_tau, "tau_x/rho [m2/s2]": shear,
                "wall_cell_yplus": wall_yplus,
                "Re_tau": station["wall_radius"] * u_tau / meta["nu"],
                "Darcy_f_from_wall": 8 * u_tau ** 2 / station["bulk_velocity"] ** 2,
                "integrated_Ubulk [m/s]": bulk,
                "bulk_error [%]": 100 * (bulk / station["bulk_velocity"] - 1),
                "max_radial_speed/Ubulk [%]": 100 * scaled.Ur.abs().max() / station["bulk_velocity"],
                "sublayer_RMS_relative [%]": 100 * np.sqrt(np.mean((viscous.uplus / viscous.yplus - 1) ** 2)),
                "log_layer_RMS_uplus": np.sqrt(np.mean((log_layer.uplus - np.log(log_layer.yplus) / .41 - 5) ** 2)),
                "sublayer_samples": len(viscous), "log_layer_samples": len(log_layer),
            })
        for step in ("expansion", "contraction"):
            target = profiles[level, f"{step}_3D"]
            for offset in (4, 2):
                rms, maximum = profile_change(target, profiles[level, f"{step}_{offset}D"])
                changes.append({"mesh": level, "step": step, "comparison": f"3D vs {offset}D",
                                "area_weighted_RMS/Ubulk [%]": rms,
                                "maximum_difference/Ubulk [%]": maximum})
    return profiles, pd.DataFrame(summary), pd.DataFrame(changes)


def plot_wall_curve(ax, profile, label, color):
    resolved = profile[~profile.below_first_cell]
    ax.semilogx(resolved.yplus, resolved.uplus, color=color, label=label)


def style_wall_axes(ax, max_yplus, title):
    """Shared wall-law references and formatting for both velocity figures."""
    viscous = np.geomspace(.1, 25, 100)
    log_y = np.geomspace(1, max(31, max_yplus), 100)
    ax.plot(viscous, viscous, "k--", label="Viscous: u+ = y+")
    ax.plot(log_y, np.log(log_y) / .41 + 5, "k:", label="Log reference (κ=0.41, B=5)")
    ax.set(xlim=(.1, None), ylim=(0, 25), title=title,
           xlabel="y+ (distance from wall)", ylabel="u+")
    ax.grid(alpha=.25)
    ax.legend(fontsize=8, frameon=False, loc="upper left")


def plot(profiles, summary, figures_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for row, step in enumerate(("expansion", "contraction")):
        ax = axes[row]
        name = f"{step}_3D"
        for level in LEVELS:
            p = profiles[level, name]
            wall_yplus = summary[(summary.mesh == level) & (summary.station == name)].wall_cell_yplus.iloc[0]
            plot_wall_curve(ax, p, f"{level} (wall cell y+={wall_yplus:.2f})", COLORS[level])
        max_log = max(profiles[level, name].yplus.max() for level in LEVELS)
        x = summary[summary.station == name]["x [m]"].iloc[0]
        style_wall_axes(ax, max_log, f"3D before {step}: x = {x:g} m")
    fig.suptitle("Wall-scaled axial velocity")
    fig.tight_layout(rect=(0, 0, 1, .95))
    fig.savefig(figures_dir / "upstream_velocity_profiles.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(3, 2, figsize=(12, 13), sharey=True)
    station_colors = {2: "#d6604d", 3: "#e69f00", 4: "#2166ac"}
    for row, level in enumerate(LEVELS):
        for col, step in enumerate(("expansion", "contraction")):
            ax = axes[row, col]
            max_yplus = 0
            for offset in (2, 3, 4):
                name = f"{step}_{offset}D"
                p = profiles[level, name]
                station = summary[(summary.mesh == level) & (summary.station == name)].iloc[0]
                label = f"{offset}D (x={station['x [m]']:g} m, wall cell y+={station.wall_cell_yplus:.2f})"
                plot_wall_curve(ax, p, label, station_colors[offset])
                max_yplus = max(max_yplus, p.yplus.max())
            style_wall_axes(ax, max_yplus, f"Before {step} — {level}")
    # Keep the common scale while retaining any excursions in the sampled data.
    axes[0, 0].set_ylim(0, max(25, 1.05 * max(
        p.loc[~p.below_first_cell, "uplus"].max() for p in profiles.values())))
    fig.suptitle("Streamwise development of wall-scaled axial velocity")
    fig.tight_layout(rect=(0, 0, 1, .97))
    fig.savefig(figures_dir / "upstream_velocity_development.png", dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=HERE / "data")
    parser.add_argument("--output-dir", type=Path, default=HERE / "tables")
    parser.add_argument("--figures-dir", type=Path, default=HERE / "figs")
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()
    try:
        profiles, summary, changes = analyse(args.data_dir)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        summary.to_csv(args.output_dir / "upstream_wall_summary.csv", index=False)
        changes.to_csv(args.output_dir / "upstream_development.csv", index=False)
        pd.concat([p.assign(mesh=level) for (level, _), p in profiles.items()]).to_csv(
            args.output_dir / "upstream_scaled_profiles.csv", index=False)
        print(summary[summary.station.str.endswith("3D")].to_string(index=False))
        print(changes.to_string(index=False))
        if not args.no_plot:
            args.figures_dir.mkdir(parents=True, exist_ok=True)
            plot(profiles, summary, args.figures_dir)
        print(f"Tables: {args.output_dir}\nFigures: {args.figures_dir}")
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(1, f"Velocity diagnostics failed: {exc}\n"
                    "Export data first with sample_wall_profiles.py --resample in an OpenFOAM shell.\n")


if __name__ == "__main__":
    main()
