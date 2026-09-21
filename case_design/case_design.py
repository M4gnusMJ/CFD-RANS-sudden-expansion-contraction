#!/usr/bin/env python3
"""
CFD Assignment 3 - Group 22
Case design + analytical baseline + parametric blockMeshDict generator.

Everything that has to be DECIDED before anyone runs OpenFOAM is here:

  1. Reynolds numbers / pipe dimensions / bulk velocities
  2. Entrance lengths  -> the three pipe-section lengths
  3. Friction factors  -> u_tau -> first-cell height for a target y+
  4. Minor-loss coefficients for the sudden expansion / contraction
  5. The analytical "Bernoulli with losses" pressure distribution p/rho(x)
     that the simulation is supposed to reproduce
  6. Three systematically refined blockMeshDict files (r = 1.5) for the GCI study

References (put these in the .bib):
  [1] Cengel & Cimbala, Fluid Mechanics: Fundamentals and Applications,
      Ch. 8 (entrance length, Colebrook/Haaland, minor losses Table 8-4).
  [2] Celik et al., "Procedure for Estimation and Reporting of Uncertainty Due to
      Discretization in CFD Applications", J. Fluids Eng. 130(7), 2008.  (GCI)

Usage:  python3 case_design.py
"""

import math
import os
import csv

OUT = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------------------
# 1. PRIMARY CHOICES  (the only things you really "pick" - everything else follows)
# ----------------------------------------------------------------------------
nu = 1.0e-5          # m^2/s   kinematic viscosity (matches transportProperties in the repo)
d1 = 0.10            # m       small-pipe diameter (inlet and outlet pipe)
d2 = 0.20            # m       large-pipe diameter (expansion section) -> area ratio 4
Re1 = 20000.0        # -       Reynolds number in the small pipe
y_plus_target = 1.0  # -       we RESOLVE the boundary layer (low-Re RANS model)
half_angle_deg = 2.0 # deg     wedge half-angle (the handout blockMeshDict uses 2 deg)

r1, r2 = d1 / 2.0, d2 / 2.0
U1 = Re1 * nu / d1                 # bulk velocity, small pipe
U2 = U1 * (d1 / d2) ** 2           # bulk velocity, large pipe (continuity)
Re2 = U2 * d2 / nu
AR = (d1 / d2) ** 2                # A1/A2 = 0.25


# ----------------------------------------------------------------------------
# 2. FRICTION FACTORS  (smooth pipe -> Haaland explicit form of Colebrook) [1]
# ----------------------------------------------------------------------------
def haaland_smooth(Re, eps_over_d=0.0):
    """Haaland (1983) explicit approximation to Colebrook, Darcy friction factor."""
    return (-1.8 * math.log10((eps_over_d / 3.7) ** 1.11 + 6.9 / Re)) ** -2


def colebrook_smooth(Re, eps_over_d=0.0, f0=0.02):
    """Colebrook-White solved by fixed-point iteration (sanity check on Haaland)."""
    f = f0
    for _ in range(200):
        f = (-2.0 * math.log10(eps_over_d / 3.7 + 2.51 / (Re * math.sqrt(f)))) ** -2
    return f


f1, f2 = haaland_smooth(Re1), haaland_smooth(Re2)
u_tau1 = U1 * math.sqrt(f1 / 8.0)   # tau_w/rho = (f/8) U^2  ->  u_tau = U sqrt(f/8)
u_tau2 = U2 * math.sqrt(f2 / 8.0)

# wall-normal distance of the FIRST CELL CENTRE for the target y+
y_c1 = y_plus_target * nu / u_tau1
y_c2 = y_plus_target * nu / u_tau2
dw1 = 2.0 * y_c1                    # -> first cell HEIGHT (centre sits at half the height)
dw2 = 2.0 * y_c2


# ----------------------------------------------------------------------------
# 3. ENTRANCE LENGTHS -> SECTION LENGTHS
#    turbulent hydrodynamic entrance length, Cengel & Cimbala: L_e/D ~ 4.4 Re^(1/6)
#    after a sudden expansion the flow first has to reattach: x_r ~ 8-10 step heights
# ----------------------------------------------------------------------------
def entrance_length(Re, D):
    return 4.4 * Re ** (1.0 / 6.0) * D


h_step = (d2 - d1) / 2.0            # step height
Le1 = entrance_length(Re1, d1)
Le2 = entrance_length(Re2, d2)
x_reattach = 10.0 * h_step

L1 = math.ceil((1.10 * Le1) / (5 * d1)) * 5 * d1          # inlet pipe, 10 % margin
L2 = 9  # expansion pipe
L3 = L1                                                   # outlet pipe, same as inlet
Ltot = L1 + L2 + L3
x_exp = L1                 # axial position of the sudden expansion
x_con = L1 + L2            # axial position of the sudden contraction


# ----------------------------------------------------------------------------
# 4. MINOR LOSSES  (Cengel & Cimbala Table 8-4) [1]
#    sudden expansion:   K_SE = (1 - A1/A2)^2 , based on the UPSTREAM velocity V1
#    sudden contraction: K_SC = 0.5 (1 - A2/A1), based on the DOWNSTREAM velocity
# ----------------------------------------------------------------------------
K_SE = (1.0 - AR) ** 2
K_SC = 0.5 * (1.0 - AR)
K_SC_alt = 0.42 * (1.0 - AR)        # alternative correlation - quote as uncertainty

# static-pressure jumps (kinematic pressure p/rho, as OpenFOAM's incompressible p)
dp_exp = (U1 ** 2 - U2 ** 2) / 2.0 - K_SE * U1 ** 2 / 2.0   # net RISE at the expansion
dp_con = (U2 ** 2 - U1 ** 2) / 2.0 - K_SC * U1 ** 2 / 2.0   # net DROP at the contraction

# friction gradients, Darcy-Weisbach:  d(p/rho)/dx = -(f/D) U^2/2
grad1 = -(f1 / d1) * U1 ** 2 / 2.0
grad2 = -(f2 / d2) * U2 ** 2 / 2.0


def p_over_rho(x):
    """Analytical kinematic pressure along the axis, normalised to p = 0 at the outlet."""
    if x <= x_exp:
        p = grad1 * x
    elif x <= x_con:
        p = grad1 * x_exp + dp_exp + grad2 * (x - x_exp)
    else:
        p = grad1 * x_exp + dp_exp + grad2 * L2 + dp_con + grad1 * (x - x_con)
    return p - _P_OUT


_P_OUT = 0.0                       # set below so that p = 0 at the outlet
_P_OUT = p_over_rho(L1 + L2 + L3)


# ----------------------------------------------------------------------------
# 5. TURBULENT INLET BOUNDARY CONDITIONS
#    fully-developed pipe flow: I ~ 0.16 Re^(-1/8), l = 0.07 D (Cengel / standard practice)
# ----------------------------------------------------------------------------
I_inlet = 0.16 * Re1 ** (-1.0 / 8.0)
l_inlet = 0.07 * d1
k_inlet = 1.5 * (U1 * I_inlet) ** 2
Cmu = 0.09
eps_inlet = Cmu ** 0.75 * k_inlet ** 1.5 / l_inlet
omega_inlet = eps_inlet / (Cmu * k_inlet)


# ----------------------------------------------------------------------------
# 6. GRADED-MESH HELPERS
# ----------------------------------------------------------------------------
def expansion_k(L, n, d_small, lo=1.0 + 1e-9, hi=1.5):
    """Per-cell growth factor k of a geometric distribution of n cells over length L
    whose smallest cell has size d_small.  Solves d_small (k^n - 1)/(k - 1) = L."""
    if n * d_small >= L:
        return 1.0
    for _ in range(300):
        k = 0.5 * (lo + hi)
        s = d_small * (k ** n - 1.0) / (k - 1.0) if k > 1.0 + 1e-12 else d_small * n
        if s < L:
            lo = k
        else:
            hi = k
    return 0.5 * (lo + hi)


def grading_one_sided(L, n, d_small, fine_at_start):
    """blockMesh expansion ratio (last cell / first cell) for a one-sided grading."""
    k = expansion_k(L, n, d_small)
    R = k ** (n - 1)
    return R if fine_at_start else 1.0 / R


def grading_two_sided(L, n, d_start, d_end):
    """Multi-grading string: fine at both ends of the block."""
    nh = max(2, n // 2)
    R1 = expansion_k(L / 2.0, nh, d_start) ** (nh - 1)
    R2 = expansion_k(L / 2.0, nh, d_end) ** (nh - 1)
    return "( (0.5 0.5 %.5f) (0.5 0.5 %.5f) )" % (R1, 1.0 / R2)


def first_cell_one_sided(L, n, d_small):
    return d_small


def grading_with_uniform_middle(L, n, d_small):
    """Graded end quarters with a uniform middle half of the cells.

    Match the last cell of each end section to the uniform cell size, avoiding
    the very long central cells of a full-length geometric progression.
    """
    n_end = max(2, n // 4)
    n_middle = n - 2 * n_end
    if n_middle < 1 or n * d_small >= L:
        raise ValueError("Axial grading needs a middle section and n*d_small < L")

    def end_length(d_max):
        ratio = d_max / d_small
        return sum(d_small * ratio ** (i / (n_end - 1)) for i in range(n_end))

    lo, hi = L / n, L / n_middle
    for _ in range(80):
        d_max = (lo + hi) / 2
        if 2 * end_length(d_max) + n_middle * d_max < L:
            lo = d_max
        else:
            hi = d_max
    d_max = (lo + hi) / 2
    end_fraction = end_length(d_max) / L
    ratio = d_max / d_small
    return ("( (%.12g %.12g %.12g) (%.12g %.12g 1) (%.12g %.12g %.12g) )"
            % (end_fraction, n_end / n, ratio,
               1 - 2 * end_fraction, n_middle / n,
               end_fraction, n_end / n, 1 / ratio))


# ----------------------------------------------------------------------------
# 7. blockMeshDict WRITER  (same wedge topology as the handout)
# ----------------------------------------------------------------------------
HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM
   \\\\    /   O peration     | Group 22 - Assignment 3
    \\\\  /    A nd           | GENERATED BY case_design.py - DO NOT EDIT BY HAND
     \\\\/     M anipulation  | mesh level: %s
\\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1;

"""


def write_blockmeshdict(path, level_name, n_ax1, n_ax2, n_ax3, n_r_core, n_r_ann,
                        dx_fine, dr_wall1, dr_wall2):
    th = math.radians(half_angle_deg)
    c, s = math.cos(th), math.sin(th)

    xs = [0.0, x_exp, x_con, Ltot]
    verts = []
    for x in xs:                                   # 0..11 : core pipe, r = r1
        verts.append((x, 0.0, 0.0))                # 3i   axis
        verts.append((x, r1 * c,  r1 * s))         # 3i+1 front (+z)
        verts.append((x, r1 * c, -r1 * s))         # 3i+2 back  (-z)
    verts.append((x_exp, r2 * c,  r2 * s))         # 12
    verts.append((x_exp, r2 * c, -r2 * s))         # 13
    verts.append((x_con, r2 * c,  r2 * s))         # 14
    verts.append((x_con, r2 * c, -r2 * s))         # 15

    # gradings
    gx1 = grading_one_sided(L1, n_ax1, dx_fine, fine_at_start=False)   # fine at expansion
    gx2 = grading_with_uniform_middle(L2, n_ax2, dx_fine)             # fine at steps, capped middle
    gx3 = grading_one_sided(L3, n_ax3, dx_fine, fine_at_start=True)    # fine at contraction
    gr_core = grading_one_sided(r1, n_r_core, dr_wall1, fine_at_start=False)
    # Bottom of the annulus is the shear layer off the step, NOT a wall, and it
    # does not have to match the core block: the shared face r = r1 is indexed
    # axially. Targeting the wall scale here forced expansion ratios > 1.4 and
    # made the grading differ between mesh levels, which breaks the systematic
    # refinement a GCI study needs. Use the shear-layer scale instead.
    gr_ann = grading_two_sided(r2 - r1, n_r_ann, dx_fine, dr_wall2)

    lines = [HEADER % level_name, "vertices\n(\n"]
    for v in verts:
        lines.append("    (%.12f %.12f %.12f)\n" % v)
    lines.append(");\n\nblocks\n(\n")
    lines.append("    hex (0 3 5 2 0 3 4 1) (%d %d 1) simpleGrading (%.5f %.5f 1)\n"
                 % (n_ax1, n_r_core, gx1, gr_core))
    lines.append("    hex (3 6 8 5 3 6 7 4) (%d %d 1) simpleGrading (%s %.5f 1)\n"
                 % (n_ax2, n_r_core, gx2, gr_core))
    lines.append("    hex (6 9 11 8 6 9 10 7) (%d %d 1) simpleGrading (%.5f %.5f 1)\n"
                 % (n_ax3, n_r_core, gx3, gr_core))
    lines.append("    hex (5 8 15 13 4 7 14 12) (%d %d 1) simpleGrading (%s %s 1)\n"
                 % (n_ax2, n_r_ann, gx2, gr_ann))
    lines.append(");\n\nedges\n(\n")
    for i, x in enumerate(xs):
        lines.append("    arc %d %d (%.12f %.12f 0)\n" % (3 * i + 1, 3 * i + 2, x, r1))
    lines.append("    arc 12 13 (%.12f %.12f 0)\n" % (x_exp, r2))
    lines.append("    arc 14 15 (%.12f %.12f 0)\n" % (x_con, r2))
    lines.append(""");

boundary
(
    inlet
    {
        type patch;
        faces ( (0 1 2 0) );
    }
    outlet
    {
        type patch;
        faces ( (9 11 10 9) );
    }
    upperWall
    {
        type wall;
        faces
        (
            (1 4 5 2)
            (4 12 13 5)
            (12 14 15 13)
            (7 8 15 14)
            (7 10 11 8)
        );
    }
    back
    {
        type wedge;
        faces ( (0 2 5 3) (3 5 8 6) (6 8 11 9) (5 13 15 8) );
    }
    front
    {
        type wedge;
        faces ( (0 3 4 1) (3 6 7 4) (6 9 10 7) (4 7 14 12) );
    }
    axis
    {
        type empty;
        faces ( (0 3 3 0) (3 6 6 3) (6 9 9 6) );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //
""")
    with open(path, "w") as fh:
        fh.writelines(lines)
    return (n_ax1 + n_ax2 + n_ax3) * n_r_core + n_ax2 * n_r_ann


# ----------------------------------------------------------------------------
# 8. THE THREE MESHES  (systematic refinement, r = 1.5 in every direction)
# ----------------------------------------------------------------------------
REF = 1.5
BASE = dict(n_ax1=200, n_ax2=250, n_ax3=200, n_r_core=40, n_r_ann=40,
            dx_fine=0.002, dr_wall1=dw1 / REF, dr_wall2=dw2 / REF)

levels = []
# The former medium resolution is now fine; retain the same refinement factor.
for name, p in [("coarse", 1.0 / REF**2), ("medium", 1.0 / REF), ("fine", 1.0)]:
    cfg = dict(
        n_ax1=int(round(BASE["n_ax1"] * p)),
        n_ax2=int(round(BASE["n_ax2"] * p)),
        n_ax3=int(round(BASE["n_ax3"] * p)),
        n_r_core=int(round(BASE["n_r_core"] * p)),
        n_r_ann=2 * int(round(BASE["n_r_ann"] * p / 2)),
        dx_fine=BASE["dx_fine"] / p,
        dr_wall1=BASE["dr_wall1"] / p,
        dr_wall2=BASE["dr_wall2"] / p,
    )
    path = os.path.join(OUT, "blockMeshDict.%s" % name)
    N = write_blockmeshdict(path, name, **cfg)
    levels.append((name, N, cfg))


# ----------------------------------------------------------------------------
# 9. REPORT
# ----------------------------------------------------------------------------
def banner(t):
    return "\n" + t + "\n" + "-" * len(t)


rep = []
rep.append("CFD Assignment 3 - Group 22 - CASE DESIGN SUMMARY")
rep.append("=" * 60)

rep.append(banner("Geometry and flow"))
rep.append(f"  nu                      = {nu:.3e} m2/s")
rep.append(f"  d1 (small pipe)         = {d1:.4f} m      r1 = {r1:.4f} m")
rep.append(f"  d2 (large pipe)         = {d2:.4f} m      r2 = {r2:.4f} m")
rep.append(f"  area ratio A1/A2        = {AR:.4f}        step height h = {h_step:.4f} m")
rep.append(f"  U1                      = {U1:.4f} m/s    Re1 = {Re1:.0f}")
rep.append(f"  U2                      = {U2:.4f} m/s    Re2 = {Re2:.0f}")

rep.append(banner("Entrance lengths -> section lengths"))
rep.append(f"  Le1 = 4.4 Re1^(1/6) d1  = {Le1:.3f} m   ({Le1/d1:.1f} d1)")
rep.append(f"  Le2 = 4.4 Re2^(1/6) d2  = {Le2:.3f} m   ({Le2/d2:.1f} d2)")
rep.append(f"  reattachment ~ 10 h     = {x_reattach:.3f} m")
rep.append(f"  L1 (inlet pipe)         = {L1:.3f} m   ({L1/d1:.1f} d1)")
rep.append(f"  L2 (expansion pipe)     = {L2:.3f} m   ({L2/d2:.1f} d2)")
rep.append(f"  L3 (outlet pipe)        = {L3:.3f} m   ({L3/d1:.1f} d1)")
rep.append(f"  total length            = {Ltot:.3f} m   expansion at x = {x_exp:.3f} m,"
           f" contraction at x = {x_con:.3f} m")

rep.append(banner("Friction, wall units, first cell"))
rep.append(f"  f1 (Haaland)            = {f1:.5f}   (Colebrook: {colebrook_smooth(Re1):.5f})")
rep.append(f"  f2 (Haaland)            = {f2:.5f}   (Colebrook: {colebrook_smooth(Re2):.5f})")
rep.append(f"  u_tau1                  = {u_tau1:.5f} m/s")
rep.append(f"  u_tau2                  = {u_tau2:.5f} m/s")
rep.append(f"  y(y+=1), small pipe     = {y_c1:.3e} m -> first cell height {dw1:.3e} m")
rep.append(f"  y(y+=1), large pipe     = {y_c2:.3e} m -> first cell height {dw2:.3e} m")
rep.append("  NOTE: with a low-Re model the first cell must stay in the viscous sublayer")
rep.append("        on ALL three meshes. The coarser family targets y+ = 1.50 -> 1.00 -> 0.67,")
rep.append("        which is still valid. With wall functions the same refinement would")
rep.append("        walk out of the 30 < y+ < 300 band, which is why the GCI study forces")
rep.append("        the low-Re / resolved-boundary-layer choice.")

rep.append(banner("Turbulent inlet boundary conditions (fully developed pipe)"))
rep.append(f"  I = 0.16 Re^(-1/8)      = {I_inlet:.4f}  ({100*I_inlet:.2f} %)")
rep.append(f"  l = 0.07 d1             = {l_inlet:.5f} m")
rep.append(f"  k     = 1.5 (U I)^2     = {k_inlet:.6f} m2/s2")
rep.append(f"  epsilon = Cmu^0.75 k^1.5 / l = {eps_inlet:.6f} m2/s3")
rep.append(f"  omega   = eps/(Cmu k)   = {omega_inlet:.4f} 1/s")
rep.append("  (Better: map a precursor fully-developed profile, or use a cyclic precursor.")
rep.append("   The values above are the fallback if you keep the long inlet pipe instead.)")

rep.append(banner("Minor losses (Cengel & Cimbala Table 8-4)"))
rep.append(f"  K_SE = (1 - A1/A2)^2            = {K_SE:.4f}   (based on U1)")
rep.append(f"  K_SC = 0.5 (1 - A2/A1)          = {K_SC:.4f}   (based on U1 downstream)")
rep.append(f"  alternative K_SC = 0.42(1-A2/A1) = {K_SC_alt:.4f}  -> quote as uncertainty band")
rep.append(f"  d(p/rho)/dx, small pipe          = {grad1:.5f} m/s2")
rep.append(f"  d(p/rho)/dx, large pipe          = {grad2:.5f} m/s2")
rep.append(f"  jump at expansion   D(p/rho)     = {dp_exp:+.5f} m2/s2  (pressure RISES)")
rep.append(f"  jump at contraction D(p/rho)     = {dp_con:+.5f} m2/s2  (pressure DROPS)")
rep.append(f"  total D(p/rho) inlet->outlet     = {p_over_rho(0.0) - p_over_rho(Ltot):+.5f} m2/s2")

rep.append(banner("Meshes (systematic refinement r = 1.5)"))
rep.append("  level    cells    n_ax1 n_ax2 n_ax3 n_r  n_ann   y+_wall(small pipe)")
for name, N, cfg in levels:
    yp = (cfg["dr_wall1"] / 2.0) * u_tau1 / nu
    rep.append("  %-8s %-8d %-5d %-5d %-5d %-4d %-6d %.3f"
               % (name, N, cfg["n_ax1"], cfg["n_ax2"], cfg["n_ax3"],
                  cfg["n_r_core"], cfg["n_r_ann"], yp))
rep.append("  -> blockMeshDict.coarse / .medium / .fine written next to this script.")
rep.append("  -> run checkMesh on each; max non-orthogonality and skewness must pass.")

report = "\n".join(rep)
print(report)
with open(os.path.join(OUT, "design_summary.txt"), "w") as fh:
    fh.write(report + "\n")

# Single source of truth for the other scripts. postprocess/plot_pressure.py and
# postprocess/gci.py read this instead of keeping their own copies, so changing a
# value at the top of THIS file updates everything downstream.
import json
params = {
    "half_angle_deg": half_angle_deg,
    "nu": nu, "d1": d1, "d2": d2, "Re1": Re1, "Re2": Re2, "U1": U1, "U2": U2,
    "L1": L1, "L2": L2, "L3": L3, "Ltot": Ltot, "x_exp": x_exp, "x_con": x_con,
    "f1": f1, "f2": f2, "u_tau1": u_tau1, "u_tau2": u_tau2,
    "I_inlet": I_inlet, "k_inlet": k_inlet, "eps_inlet": eps_inlet,
    "omega_inlet": omega_inlet,
    "K_SE": K_SE, "K_SC": K_SC, "dp_exp": dp_exp, "dp_con": dp_con,
    "cells": {name: N for name, N, _ in levels},
    "y_plus": {name: (cfg["dr_wall1"] / 2.0) * u_tau1 / nu
               for name, _, cfg in levels},
}
with open(os.path.join(OUT, "case_params.json"), "w") as fh:
    json.dump(params, fh, indent=2)
print("\nWrote case_params.json (read by postprocess/)")

# The inlet values also have to reach the solver. Write them as an OpenFOAM
# include file that base_case/0/{U,k,epsilon} read, so changing a number at the
# top of THIS file and re-running is all it takes - no hand-editing of the case.
INC = os.path.join(OUT, "..", "base_case", "0", "include")
os.makedirs(INC, exist_ok=True)
with open(os.path.join(INC, "initialConditions"), "w") as fh:
    fh.write(f"""// -*- C++ -*-
// GENERATED BY case_design/case_design.py - DO NOT EDIT BY HAND
// Re-run that script after changing any primary choice at the top of it.
//
//   Re1 = {Re1:.0f},  d1 = {d1} m,  nu = {nu:.1e} m2/s
//   I   = {I_inlet:.4f} = 0.16 Re^(-1/8)
//   l   = {l_inlet:.5f} m = 0.07 d1

Uinlet          ({U1:.6f} 0 0);
kInlet          {k_inlet:.6f};
epsilonInlet    {eps_inlet:.6f};
omegaInlet      {omega_inlet:.4f};

// ************************************************************************* //
""")
print("Wrote base_case/0/include/initialConditions (read by 0/U, 0/k, 0/epsilon)")


# ----------------------------------------------------------------------------
# 10. ANALYTICAL PRESSURE CURVE -> CSV (+ plot if matplotlib is available)
# ----------------------------------------------------------------------------
xs_plot, ps_plot = [], []
n_pts = 2000
for i in range(n_pts + 1):
    x = Ltot * i / n_pts
    xs_plot.append(x)
    ps_plot.append(p_over_rho(x))

with open(os.path.join(OUT, "analytical_pressure.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["x_m", "p_over_rho_m2_s2"])
    for x, p in zip(xs_plot, ps_plot):
        w.writerow(["%.6f" % x, "%.6f" % p])

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(xs_plot, ps_plot, lw=2, color="#2166ac",
            label="Bernoulli with losses (analytical)")
    ax.axvline(x_exp, color="#4d9221", ls="--", lw=1, label="sudden expansion")
    ax.axvline(x_con, color="#b2182b", ls="--", lw=1, label="sudden contraction")
    ax.set_xlabel("x  [m]")
    ax.set_ylabel(r"$p/\rho$  [m$^2$/s$^2$]   (p = 0 at outlet)")
    ax.set_title("Group 22 - analytical target for the axial pressure distribution")
    ax.grid(alpha=0.3)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "analytical_pressure.png"), dpi=160)
    print("\nWrote analytical_pressure.png")
except Exception as exc:          # matplotlib not installed -> CSV is still there
    print("\n(plot skipped: %s)" % exc)
