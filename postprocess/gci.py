
import math


def representative_h(N, dim=2, extent=1.0):
    """Representative grid size. For a 2D (axisymmetric wedge) study the standard
    choice is h ~ (A/N)^(1/2); for 3D, (V/N)^(1/3)."""
    return (extent / N) ** (1.0 / dim)


def apparent_order(e21, e32, r21, r32, tol=1e-10, itmax=500):
    ratio = e32 / e21
    s = 1.0 if ratio > 0 else -1.0
    p = 2.0
    for _ in range(itmax):
        q = math.log((r21 ** p - s) / (r32 ** p - s))
        p_new = abs(math.log(abs(ratio)) + q) / math.log(r21)
        if abs(p_new - p) < tol:
            return p_new, s
        p = p_new
    return p, s


def gci(f1, f2, f3, N1=None, N2=None, N3=None,
        h1=None, h2=None, h3=None, dim=2, Fs=1.25):
    """
    f1 finest ... f3 coarsest. Give eiated value,
    the relative errors and GCI on the fine and medium meshes.ther cell counts (N) or grid sizes (h).
    Returns a dict with the apparent order, the Richardson-extrapol
    """
    if h1 is None:
        h1 = representative_h(N1, dim)
        h2 = representative_h(N2, dim)
        h3 = representative_h(N3, dim)

    r21, r32 = h2 / h1, h3 / h2
    e21, e32 = f2 - f1, f3 - f2

    if abs(e21) < 1e-14:
        return {"note": "f1 == f2 to machine precision; GCI undefined (already converged)",
                "r21": r21, "r32": r32}

    p, s = apparent_order(e21, e32, r21, r32)

    f_ext21 = (r21 ** p * f1 - f2) / (r21 ** p - 1.0)
    e_a21 = abs((f1 - f2) / f1)
    e_ext21 = abs((f_ext21 - f1) / f_ext21)
    gci21 = Fs * e_a21 / (r21 ** p - 1.0)

    e_a32 = abs((f2 - f3) / f2)
    gci32 = Fs * e_a32 / (r32 ** p - 1.0)

    monotone = (e21 * e32) > 0

    return {
        "h1": h1, "h2": h2, "h3": h3,
        "r21": r21, "r32": r32,
        "e21": e21, "e32": e32,
        "monotone_convergence": monotone,
        "apparent_order_p": p,
        "extrapolated_f": f_ext21,
        "approx_rel_error_fine_pct": 100 * e_a21,
        "extrap_rel_error_fine_pct": 100 * e_ext21,
        "GCI_fine_pct": 100 * gci21,
        "GCI_medium_pct": 100 * gci32,
        "f_fine_with_bar": (f1, 100 * gci21),
    }


def report(name, res):
    print(f"\n{name}")
    print("-" * len(name))
    if "note" in res:
        print("  " + res["note"])
        return
    print(f"  refinement ratios      r21 = {res['r21']:.4f}   r32 = {res['r32']:.4f}")
    print(f"  differences            e21 = {res['e21']:+.6g}  e32 = {res['e32']:+.6g}")
    print(f"  monotone convergence   {res['monotone_convergence']}")
    print(f"  apparent order p       {res['apparent_order_p']:.3f}")
    print(f"  extrapolated value     {res['extrapolated_f']:.6g}")
    print(f"  approx. rel. error     {res['approx_rel_error_fine_pct']:.3f} %")
    print(f"  extrap. rel. error     {res['extrap_rel_error_fine_pct']:.3f} %")
    print(f"  GCI (fine)             {res['GCI_fine_pct']:.3f} %   <- the error bar to plot")
    print(f"  GCI (medium)           {res['GCI_medium_pct']:.3f} %")
    if not res["monotone_convergence"]:
        print("  WARNING: non-monotone convergence. Either the meshes are not in the")
        print("           asymptotic range, or the runs are not converged to the same")
        print("           residual level. Check residuals before trusting this GCI.")


# cell counts from case_design.py
N_FINE, N_MEDIUM, N_COARSE = 81000, 36000, 16033

if __name__ == "__main__":
    print("GCI self-test with the  mesh sizes and PLACEHOLDER monitor values.")
    print("Replace the f-values with the real pressure drops once the runs finish.")
    print(f"cells: fine {N_FINE}, medium {N_MEDIUM}, coarse {N_COARSE}")

    # Placeholder: plausible dp/rho across the contraction, converging from below
    report("Monitor 1: dp/rho across the contraction  [m2/s2]  (PLACEHOLDER)",
           gci(f1=2.601, f2=2.588, f3=2.551,
               N1=N_FINE, N2=N_MEDIUM, N3=N_COARSE, dim=2))

    report("Monitor 2: dp/rho across the expansion  [m2/s2]  (PLACEHOLDER)",
           gci(f1=0.742, f2=0.735, f3=0.713,
               N1=N_FINE, N2=N_MEDIUM, N3=N_COARSE, dim=2))
