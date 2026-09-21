"""Outlet-anchored Bernoulli pressure with Darcy friction and step losses."""

import numpy as np


def pressure_profile(params, alpha=1.05, roughness=0.0,
                     contraction_prefactor=0.5, outlet_pressure=0.0, points=300):
    """Return x and kinematic static pressure, retaining both sides of each step.

    Geometry and fluid inputs come from case_params.json. Model defaults follow
    the supplied analytical calculation; outlet_pressure is p/rho in m2/s2.
    """
    boundaries = np.array([0.0, params['x_exp'], params['x_con'], params['Ltot']])
    diameters = np.array([params['d1'], params['d2'], params['d1']])
    velocities = params['U1'] * (diameters[0] / diameters)**2
    reynolds = velocities * diameters / params['nu']
    friction = (-1.8 * np.log10((roughness / diameters / 3.7)**1.11
                               + 6.9 / reynolds))**-2
    gradients = -friction * velocities**2 / (2 * diameters)
    k_exp = (1 - (diameters[0] / diameters[1])**2)**2
    k_cont = contraction_prefactor * (1 - (diameters[2] / diameters[1])**2)
    # Jump = downstream minus upstream pressure; alpha corrects kinetic energy.
    jumps = alpha * (velocities[:-1]**2 - velocities[1:]**2) / 2
    jumps -= np.array([k_exp * velocities[0]**2, k_cont * velocities[2]**2]) / 2
    starts = np.empty(3)
    ends = np.empty(3)
    ends[2] = outlet_pressure
    for i in range(2, -1, -1):
        starts[i] = ends[i] - gradients[i] * (boundaries[i+1] - boundaries[i])
        if i:
            ends[i-1] = starts[i] - jumps[i-1]
    xs = [np.linspace(boundaries[i], boundaries[i+1], points) for i in range(3)]
    ps = [starts[i] + gradients[i] * (xs[i] - boundaries[i]) for i in range(3)]
    return np.concatenate(xs), np.concatenate(ps)
