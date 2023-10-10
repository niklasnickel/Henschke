import os

import numpy as np
import pandas as pd

g = 9.81  # Acceleration due to gravity [m/s^2]
H_cd = 1e-20  # Hamaker coefficient [Nm] (set to 1*10^-20 by default)


def update_params(file):
    print("Updating parameters...")
    params = pd.read_excel(os.path.join('in', file), 'Parameters', index_col=0)

    global H_0
    H_0 = params["H_0"]["Value"]  # Total height [m]
    global H_d
    H_d = params["H_d"]["Value"]  # Total height [m]
    global ε_0
    ε_0 = H_d / H_0  # Holdup at t=0 [-]

    global ρ_c
    ρ_c = params["ρ_c"]["Value"]  # Density of continuous phase [kg/m^3]
    global ρ_d
    ρ_d = params["ρ_d"]["Value"]  # Density of dispers phase [kg/m^3]
    global Δρ
    Δρ = np.abs(ρ_c - ρ_d)  # Density difference [kg/m^3]

    global σ
    σ = params["σ"]["Value"]  # Surface tension [N/m]

    global η_c
    η_c = params["η_c"]["Value"]  # Viscosity of continuous phase [Pas]
    global η_d
    η_d = params["η_d"]["Value"]  # Viscosity of dispers phase [Pas]
    global η_v
    η_v = params["η_v"]["Value"]  # Correction Viscosity for surface-active components [Pas]


# Numerical parameters
N_t = 500
N_h = 500
ε_di = 1

n_c_low = 3
n_c_high = 11

bias = 0.5

π = np.pi


def τ(h_p, Φ, r_s, ID):
    La_mod = (g * Δρ / σ) ** 0.6 * Φ * h_p ** 0.2
    R_F = Φ * np.sqrt(1 - 4.7 / (4.7 + La_mod))
    if ID == 'd':
        R_F = 0.3025 * R_F
    else:
        R_F = 0.5240 * R_F
    R_a = 0.5 * Φ * (1 - np.sqrt(1 - 4.7 / (4.7 + La_mod)))
    τ = 7.65 * η_c * R_a ** (7 / 3) / (H_cd ** (1 / 6) * σ ** (5 / 6) * R_F * r_s)
    return τ
