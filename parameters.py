import numpy as np
import pandas as pd

print("Importing data...")
params = pd.read_excel('in/data.xlsx', 'Parameters', index_col=0)

g = 9.81  # Acceleration due to gravity [m/s^2]
H_cd = 1e-20  # Hamaker coefficient [Nm] (set to 1*10^-20 by default)

H_0 = params["H_0"]["Value"]  # Total height [m]
H_d = params["H_d"]["Value"]  # Total height [m]
ε_0 = H_d / H_0  # Holdup at t=0 [-]

ρ_c = params["ρ_c"]["Value"]  # Density of continuous phase [kg/m^3]
ρ_d = params["ρ_d"]["Value"]  # Density of dispers phase [kg/m^3]
Δρ = np.abs(ρ_c - ρ_d)  # Density difference [kg/m^3]

σ = params["σ"]["Value"]  # Surface tension [N/m]

η_c = params["η_c"]["Value"]  # Viscosity of continuous phase [Pas]
η_d = params["η_d"]["Value"]  # Viscosity of dispers phase [Pas]
η_v = params["η_v"]["Value"]  # Correction Viscosity for surface-active components [Pas]

# Numerical parameters
N_t = 500
N_h = 500
ε_di = 1

n_c_low = 3
n_c_high = 11

bias = 0.9

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
