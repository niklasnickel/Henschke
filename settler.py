import numpy as np
import pandas as pd

print("Importing data...")
params = pd.read_excel('data.xlsx', 'Parameters', index_col=0)

V_ges = 0.1  # Total volume [m^3/s]
L_A = 1  # Settler length [m]
D_A = 0.3  # Settler diameter [m]
D_in = 0.1  # Inlet diameter [m]
ε_0 = 0.5  # Intake holdup [-]

Φ_32_0 = 2 / 1000  # Sauter mean diameter [m]

# Numerical values
N_h = 50
N_l = 100

Δl = L_A / N_l
ε_p = 0.9

ρ_c = params["ρ_c"]["Value"]  # Density of continuous phase [kg/m^3]
ρ_d = params["ρ_d"]["Value"]  # Density of dispers phase [kg/m^3]
Δρ = np.abs(ρ_c - ρ_d)  # Density difference [kg/m^3]

σ = params["σ"]["Value"]  # Surface tension [N/m]

η_c = params["η_c"]["Value"]  # Viscosity of continuous phase [Pas]
η_d = params["η_d"]["Value"]  # Viscosity of dispers phase [Pas]
η_v = params["η_v"]["Value"]  # Viscosity of dispers phase [Pas]

π = np.pi
ε_di = 1
g = 9.814
H_cd = 1e-20  # Hamaker coefficient [Nm] (set to 1*10^-20 by default)


def η_dis():
    α = η_c / (η_d + η_v)
    Ω_1 = 4 * ε_p ** (7 / 3) + 10 - (84 / 11) * ε_p ** 2 / 3 + 4 * α * (1 - ε_p ** (7 / 3))
    Ω_2 = 10 * (1 - ε_p ** (10 / 3)) - 25 * ε_p * (1 - ε_p ** (4 / 3)) + 10 * α * (1 - ε_p) * (1 - ε_p ** (7 / 3))
    Ω = Ω_1 / Ω_2
    return η_c * (1 + 5.5 * Ω * ε_p)


def slip():
    σ_0 = 1
    return 1 - np.exp(-12300 * ρ_c / Δρ * (σ / σ_0) ** 3)


def τ(h_p, Φ, ID, r_s):
    La_mod = (g * Δρ / σ) ** 0.6 * Φ * h_p ** 0.2
    R_F = Φ * np.sqrt(1 - 4.7 / (4.7 + La_mod))
    if ID == 'd':
        R_F = 0.3025 * R_F
    else:
        R_F = 0.5240 * R_F
    R_a = 0.5 * Φ * (1 - np.sqrt(1 - 4.7 / (4.7 + La_mod)))
    τ = 7.65 * η_c * R_a ** (7 / 3) / (H_cd ** (1 / 6) * σ ** (5 / 6) * R_F * r_s)
    return τ


def L_in():
    v_in = (4 * V_ges) / (π * D_in ** 2)
    v_A = (4 * V_ges) / (π * D_A ** 2)
    ρ = ε_0 * ρ_d + (1 - ε_0) * ρ_c

    Ar_Dis = (g * Φ_32_0 ** 3 * Δρ * ρ) / (η_dis() ** 2)
    Re_A = (ρ * v_A * D_A) / η_dis()
    Re_in = (ρ * v_in * D_in) / η_dis()

    Ψ_Φ = Φ_32_0 / (Φ_32_0 + H_p_0)
    Ψ_Re = (Re_in * Re_A) / Ar_Dis
    Ψ_ε = 1 / (1 - ε_0)
    Ψ_ρ = Δρ / ρ
    Ψ_D = Φ_32_0 / D_A

    L_in = 43.7 * Φ_32_0 * (Ψ_Φ ** 0.4) * (Ψ_Re ** 0.5) * (Ψ_ε ** 0.2) * (Ψ_ρ ** 0.2) * (Ψ_D ** 0.1)
    return L_in


# Initial values
H_p_0 = D_A / 2
ΔH_p_0 = H_p_0 / 2
h_p = [N_h]

for k in range(20):
    l = 0
    i_low = 1
    V_dis = V_ges * ε_0 / ε_p
    Δh = H_p_0 / N_h
    h_p[0] = H_p_0
    L_eff = L_A - L_in()
    Φ_32 = np.full(N_h, Φ_32_0)
    j = 0
    while True:
        Δt = h_p[j] * D_A * Δl / V_dis
        j += 1
        l += Δl
        h_p = (N_h - i_low + 0.5) * Δh
        for i in range(i_low, N_h):
            h_py = (N_h - i + 0.5) * Δh
            τ_dd = τ(h_py, Φ_32[i], 'd')
            Φ_32[i] += Δt * Φ_32[i] / (6 * τ_dd)
            if i == i_low:
                if h_p < Φ_32[i_low] / 2:
                    h_p = Φ_32[i_low] / 2
                τ_di = τ(h_p, Φ_32[i], 'i')
                ΔV_dis = Δl * 2 * D_A * ε_di * Φ_32[i_low] / (6 * τ_dd)
        Δh_p = (126 * (η_c + η_d) + 11.3 * slip() * η_dis()) * V_dis * Δl
        Δh_p = Δh_p / (Δh_p[j - 1] * D_A ** 3 * g * ε_p * (1 - ε_p) * Δρ)
        h_p[j] = h_p[j - 1] - Δh_p
        V_dis -= ΔV_dis
        i_low = round(N_h - (N_h - 1) * V_dis * ε_p / (V_ges * ε_0))
        Δh = h_p[j] * V_ges * ε_0 / (N_h * V_dis * ε_p)
        if Δh > h_p[j]:
            Δh = h_p[j]
        if V_dis <= 0 or h_p <= 0 or l > L_eff:
            break
    if V_dis >= 0:
        H_p_0 += ΔH_p_0
    else:
        H_p_0 -= ΔH_p_0
    ΔH_p_0 = ΔH_p_0 / 2
