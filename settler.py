from matplotlib import pyplot as plt

from parameters import *
import parameters as p

V_tot = 500 / 1e6  # Total volume [m^3/s]
L_A = 4  # Settler length [m]
D_A = 0.8  # Settler diameter [m]
D_in = 0.1  # Inlet diameter [m]
ε_in = 0.1  # Intake holdup [-]

Φ_32_0 = 1 / 1000  # Sauter mean diameter [m]
r_s = 0.0030250976562500023

# Numerical values
N_h = 50
N_l = 1000

Δl = L_A / N_l
ε_p = 0.9


def get_η_dis():
    α = p.η_c / (p.η_d + p.η_v)
    Ω = 4 * ε_p ** (7 / 3) + 10 - (84 / 11) * ε_p ** (2 / 3) + 4 * α * (1 - ε_p ** (7 / 3))
    Ω /= 10 * (1 - ε_p ** (10 / 3)) - 25 * ε_p * (1 - ε_p ** (4 / 3)) + 10 * α * (1 - ε_p) * (1 - ε_p ** (7 / 3))
    η_dis = p.η_c * (1 + 5.5 * Ω * ε_p)
    print(f"η_dis = {η_dis}")
    return η_dis


η_dis = get_η_dis()


def get_slip():
    σ_0 = 1
    slip = 1 - np.exp(-12300 * p.ρ_c / p.Δρ * (p.σ / σ_0) ** 3)
    print(f"slip = {slip}")
    # slip = 0.15
    return slip


slip = get_slip()


def L_in(H_p_0):
    v_in = (4 * V_tot) / (π * D_in ** 2)
    v_A = (4 * V_tot) / (π * D_A ** 2)
    ρ = ε_in * p.ρ_d + (1 - ε_in) * p.ρ_c

    Ar_Dis = (g * Φ_32_0 ** 3 * p.Δρ * ρ) / (η_dis ** 2)
    Re_A = (ρ * v_A * D_A) / η_dis
    Re_in = (ρ * v_in * D_in) / η_dis

    Ψ_Φ = Φ_32_0 / (Φ_32_0 + H_p_0)
    Ψ_Re = (Re_in * Re_A) / Ar_Dis
    Ψ_ε = 1 / (1 - ε_in)
    Ψ_ρ = Δρ / ρ
    Ψ_D = Φ_32_0 / D_A

    L_in = 43.7 * Φ_32_0 * (Ψ_Φ ** 0.4) * (Ψ_Re ** 0.5) * (Ψ_ε ** 0.2) * (Ψ_ρ ** 0.2) * (Ψ_D ** 0.1)
    return L_in


# Initial values
H_p_0 = D_A / 2
ΔH_p_0 = H_p_0 / 2
h_p = np.full(N_l, 0.0)

print("Calculating separator...")

for k in range(20):
    l = 0
    i_low = 1
    V_dis = V_tot * ε_in / ε_p
    Δh = H_p_0 / N_h
    h_p[0] = H_p_0
    L_eff = L_A - L_in(H_p_0)
    Φ_32 = np.full(N_h, Φ_32_0)
    j = 0
    while True:
        Δt = h_p[j] * D_A * Δl / V_dis
        j += 1
        l += Δl
        h_p_l = (N_h - i_low + 0.5) * Δh
        for i in range(i_low, N_h):
            h_py = (N_h - i + 0.5) * Δh
            τ_dd = τ(h_py, Φ_32[i], r_s, 'd')
            Φ_32[i] += Δt * Φ_32[i] / (6 * τ_dd)
            if i == i_low:
                if h_p_l < Φ_32[i_low] / 2:
                    h_p_l = Φ_32[i_low] / 2
                τ_di = τ(h_p_l, Φ_32[i_low], r_s, 'i')
                ΔV_dis = Δl * 2 * D_A * ε_di * Φ_32[i_low] / (3 * τ_di * ε_p)
        Δh_p = (126 * (p.η_c + p.η_d) + 11.3 * slip * η_dis) * V_dis * Δl
        Δh_p /= (h_p[j - 1] * D_A ** 3 * g * ε_p * (1 - ε_p) * p.Δρ)
        h_p[j] = h_p[j - 1] - Δh_p
        V_dis -= ΔV_dis
        i_low = round(N_h - (N_h - 1) * V_dis * ε_p / (V_tot * ε_in))
        Δh = h_p[j] * V_tot * ε_in / (N_h * V_dis * ε_p)
        if Δh > h_p[j]:
            Δh = h_p[j]
        if V_dis <= 0 or h_p[j] <= 0 or l > L_eff:
            break
    if V_dis >= 0:
        H_p_0 += ΔH_p_0
    else:
        H_p_0 -= ΔH_p_0
    ΔH_p_0 = ΔH_p_0 / 2

    # Plotting
    plt.clf()
    x_values = np.linspace(L_in(H_p_0), L_A, N_l)
    x_values = np.insert(x_values, 0, 0)

    line_h_p = np.insert(h_p, 0, h_p[0]) * 1000

    # line_Φ_32 = np.insert(Φ_32, 0, Φ_32[0]) * 1000

    plt.plot(x_values, line_h_p, label='$h_p$')
    # plt.plot(x_values, line_Φ_32, label='$\Phi_{32}$')

    plt.legend()
    plt.xlabel("Length [m]")
    plt.ylabel("Packed layer [mm]")

    plt.draw()
    plt.pause(0.0001)

plt.show()
