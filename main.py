import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


class SettlingCurve:
    h_c_exp = None,
    h_c_calc = None
    h_d_exp = None
    h_d_calc = None
    t_exp = None
    t_calc = None

    def __init__(self):
        data = pd.read_excel('data.xlsx', 'Data')
        self.h_c_exp = data['h_c']
        self.h_c_calc = list()

        self.h_d_exp = data['h_d']
        self.h_d_calc = list()

        self.t_exp = data['t']

    def get_settling_curve(self):
        # Set initial conditions
        h_c = H_0
        h_d = 0
        Δh_d = 0
        h_p = 1
        Δh = H_0 / N_h

        t = 0
        Δt = t_E / N_t
        t_x = 10 * t_E

        ε_p = (ε_0 + ε_di) / 2

        Φ_32 = np.full(N_h + 1, Φ_32_0)

        while h_p > 0:
            t = t + Δt
            h_d = h_d + Δh_d

            if t < t_x:
                h_c = h_c + v_s * Δt
                # h_p = 0.001
                h_p = ((H_0 - h_c) * ε_0 - (1 - ε_0) * h_d) / (ε_p - ε_0)
                if h_d + h_p >= h_c:
                    t_x = t
                    dh_d = Δh_d / Δt
                    C_1 = ((v_s - dh_d) * ε_p ** 2 + ε_p * dh_d) / ((h_d - H_0 * ε_0) * (ε_di - ε_p))
                    C_2 = - C_1 * t_x - np.log(ε_di - ε_p)
            if t >= t_x:
                ε_p = ε_di - np.exp(-C_1 * t - C_2)
                h_p = (H_0 * ε_0 - h_d) / ε_p
                h_c = h_d + h_p
            i_low = round(h_d / (Δh * ε_0))
            i_high = round((h_d + h_p * ε_p) / (Δh * ε_0))
            if h_p < Φ_32[i_low] / 2:
                h_eff = Φ_32[i_low] / 2
            else:
                h_eff = h_p
            τ_di = τ(h_eff, Φ_32[i_low], 'i')
            Δh_d = 2 * ε_di * Φ_32[i_low] * Δt / (3 * τ_di)
            for i in range(i_low, i_high - 1):
                h_py = (h_d + h_p * ε_p - i * Δh * ε_0) / ε_p
                τ_dd = τ(h_py, Φ_32[i], 'd')
                Φ_32[i] = Φ_32[i] + Δt * Φ_32[i] / (6 * τ_dd)

            # Plotting
            packed_layer_curve.append(h_d + h_p)
            self.h_d_calc.append(h_d)
            self.h_c_calc.append(h_c)

        t_E_calc = t
        self.t_calc = np.linspace(0, t_E_calc, len(packed_layer_curve))

    def plot(self):
        plt.subplots()

        plt.plot(self.t_calc, packed_layer_curve, label='h_d + h_p')
        plt.plot(self.t_calc, self.h_d_calc, label='h_d')
        plt.plot(self.t_calc, self.h_c_calc, label='h_c')

        plt.scatter(self.t_exp, self.h_d_exp)
        plt.scatter(self.t_exp, self.h_c_exp)

        plt.legend()
        plt.draw()

        print(f"{self.t_calc[-1]}")
        plt.show()


settling_curve = SettlingCurve()

g = 9.81  # Acceleration due to gravity [m/s^2]
Δρ = 81.2  # Density difference [kg/m^3]
σ = 30 / 1000  # Surface tension [N/m]
η_c = 1 / 1000  # Viscosity of continuous phase [Pas]

H_cd = 1e-20  # Hamaker coefficient [Nm] (set to 1*10^-20 by default)

# Measured values
H_0 = 1  # Total height [m]
t_E = 1600  # Time at which the separating surface is half covered with droplets [s]
ε_0 = 0.37  # Holdup at t=0 [-]

# Numerical parameters
N_t = 500
N_h = 500
ε_di = 1

r_s = 0.002  # Initial guess of coalescence parameter [-]

# Derived values
Φ_32_0 = 2 / 1000  # Sauter mean diameter [m]
v_s = -3 / 1000  # Sedimentation velocity [m/s]

packed_layer_curve = list()


def τ(h_p, Φ, ID):
    La_mod = (g * Δρ / σ) ** 0.6 * Φ * h_p ** 0.2
    R_F = Φ * np.sqrt(1 - 4.7 / (4.7 + La_mod))
    if ID == 'd':
        R_F = 0.3025 * R_F
    else:
        R_F = 0.5240 * R_F
    R_a = 0.5 * Φ * (1 - np.sqrt(1 - 4.7 / (4.7 + La_mod)))
    τ = 7.65 * η_c * R_a ** (7 / 3) / (H_cd ** (1 / 6) * σ ** (5 / 6) * R_F * r_s)
    return τ


# def Q():
#     sum = 0
#     for t, h_d in t_exp, h_d_exp:
#         h_d_calc = np.interp(t, h_d_calc, fp)
#         pass


def import_data():
    # excel = pd.ExcelFile('data.xlsx')
    params = pd.read_excel('data.xlsx', 'Parameters')


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    import_data()
    settling_curve.get_settling_curve()
    settling_curve.plot()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
