import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from parameters import *
from scipy import optimize, stats


class SettlingCurve:

    def reset(self):
        self.h_c_calc = list()
        self.h_d_calc = list()
        self.h_pl_calc = list()
        self.t_calc = list()

    def __init__(self):
        data = pd.read_excel('data.xlsx', 'Data')
        self.h_c_exp = data['h_c']
        self.h_c_calc = list()

        self.h_d_exp = data['h_d']
        self.h_d_calc = list()

        self.h_pl_calc = list()

        self.t_exp = data['t']
        self.t_E_exp = self.t_exp.values[-1]  # Time at which the separating surface is half covered with droplets [s]
        self.t_calc = list()
        self.t_E_calc = self.t_E_exp

        self.v_s = self.get_v_s()  # Sedimentation velocity [m/s]
        self.Φ_32_0 = 0  # Sauter mean diameter [m]
        self.get_Φ_32_0()

    def get_settling_curve(self, r_s, plotting=False):
        self.reset()
        # Set initial conditions
        h_c = H_0
        h_d = 0
        Δh_d = 0
        h_p = 1
        Δh = H_0 / N_h

        t = 0
        Δt = self.t_E_exp / N_t
        t_x = 10 * self.t_E_exp

        ε_p = (ε_0 + ε_di) / 2

        Φ_32 = np.full(N_h + 1, self.Φ_32_0)

        while h_p > 0:
            t = t + Δt
            h_d = h_d + Δh_d

            if t < t_x:
                h_c = h_c + self.v_s * Δt
                # h_p = 0.001
                h_p = ((H_0 - h_c) * ε_0 - (1 - ε_0) * h_d) / (ε_p - ε_0)
                if h_d + h_p >= h_c:
                    t_x = t
                    dh_d = Δh_d / Δt
                    C_1 = ((self.v_s - dh_d) * ε_p ** 2 + ε_p * dh_d) / ((h_d - H_0 * ε_0) * (ε_di - ε_p))
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
            τ_di = τ(h_eff, Φ_32[i_low], r_s, 'i')
            Δh_d = 2 * ε_di * Φ_32[i_low] * Δt / (3 * τ_di)
            for i in range(i_low, i_high - 1):
                h_py = (h_d + h_p * ε_p - i * Δh * ε_0) / ε_p
                τ_dd = τ(h_py, Φ_32[i], r_s, 'd')
                Φ_32[i] = Φ_32[i] + Δt * Φ_32[i] / (6 * τ_dd)

            self.h_pl_calc.append(h_d + h_p)
            self.h_d_calc.append(h_d)
            self.h_c_calc.append(h_c)

        self.t_E_calc = t
        self.t_calc = np.linspace(0, self.t_E_calc, len(self.h_c_calc))
        if plotting: self.plot()
        return self

    def get_error(self, bias=1.0):
        point_error = 0
        N = len(self.t_exp)
        for t, h_d in zip(self.t_exp, self.h_d_exp):
            h_d_calc = np.interp(t, self.t_calc, self.h_d_calc)
            point_error += ((h_d_calc - h_d) / H_0) ** 2
        time_error = ((self.t_E_exp - self.t_calc[-1]) / (2 * self.t_calc[-1])) ** 2
        err = bias * np.sqrt(point_error / N) + (1 - bias) * time_error
        # print(f"Q = {err}")
        return err

    def plot(self):
        plt.clf()

        plt.plot(self.t_calc, self.h_d_calc, label='$h_d$')
        plt.scatter(self.t_exp, self.h_d_exp, marker='+')

        plt.plot(self.t_calc, self.h_c_calc, label='$h_c$')
        plt.scatter(self.t_exp, self.h_c_exp, marker='+')

        plt.plot(self.t_calc, self.h_pl_calc, label='$h_d + h_p$')

        plt.legend()
        plt.xlabel("Time [s]")
        plt.ylabel("Height [m]")

        plt.draw()
        plt.pause(0.0001)

    def get_v_s(self):
        print(f"Calculating sedimentation velocity at t₀")
        t = self.t_exp[n_c_low:n_c_high]
        h_c = self.h_c_exp[n_c_low:n_c_high]
        v_s, intercept, r, p, std_err = stats.linregress(t, h_c)
        print(f"v_s = {v_s * 1000} mm/s, σ = {std_err * 1000}\n")
        return v_s

    def get_Φ_32_0(self):
        print(f"Calculating Sauter mean diameter at t₀")

        K_HR = (1 + η_d / η_c) / (2 / 3 + η_d / η_c)

        def Ar(Φ_32_0):
            return (g * Φ_32_0 ** 3 * Δρ * ρ_c) / (η_c ** 2)

        def get_Φ_32_0_new(Re):
            return (Re * η_c) / (ρ_c * abs(self.v_s))

        def step(Φ_32_0):
            if Ar(Φ_32_0) <= 1:
                Re_rs = ((1 - ε_0) * Ar(Φ_32_0) * K_HR) / (18 * np.exp((2.5 * ε_0) / (1 - 0.61 * ε_0)))
                return Φ_32_0 - get_Φ_32_0_new(Re_rs)
            else:
                Re_inf = 9.72 * ((1 + 0.01 * Ar(Φ_32_0)) ** (4 / 7) - 1)
                c_w = Ar(Φ_32_0) / (6 * Re_inf ** 2) - 3 / (K_HR * Re_inf)
                zq2 = (1 - ε_0) / (2 * ε_0 * K_HR) * np.exp((2.5 * ε_0) / (1 - 0.61 * ε_0))
                q3 = 5 * (ε_0 / (1 - ε_0)) ** 0.45 * K_HR ** (- 3 / 2)

                root = np.sqrt(1 + (c_w * q3 * (1 - ε_0) ** 3) / (54 * zq2 ** 2 * ε_0 ** 2) * Ar(Φ_32_0))
                Re_rs = (3 * zq2 * ε_0) / (c_w * q3 * (1 - ε_0)) * (root - 1)
                return Φ_32_0 - get_Φ_32_0_new(Re_rs)

        # Initial estimate
        self.Φ_32_0 = np.sqrt((18 * η_c * abs(self.v_s)) / (g * Δρ))
        print(f"Initial estimate: Φ_32_0 = {self.Φ_32_0 * 1000} mm, Ar = {Ar(self.Φ_32_0)}")

        helping_factor = 10 # Improve numerical stability and velocity
        Φ_32_0 = optimize.fsolve(step, self.Φ_32_0 * helping_factor)
        self.Φ_32_0 = Φ_32_0[0]

        print(f"Φ_32_0 = {self.Φ_32_0 * 1000} mm, Ar = {Ar(self.Φ_32_0)}\n")
        # return self.Φ_32_0


settling_curve = SettlingCurve()



# def import_data():
#     # excel = pd.ExcelFile('data.xlsx')
#     params = pd.read_excel('data.xlsx', 'Parameters')

def get_r_s(plotting=False):
    print(f"Calculating settling curve.")
    r_s_0 = 0.002  # Initial guess of coalescence parameter [-]
    fun = lambda r_s: settling_curve.get_settling_curve(r_s[0], plotting=plotting).get_error(bias=bias)
    result = optimize.minimize(fun, [r_s_0], method='Nelder-Mead', bounds=((0.0001, 0.5),))
    # print(f"Minimization terminated with Result \n{result}")
    print(f"r_s = {result.x[0]}")
    return result.x[0]


if __name__ == '__main__':
    # import_data()
    plt.show()
    settling_curve.get_settling_curve(get_r_s(plotting=True))
    settling_curve.plot()
    plt.show()
    # Q = settling_curve.get_error()
    # print(Q)
