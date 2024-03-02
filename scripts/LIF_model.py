# The leaky integrate and fire (LIF) model is based on the simple RC model of a
# cell, with an extra rule: the model neuron spikes whenever the membrane
# potential passes above a threshold.
# However, real neurons produce spikes due to the presence of voltage-dependent
# ion channels; this dynamics is captured by the Hodgkin-Huxley model
import numpy as np
import matplotlib.pyplot as plt

r = 0.05  # mm
L = 0.1   # mm

g_L = 5   # uS/mm2 - conductance per unit area of the membrane
c_m = 100    # nF/mm2 - capciatnce per unit area of the membrane

E_L = - 70   # mV  equilibium potential for potassium-like ions
# mV this is the reset potential in the LIF mode (in in HH V_reset < EL)
V_reset = - 65  # mV
area = L * np.pi * (2 * r) + 2 * np.pi * r**2 # mm2
G_L = g_L * area    # units: uS/mm^2 * mm^2 = uS
R_L = 1/G_L   # units: 1/uS = 1/10^(-6) S = 10^6 Ohms = 1 MOhm
C_m = c_m * area    # nF/mm^2*mm^2 = nF = 10^(-9) Farad
tau_m = R_L * C_m   # 20 ms

print("tau = {:.2f} ms".format(tau_m))
print("resistance = {:.2f} MOhm".format(R_L))

def simulate_LIF(
    dt,
    t_dur,
    I_start,
    I_stop,
    I_amplitude,
    V_threshold,
    V_reset,
    reset=True
):
    num_steps = int(np.ceil(t_dur/dt))
    V = np.zeros(num_steps + 1)
    V[0] = E_L
    time_vector = dt * np.arange(0, num_steps + 1)
    I_e = np.zeros(num_steps + 1)
    I_e[np.logical_and(
        time_vector>=I_start, time_vector<=I_stop)] = I_amplitude

    # integrate using the Euler method
    for t in range(0, num_steps):
        V_inf = E_L + R_L * I_e[t]
        V[t+1] = V_inf + (V[t] - V_inf) * np.exp(-dt/tau_m)
        if V[t+1] > V_threshold and reset == True:
            V[t+1] = V_reset

    theoretical_firing_rate = (1000/tau_m) * \
            (0.5 + \
            (E_L - V_threshold + R_L * I_amplitude)/(V_threshold - V_reset))

    arg = (E_L - V_reset + R_L * I_amplitude)/(E_L - V_threshold + R_L * I_amplitude)
    I_c = (V_threshold - E_L)/R_L
    if I_amplitude > I_c:
        theoretical_firing_rate = (1000/tau_m) / np.log(arg)
    else:
        theoretical_firing_rate = 0

    # compute numerical firing rates
    spike_idx = np.where(V == V_reset)[0]
    spike_interval = 0
    if len(spike_idx) > 2:
        spike_interval = dt * (spike_idx[2] - spike_idx[1])
    if spike_interval > 0:
        numerical_firing_rate = 1000/spike_interval
    else:
        numerical_firing_rate = 0

    return {"V": V,
            "I_e": I_e,
            "time": time_vector,
            "f_numerical": numerical_firing_rate,
            "f_theoretical": theoretical_firing_rate
    }


def firing_rates_vs_curent(dt, t_dur, I_amplitudes, V_threshold, V_reset,
        I_start, I_stop):
    """
    compute firing rates as a function of various current amplitudes.
    """
    # 1. compute V(I_e, ...)
    # 2. compute delta_T
    # 3. compute f
    num_steps = int(np.ceil(t_dur/dt))
    I_vector = np.ones(num_steps + 1)
    firing_rates_numerical = np.array([])
    firing_rates_theoretical = np.array([])

    for I_amplitude in I_amplitudes:
        outputs = simulate_LIF(dt, t_dur, I_start, I_stop, I_amplitude,
            V_threshold, V_reset)
        firing_rates_numerical = np.append(
            firing_rates_numerical,
            outputs["f_numerical"]
        )
        firing_rates_theoretical = np.append(
            firing_rates_theoretical,
            outputs["f_theoretical"]
        )

    return firing_rates_numerical, firing_rates_theoretical


if __name__ == "__main__":
    dt = 0.01
    t_dur = 500
    I_start = 20   # ms
    I_stop = 400    # ms
    I_amplitude = 10.0   # nA
    V_threshold = -50
    V_reset = -65

    outputs = simulate_LIF(dt, t_dur, I_start, I_stop, I_amplitude,
            V_threshold, V_reset)
    I_amplitudes = np.arange(1, 15, 0.1)

