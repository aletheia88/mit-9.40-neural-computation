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
# mV this is the reset potential in the LIF mode (in in HH V_reset < E_L)
V_reset = - 65  # mV
area = L * np.pi * (2 * r) + 2 * np.pi * r**2 # mm2
G_L = g_L * area    # units: uS/mm^2 * mm^2 = uS
R_L = 1/G_L   # units: 1/uS = 1/10^(-6) S = 10^6 Ohms = 1 MOhm
C_m = c_m * area    # nF/mm^2*mm^2 = nF = 10^(-9) Farad
tau_m = R_L * C_m   # 20 ms

print("tau = {:.2f} ms".format(tau_m))
print("resistance = {:.2f} MOhm".format(R_L))

def simulate_basic_LIF(
    dt,
    t_dur,
    I_start,
    I_stop,
    I_amplitude,
    V_threshold,
    V_reset,
    reset=True
):
    """
    Simulate a basic version of LIF: after neuron spikes and reaches
    `V_threshold`, it will be reset to `V_reset`
    """
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
        outputs = simulate_basic_LIF(dt, t_dur, I_start, I_stop, I_amplitude,
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


def simulate_extended_LIF(tau, tau_ref, E_L, E_K, V_threshold_0, V_reset, C_m,
        dt, I_amplitude, potassium_on=True):

    """
     Instead of resetting the voltage each time it reaches `V_threshold`, we
     will increase the conductance of potassium out of the cell right. The
     model now has two equations
        C_m dV/dt = -G_L * (V - E_L) - G_K(t) * (V - E_K) + I_e
        τ_r dG_K/dt = -G_K, G_K -> G_K + ΔG_K
    """
    G_L = C_m/(tau/1000.)    # S
    delta_G = 500e-9    #S

    time_vector = np.arange(0, 300, dt) # simulate from 0 to 300 ms
    time_Ion = 100      # time to begin applied current (onset)
    time_Ioff = 200     # time to end applied current (offset)
    index_Ion = round(time_Ion/dt)   # time-point index of current onset
    index_Ioff = round(time_Ioff/dt)     # time-point index of current offset

    # initialize
    G_K = np.zeros(len(time_vector) + 1)
    spikes = np.zeros(len(time_vector) + 1)
    V = E_L * np.ones(len(time_vector) + 1)
    V_thresholds = V_threshold_0 * np.ones(len(time_vector) + 1)
    I = np.zeros(len(time_vector) + 1)
    I[index_Ion:index_Ioff] = I_amplitude * np.ones(index_Ioff - index_Ion)

    def update_for_spike(t, G_K, V, spikes):
        if V[t+1] > V_thresholds[t+1]:
            spikes[t+1] = 1
            if potassium_on:
                G_K[t+1] += delta_G
            else:
                V[t+1] = V_reset
        return G_K, V, spikes

    # iteratively update over all the steps
    for t in range(0, len(time_vector)):
        if potassium_on:
            G_K[t+1] = G_K[t] * np.exp(-dt/tau_ref)
        else:
            G_K[t+1] = 0

        # update Vinf or Vss
         # - G_L * E_L -> mA -> *1e-3 -> A
         # V -> *1e+3 -> mV
        Vinf = 1e+3 * (I[t+1] + G_L * E_L * 1e-3 + \
                G_K[t+1] * E_K * 1e-3)/(G_L + G_K[t+1])

        # update the effective time constant
        tau_eff = 1e+3 * C_m/(G_L + G_K[t+1])

        # update the membrane potential
        V[t+1] = Vinf + (V[t] - Vinf) * np.exp(-dt/tau_eff)

        # update the membrane potential after a potential spike
        G_K, V, spikes = update_for_spike(t, G_K, V, spikes)

    return G_K, V, spikes


if __name__ == "__main__":
    ### parameters for simulating basic LIF model (once reaching V_threshold
    ### reset V to V_reset
    dt = 0.01
    t_dur = 500
    I_start = 20   # ms
    I_stop = 400    # ms
    I_amplitude = 10.0   # nA
    V_threshold = -50
    V_reset = -65

    outputs = simulate_basic_LIF(dt, t_dur, I_start, I_stop, I_amplitude,
            V_threshold, V_reset)
    I_amplitudes = np.arange(1, 15, 0.1)

    ### parameters for simulating extended LIF model
    tau = 10.    # ms
    tau_ref = 2.0  # ms
    E_L = -70   # mV
    E_K = -80.  # mV
    V_threshold_0 = -50.  # mV
    V_reset = E_K   # note that `V_reset` does not have to equal `E_K`
    C_m = 100e-12    # Farad (F)
    dt = 0.01   # ms
    I_low = 200e-12     # 0.2 nA
    I_high = 400e-12    # 0.4 nA

