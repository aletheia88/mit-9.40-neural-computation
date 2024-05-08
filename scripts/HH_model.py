# lecture 6, 7, 8 in class activities
# implementations to `notes/HH-model.md`
# author: Alicia KY. Lu

# The Hodgkin-Huxley (HH) model describes the interplay between the sodium and
# potassium flow in and out of the cell to create a spike event. The main
# ingredient is the equations that describe the time-dependence of the sodium
# and potassium channels

import numpy as np
import matplotlib.pyplot as plt

# from HH model: the potassium and sodium channel parameters
EL = -54.4  # mV
EL_K = -77  # mV potassium equilibrium potential
EL_Na = 50  # mV sodium equilibrium potential
sigma_K = 25
sigma_Na = 5
V_offset_K = -10
V_offset_Na = -30

g_L = 0.3   # mS /cm^2, specific leak conductance
g_Na = 120  # mS /cm^2, specific sodium conductance
g_K = 36    # mS /cm^2, specific potassium conductance
c_m = 1.    # uF / cm^2, specific membrane capacitance for most neurons

# initial voltage and gating variables
V_0 = -70
n_0 = 0.3   # potassium activation probability
m_0 = 0.05  # sodium activation probability
h_0 = 0.6   # sodium inactivation


def sigmoid(V, V_offset, sigma, g_max = 0.01):
    """
    Models neural systems as sigmoidal functions because they are bounded and
    positive.
    """
    return g_max/(1 + np.exp(-(V - V_offset)/sigma))

def simulate_ion_currents(
    sigma_Na=5,
    sigma_K=25,
    gK=0.05,
    gNa=0.8,
    time_steps=200,
    plot=True,
):
    voltage = np.linspace(-100, 80, time_steps)
    # potassium conductance
    G_K = sigmoid(voltage, V_offset_K, sigma_K)
    K_current =  G_K * (voltage - EL_K)
    # sodium conductance
    G_Na = sigmoid(voltage, V_offset_Na, sigma_Na)
    Na_current = G_Na * (voltage - EL_Na)

    if plot:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(voltage, Na_current, 'r-')
        ax.set_xlabel("voltage (mV)")
        ax.set_ylabel("current (mA)")
        ax.grid()
        ax1 = ax.twinx()
        ax1.plot(voltage, G_Na, 'b-')
        ax1.set_ylabel("conductance (S)")
        plt.title("sodium")
        plt.show()

        _, ax2 = plt.subplots(figsize=(10, 6))
        ax2.plot(voltage, K_current, 'r-')
        ax2.set_xlabel("voltage (mV)")
        ax2.set_ylabel("current (mA)")
        ax2.grid()
        ax3 = ax2.twinx()
        ax3.plot(voltage, G_K, 'b-')
        ax3.set_ylabel("conductance (S)")
        plt.title("potassium")
        plt.show()


def current_pulse(
    time_Ion=5,
    time_Ioff=6,
    max_time=20, # ms
    dt=0.01,
    I_e=20, # (uA/cm^2)
    plot=True
):
    """
    Models the current pulse as a delta function.
    """
    time_vector = np.arange(0, max_time, dt)
    V = np.zeros(len(time_vector))
    I_applied = np.zeros(len(time_vector))

    pulse_width = time_Ioff - time_Ion
    time_middle = (time_Ioff + time_Ion)/2
    tau = time_middle / 20

    time_vector_2 = np.arange(time_Ion, time_Ioff, dt)
    pulse = np.exp(-(time_vector_2 - time_middle) / tau) * \
            (1 - np.exp(-(time_vector_2 - time_Ion) / tau))
    pulse = pulse/np.max(pulse)
    I_applied[round(time_Ion/dt):round(time_Ioff/dt)] = I_e * pulse

    if plot:
        plt.plot(time_vector, I_applied)
        plt.show()

    return I_applied

####################################################
### helper functions for simulating the HH model ###
####################################################

def alpha_m(V):
    """Probability sodium channel activated and opens"""
    if V == -40:
        alpha_m = 1
    else:
        alpha_m = (0.1 * (V + 40))/(1 - np.exp(- 0.1 * (V + 40)))
    return alpha_m

def beta_m(V):
    """Probability sodium channel closes"""
    return 4 * np.exp(-(V + 65)/18)

def alpha_n(V):
    """Probability potassium channel opens"""
    if V == -55:
        alpha_n = 0.01/0.1
    else:
        alpha_n = (0.01 * (V + 55)) / (1 - np.exp(-0.1 * (V + 55)))
    return alpha_n

def beta_n(V):
    """Probability potassium channel closes"""
    return 0.125 * np.exp(-0.0125 * (V + 65))

def alpha_h(V):
    return 0.07 * np.exp(-(V + 65)/20)

def beta_h(V):
    return 1 / (1 + np.exp(-(V + 35) / 10))

def tau(x, y):
    return 1 / (x + y)

def gv_inf(x, y):
    """Instanteneous steady-state"""
    return x / (x + y)

def simulate_HH(
    max_time=20,
    dt=0.01,
    time_Ion=5,
    time_Ioff=6,
    I_e=20,
    plot_voltage=True,
    plot_gates=True,
):
    time_vector = np.arange(0, max_time, dt)
    V = np.zeros(len(time_vector))
    V[0] = V_0
    I_applied = current_pulse(
        time_Ion,
        time_Ioff,
        max_time,
        dt,
        I_e,
        plot=False
    )
    # n: potassium activation gating variable
    n = np.zeros(len(time_vector))
    n_inf = np.zeros_like(V)
    n[0] = n_0
    # m: sodium activation gating variable
    m = np.zeros(len(time_vector))
    m_inf = np.zeros_like(V)
    m[0] = m_0
    # h: sodium inactivation gating variable
    h = np.zeros(len(time_vector))
    h_inf = np.zeros_like(V)
    h[0] = h_0

    I_total = np.zeros(len(time_vector))
    I_Na = np.zeros(len(time_vector))
    I_K = np.zeros(len(time_vector))
    I_L = np.zeros(len(time_vector))

    for t in range(0, len(time_vector)):

        # for each gating variable we find the steady state values (_inf) and
        # the time constants (tau_) for each m, h and n
        V_m = V[t]
        tau_m = tau(alpha_m(V_m), beta_m(V_m))
        m_inf[t] = gv_inf(alpha_m(V_m), beta_m(V_m))

        tau_n = tau(alpha_n(V_m), beta_n(V_m))
        n_inf[t] = gv_inf(alpha_n(V_m), beta_n(V_m))

        tau_h = tau(alpha_h(V_m), beta_h(V_m))
        h_inf[t] = gv_inf(alpha_h(V_m), beta_h(V_m))

        # update currents
        I_L[t] = g_L * (V[t] - EL)
        I_Na[t] = g_Na * m[t]**3 * h[t] * (V[t] - EL_Na)
        I_K[t] = g_K * n[t]**4 * (V[t] - EL_K)

        I_total[t] = -(I_L[t] + I_Na[t] + I_K[t]) + I_applied[t]

        if (t < len(time_vector) - 1):
            # update the gating variables
            m[t+1] = m[t] + (m_inf[t] - m[t]) * dt / tau_m
            h[t+1] = h[t] + (h_inf[t] - h[t]) * dt / tau_h
            n[t+1] = n[t] + (n_inf[t] - n[t]) * dt / tau_n
            V[t+1] = V[t] + I_total[t] * (dt / c_m)

    if plot_voltage:
        _, ax1 = plt.subplots(figsize=(4, 4))
        ax1.plot(time_vector, V, 'k-')
        ax1.set_ylabel("V_m (mV)")
        ax2 = ax1.twinx()
        ax2.plot(time_vector, I_applied, 'r--', label="Iapp")
        ax2.set_ylabel('current (uA/cm^2)')
        ax1.legend()
        plt.show()

    if plot_gates:
        fig3, ax3 = plt.subplots(2,1, sharex = True, figsize=(6,4))

        ax3[0].plot(time_vector, m,'k-', label="m (Na)", linewidth = 3)
        ax3[0].plot(time_vector, h, 'g-', label="h (Na)" , linewidth = 3)
        ax3[0].plot(time_vector, n, 'r-.', label='n (K)', linewidth = 3)
        ax3[0].set_ylabel("Gating variable")
        ax3[0].legend()

        ax3[1].plot(time_vector, V, 'b', linewidth = 3)
        ax3[1].set_ylabel("Membrane potential (mV)")
        ax3[1].set_xlabel("time (ms)")
        plt.show()

def simulate_voltage_clamp_experiments(
    E_Na=50,
    E_K=-77,
    g_Na=1.2,
    g_K=0.36,
    V0=-80,
    t_start=10,
    t_end=50,
    dt=0.01,
    plot=True
):
    num_intervals = round(t_end/dt)
    intervals = np.arange(0, num_intervals)
    time_vector = intervals * dt

    def simulate(V_pulse):

        V = V_pulse * np.ones(num_intervals)
        time_off_indices = intervals[intervals < round(t_start/dt)]
        V[time_off_indices] = V0

        n = np.zeros(num_intervals)
        h = np.zeros(num_intervals)
        m = np.zeros(num_intervals)

        # initialize the gating variables at their steady states
        n[0] = gv_inf(alpha_n(V0), beta_n(V0))
        h[0] = gv_inf(alpha_h(V0), beta_h(V0))
        m[0] = gv_inf(alpha_m(V0), beta_m(V0))

        tau_n = np.zeros(num_intervals)
        tau_h = np.zeros(num_intervals)
        tau_m = np.zeros(num_intervals)

        # initialize the gating variable time constants
        tau_n[0] = tau(alpha_n(V0), beta_n(V0))
        tau_h[0] = tau(alpha_h(V0), beta_h(V0))
        tau_m[0] = tau(alpha_m(V0), beta_m(V0))

        # arrays to store the transition rates of different gatings
        alpha_n_values = np.zeros(num_intervals)
        alpha_m_values = np.zeros(num_intervals)
        alpha_h_values = np.zeros(num_intervals)

        beta_n_values = np.zeros(num_intervals)
        beta_m_values = np.zeros(num_intervals)
        beta_h_values = np.zeros(num_intervals)

        n_inf = np.zeros(num_intervals)
        m_inf = np.zeros(num_intervals)
        h_inf = np.zeros(num_intervals)

        iterations = np.arange(0, num_intervals - 1)
        for i in iterations:
            n[i+1] = n[i] + (-n[i] + n_inf[i]) * dt/tau_n[i]
            m[i+1] = m[i] + (-m[i] + m_inf[i]) * dt/tau_m[i]
            h[i+1] = h[i] + (-h[i] + h_inf[i]) * dt/tau_h[i]

            alpha_n_values[i+1] = alpha_n(V[i+1])
            alpha_m_values[i+1] = alpha_m(V[i+1])
            alpha_h_values[i+1] = alpha_h(V[i+1])

            beta_n_values[i+1] = beta_n(V[i+1])
            beta_m_values[i+1] = beta_m(V[i+1])
            beta_h_values[i+1] = beta_h(V[i+1])

            n_inf[i+1] = gv_inf(alpha_n_values[i+1] , beta_n_values[i+1])
            m_inf[i+1] = gv_inf(alpha_m_values[i+1] , beta_m_values[i+1])
            h_inf[i+1] = gv_inf(alpha_h_values[i+1] , beta_h_values[i+1])

            tau_n[i+1] = tau(alpha_n_values[i+1], beta_n_values[i+1])
            tau_m[i+1] = tau(alpha_m_values[i+1], beta_m_values[i+1])
            tau_h[i+1] = tau(alpha_h_values[i+1], beta_h_values[i+1])

        I_Na = g_Na * m**3 * h * (V - E_Na)
        I_K = g_K * n**3 * (V - E_K)

        return {
            "V": V,
            "I_Na": I_Na,
            "I_K": I_K,
            "m": m,
            "n": n,
            "h": h,
        }

    if plot:
        V_pulses = np.array([-40, 0, 30, 70])
        fig, ax = plt.subplots(2, 3, figsize=(12, 6), sharex=True)

        for V_pulse in V_pulses:

            simulation_outcomes = simulate(V_pulse)
            ax[0,0].plot(
                time_vector,
                simulation_outcomes["I_Na"],
                label='V_step = %.2f'%V_pulse
            )
            ax[0,1].plot(time_vector, simulation_outcomes["V"])
            ax[0,2].plot(time_vector, simulation_outcomes["I_K"])
            ax[0,2].set_xlim([5, 25])

            ax[1,0].plot(time_vector, simulation_outcomes["n"])
            ax[1,1].plot(time_vector, simulation_outcomes["m"])
            ax[1,2].plot(time_vector, simulation_outcomes["h"])

        ax[0,0].set_title(
            r"$I_{Na} \,\, \, (\mu A/mm^2)$",
            color='green',
            rotation='horizontal'
        )
        ax[0,0].legend()
        ax[0,1].set_title(
            "V  (mV)",
            rotation='horizontal',
            color='green'
        )
        ax[0,1].set_xlabel('time (ms)')
        ax[0,2].set_title(
            r"$I_K \,\, \, (\mu A/mm^2)$",
            color='green'
        )
        ax[1,0].set_title("n(t)", color='green')
        ax[1,1].set_title("m(t)", color='green')
        ax[1,2].set_title("h(t)", color ='green')
        ax[1,1].set_xlabel("time (ms)")
        plt.tight_layout(pad=2)
        plt.show()

#simulate_HH()
simulate_voltage_clamp_experiments()