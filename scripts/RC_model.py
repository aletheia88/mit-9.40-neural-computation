# lecture 2 in class activities
# implementations to `notes/neuron_circuit_model.md`
# author: Alicia KY. Lu

import numpy as np
import matplotlib.pyplot as plt

E_L = -75   # mV
C_m = 0.1   # nF = 1e-9 F
R = 100     # 1e6 Ohms

def simulate_RC_neuron(total_time=700, dt=1):

    tau = R * C_m
    time_vector = np.arange(total_time)
    V = np.zeros(total_time)
    V_inf = np.zeros(total_time)
    V[0] = E_L
    V_inf[0] = V[0]

    I_e = np.zeros(total_time)
    currents = np.arange(100, total_time, 100)

    I_e[np.logical_and(
        time_vector >= currents[0],
        time_vector <= currents[1])
    ] = 0.25

    I_e[np.logical_and(
        time_vector >= currents[2],
        time_vector <= currents[3])
    ] = 0.5

    I_e[np.logical_and(
        time_vector >= currents[4],
        time_vector <= currents[5])
    ] = 1.0

    for t in range(len(time_vector) - 1):
        V_inf[t + 1] = E_L + R * I_e[t]
        V[t + 1] = V[t] + (dt/tau) * (- V[t] + V_inf[t + 1])

    fig1, ax1 = plt.subplots(2, 1)

    ax1[0].plot(time_vector, I_e)
    ax1[0].set(
        xlabel='Time (ms)',
        ylabel='I_e (nA)',
        title='Injected current'
    )
    ax1[1].plot(
        time_vector,
        V_inf,
        'k',
        label='V_inf'
    )
    ax1[1].plot(time_vector, V, 'r', label='V')
    ax1[1].set(
        xlabel='Time (ms)',
        ylabel='V (mV)',
        title='Membrane Voltage'
    )
    ax1[1].legend()
    plt.tight_layout()
    plt.show()

simulate_RC_neuron()
