# lecture 10 in class activities
# implementations to `notes/synaptic_dynamics.md`
# author: Alicia KY. Lu

from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt


def simulate_single_neuron_dynamics(tau, v0, V_init, t_s_vector, dt, t_start, t_end):
    """
    Solve: V(t + dt) = V(t) - dt/tau V(t) + v_0 * tau*n(t)
    """
    time_vector = np.arange(t_start, t_end, dt)
    current = np.zeros_like(time_vector)
    t_s_index = (t_s_vector/dt).astype(int)
    # initialize current
    current[t_s_index] = v0 * tau

    V = np.zeros_like(time_vector)
    V[0] = V_init
    for n in range(time_vector.size - 1):
        V[n+1] = V[n] -dt/tau*V[n] + current[n]

    return time_vector, V

def plot_single_neuron_dynamic(time_vector, V):

    tau = 5     # in ms
    v0 = 1      # in mV
    V_init = 0
    t_s_vector = np.array([10, 40, 50, 55])     # current pulses occuring time in ms
    dt = 0.01
    t_start = 0     # in ms
    t_end = 80      # in ms
    time_vector, V = simulate_single_neuron_dynamics(tau, v0, V_init, t_s_vector, dt, t_start, t_end)

    plt.plot(time_vector, V, 'b')
    plt.xlabel('time (ms)')
    plt.ylabel('V(t) (mV)')
    plt.title(r'$v_0 \cdot \tau$ = %.2f'%(v0 * tau) + " mV")
    plt.grid()
    plt.show()


def simulate_modified_LIF_neuron(time_duration=1000, dt=0.1):

    ### set cell parameters
    V_reset = -80   # reset potential after spike in mV
    tau_membrane = 10      # membrane time constant in ms
    g_L = 10        # leak conductance in nS
    E_L = -75       # leak reversal potential in mV
    V_init = -65    # in mV
    t_refractory = 2    # in ms

    ### set simulation parameters
    time_vector = np.arange(0, time_duration, dt)

    ### set synapse parameters
    gbar_exc = 2.2    # max post-conductance for excitatory currents in nS
    E_exc = 0   # excitatory reversal potential in mV
    tau_synapse_exc = 2    # excitatory synapse time constant in ms
    gbar_inh = 3.5  # max post-conductance for inhibitory currents in nS
    E_inh = -80     # inhibitory reversal potential in mV
    tau_synapse_inh = 10   # inhibitory synapse time constant in ms

    ### set external current input
    I_injected = 0      # in pA

    V_thresholds = np.array([1e3, -55, 1e3, -55])   # spike threshold in mV
    inh_to_exc_ratios = np.array([0.1, 0.1, 0.35, 0.35]) # ratio of inhibition to excitation firing rates

    def simulate_presynaptic_spikes(random_seed=11):
        # model (exc. & inh.) presynaptic neuronal firings as a Poisson event
        np.random.seed(random_seed)
        exc_presynaptic_spikes = np.random.poisson(
            lam = 10*(dt/time_vector.size),
            size=(1000, time_vector.size)).sum(axis=0)

        inh_presynaptic_spikes = np.random.poisson(
            lam = inh_to_exc_ratio * 10*(dt/time_vector.size),
            size=(200, time_vector.size)).sum(axis=0)

        return exc_presynaptic_spikes, inh_presynaptic_spikes

    def simulate_membrane_potential(
            V_threshold,
            exc_presynaptic_spikes,
            inh_presynaptic_spikes
    ):
        V_membrane = np.zeros(time_vector.size)
        V_membrane[0] = V_init
        g_exc = np.zeros(time_vector.size)
        g_inh = np.zeros(time_vector.size)
        I_injected_vector = I_injected * np.ones(time_vector.size)

        refractory_period_time_index = 0
        spikes_time_indices = []   # time indices at which a spike occurs

        for time_index in range(time_vector.size - 1):

            if refractory_period_time_index > 0:
                V_membrane[time_index] = V_reset
                refractory_period_time_index -= 1

            elif V_membrane[time_index] >= V_threshold:
                spikes_time_indices.append(time_index)
                V_membrane[time_index] = V_reset
                refractory_period_time_index = t_refractory / dt

            ### update the synaptic conductance
            g_exc[time_index + 1] = g_exc[time_index] - (dt / tau_synapse_exc) \
                * g_exc[time_index] \
                + gbar_exc * exc_presynaptic_spikes[time_index + 1]
            g_inh[time_index + 1] = g_inh[time_index] - (dt / tau_synapse_inh) \
                * g_inh[time_index] \
                + gbar_inh * inh_presynaptic_spikes[time_index + 1]

            ### calculate the change of membrane potential
            dV = (dt / tau_membrane) * (-(V_membrane[time_index] - E_L)) \
                - (g_exc[time_index + 1] / g_L) * (V_membrane[time_index] - E_exc) \
                - (g_inh[time_index + 1] / g_L) * (V_membrane[time_index] - E_inh) \
                + (I_injected_vector[time_index] / g_L)
            V_membrane[time_index + 1] = V_membrane[time_index] + dV

        spikes_times = np.array(spikes_time_indices) * dt

        return V_membrane, time_vector, g_exc, g_inh


    for inh_to_exc_ratio, V_threshold in zip(inh_to_exc_ratios, V_thresholds):

        exc_presynaptic_spikes, inh_presynaptic_spikes = simulate_presynaptic_spikes()
        plot_presynaptic_spikes(
            exc_presynaptic_spikes,
            inh_presynaptic_spikes,
            inh_to_exc_ratio,
            V_threshold,
        )
        V_membrane, time_vector, g_exc, g_inh = simulate_membrane_potential(
            V_threshold,
            exc_presynaptic_spikes,
            inh_presynaptic_spikes
        )
        plot_membrane_potential(
            V_membrane,
            time_vector,
            g_exc,
            g_inh,
            V_threshold
        )
        print('Average potential = %.2f'%np.mean(V_membrane))


def plot_presynaptic_spikes(
    exc_presynaptic_spikes,
    inh_presynaptic_spikes,
    inh_to_exc_ratio,
    V_threshold,
):
    ### initialize plots
    fig, ax = plt.subplots(1,2)
    if inh_to_exc_ratio == 0.1 and V_threshold == 1e3:

        ax[0].hist(
            exc_presynaptic_spikes,
            color="b",
            label="excitatory",
            alpha=0.8)
        ax[0].hist(
            inh_presynaptic_spikes,
            color="cyan",
            label="inhibitory",
            alpha=0.4)
        ax[0].legend()

    if inh_to_exc_ratio == 0.35 and V_threshold == 1e3:

        ax[1].hist(
            exc_presynaptic_spikes,
            color="b",
            label="excitatory",
            alpha=0.8)
        ax[1].hist(
            inh_presynaptic_spikes,
            color="cyan",
            label="inhibitory",
            alpha=0.4)
        ax[1].legend()
    plt.show()


def plot_membrane_potential(
        V_membrane,
        time_vector,
        g_exc,
        g_inh,
        V_threshold
):
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    ax.plot(
        time_vector,
        V_membrane,
        lw=1.,
        label="membrane voltage",
        zorder=2)
    ax.axhline(
        V_threshold, 0, 1,
        color='k',
        lw=2.,
        ls='--',
        label='Voltage threshold',
        zorder=1)
    ax.axhline(
        np.mean(V_membrane), 0, 1,
        color='g',
        lw=2.,
        ls='--',
        label='Average membrane potential',
        zorder=1)
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("V_m (mv)")
    ax.legend(loc=[1.02, 0.68])
    ax.set_ylim(-80, 10)

    ax.plot(time_vector, g_exc, color='r', lw =1, label="g_exc")
    ax.plot(time_vector, g_inh, color='b', lw=1, label="g_inh")
    fig.tight_layout(pad=2)


def simulate_synchrony(time_duration=200, dt=0.01):

    time_vector = np.arange(0, time_duration, dt)

    ### set cell parameters
    tau_membrane = 20   # in ms
    V_threshold = -54
    V_reset = -80
    E_L = -70
    RI = 25     # voltage due to external source of current in mV
    V_membrane_1 = np.zeros(time_vector.size)    # cell 1 voltage
    V_membrane_2 = np.zeros(time_vector.size)    # cell 2 voltage

    ### set synapse parameters
    tau_synapse = 10    # synaptic time constant in ms

    E_synapse_exc = 0   # in mV
    E_synapse_inh = -80   # in mV

    spike_times_1 = []   # spike times for neuron 1
    spike_times_2 = []   # spike times for neuron 2

    def compute_g_synapse(plot=False):

        g_synapse = (0.05/tau_synapse) * time_vector * np.exp(1 - time_vector/tau_synapse)

        if plot:
            fig, ax = plt.subplots(figsize=(20, 4))
            ax.plot(time_vector, g_synapse, 'b', lw=1)
            ax.set_xlabel('time (ms)')
            ax.set_ylabel('conductance g_s (t)')
            ax.set_xlim(0, 100)
            ax.grid()
            plt.show()

        return g_synapse

    def compute_V_membrane(g_synapse, g_c=1):

        # g_c is set to 1 for excitatory neurons; -1 for inhibitory neurons
        if g_c > 0:
            E_synapse = E_synapse_exc
        else:
            E_synapse = E_synapse_inh

        V_membrane_1[0] = E_L
        V_membrane_2[0] = V_reset

        g_synapse_1 = np.zeros(time_vector.size)
        g_synapse_2 = np.zeros(time_vector.size)

        for t in tqdm(range(time_vector.size - 1)):

            # update synapse conductances
            for spike_t in range(1, 1 + len(spike_times_2)):
                g_synapse_1[t] = g_synapse_1[t] + g_synapse[t - spike_times_2[spike_t - 1]]

            for spike_t in range(1, 1 + len(spike_times_1)):
                g_synapse_2[t] = g_synapse_2[t] + g_synapse[t - spike_times_1[spike_t - 1]]

            # update cell membrane voltages
            V_membrane_1[t + 1] = V_membrane_1[t] \
                + (dt / tau_membrane) * (-(V_membrane_1[t] - E_L)) \
                - g_synapse_1[t] * (V_membrane_1[t] - E_synapse) + RI

            V_membrane_2[t + 1]  = V_membrane_2[t] \
                + (dt / tau_membrane) * (-(V_membrane_2[t]- E_L)) \
                - g_synapse_2[t] * (V_membrane_2[t] - E_synapse) + RI

            if V_membrane_1[t + 1] > V_threshold:
                spike_times_1.append(t + 1)
                V_membrane_1[t + 1] = V_reset

            if V_membrane_2[t + 1] > V_threshold:
                spike_times_2.append(t + 1)
                V_membrane_2[t + 1] = V_reset

        return V_membrane_1, V_membrane_2, g_synapse_1, g_synapse_2

    g_synapse = compute_g_synapse()
    V_membrane_1, V_membrane_2, g_synapse_1, g_synapse_2 = compute_V_membrane(g_synapse)

    ### plot synaptic conductances
    fig, ax1 = plt.subplots(figsize=(20, 4))
    ax1.plot(time_vector, g_synapse_1, color='b', lw=1, label="g_synapse_1")
    ax1.plot(time_vector, g_synapse_2, color='r', lw=1, alpha = 0.6, label="g_synapse_2")
    ax1.grid()
    plt.legend()
    plt.show()

    ### plot spikes
    V_membrane_1_copy = np.copy(V_membrane_1)
    V_membrane_1_copy[spike_times_1] = 10
    V_membrane_2_copy = np.copy(V_membrane_2)
    V_membrane_2_copy[spike_times_2] = 10

    fig, ax = plt.subplots(figsize=(20, 4))
    ax.plot(
        time_vector,
        V_membrane_1_copy,
        color='b',
        zorder=2,
        lw=1,
        label="V_membrane_1"
    )
    ax.plot(
        time_vector,
        V_membrane_2_copy,
        color='r',
        zorder=1,
        lw=2,
        alpha=0.4,
        label="V_membrane_2"
    )
    ax.grid()
    plt.show()

simulate_synchrony()