# The Hodgkin-Huxley (HH) model describes the interplay between the sodium and
# potassium flow in and out of the cell to create a spike event. The main
# ingredient is the equations that describe the time-dependence of the sodium
# and potassium channels
import numpy as np

# from HH model: the potassium and sodium channel parameters
EL = -54.4  # mV
EL_K = -77  # mV
EL_Na = 50  # mV
sigma_K = 25
sigma_Na = 5
V_offset_K = -10
V_offset_Na = -30
g_max_K = 0.05  # scaling factor for sigmoidal G_K
g_max_Na = 0.08 # scaling factor for sigmoidal G_Na
g_L = 0.3   # mS /cm^2, specific leak conductance
g_Na = 120  # mS /cm^2, specific sodium conductance
g_K = 36    # mS /cm^2, specific potassium conductance
c_m = 1.    # uF / cm^2, specific membrane capacitance for most neurons

# initial voltage and gating variables
V_0 = -70
n_0 = 0.3   # potassium activation probability
m_0 = 0.05  # sodium activation probability
h_0 = 0.6   # sodium inactivation

def current_pulse(time_Ion, time_Ioff, max_time, dt, I_e):
    """
    Models the current pulse as a delta function.
    """
    time_vector = np.arange(0, max_time, dt)
    V = np.zeros(len(time_vector))
    I_applied = np.zeros(len(time_vector))

    pulse_width = time_Ioff - time_Ion
    tc = (time_Ioff + time_Ion)/2
    taus = tc / 20
    ts = np.arange(time_Ion, time_Ioff, dt)
    pulse = np.exp(-(ts - tc) / taus)*(1 - np.exp(-(ts - t_on) / taus))
    pulse = pulse/np.max(pulse)
    I_applied[round(time_Ion/dt):round(time_Ioff/dt)] = I_e * pulse 

    return I_applied

def sigmoid(V, V_offset, sigma, g_max = 0.01):
    """
    Models neural systems as sigmoidal functions because they are bounded and
    positive.
    """
    return g_max/(1 + np.exp(-(V - V_offset)/sigma))

def simulate_HH(max_time, dt, time_Ion, time_Ioff, I_e):

    time_vector = np.arange(0, max_time, dt)
    V = np.zeros(len(time_vector))
    I_applied = current_pulse(time_Ion, time_Ioff, max_time, dt, I_e)

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
        I_L[t] = g_L * (V[t] - E_L)
        I_Na[t] = g_Na * m[t]**3 * h[t] * (V[t] - E_Na)
        I_K[t] = g_K * n[t]**4 * (V[t] - E_K)
        I_total = -(I_L[t] + I_Na[t] + I_K[t]) + I_applied[t]

        if (t < len(time_vector - 1)):
            # update the gating variables
            m[t+1] = m[t] + (m_inf[t] - m[t]) * dt / tau_m
            h[t+1] = h[t] + (h_inf[t] - h[t]) * dt / tau_h
            n[t+1] = n[t] + (n_inf[t] - n[t]) * dt / tau_n
            V[t+1] = V[t] + I_total[t] / c_m * dt

    return time_vector, V

if __name__ == "__main__":

    dt = 0.01   # ms
    max_time = 20   # ms, max time for simulation
    time_vector = np.arange(0, max_time, dt)
    time_Ion = 5     # ms
    time_Ioff = 6    # ms
    I_e = 20    # uA/cm^2, external current density during pulse

    v_steps = 200
    V = np.linspace(-100., 80., v_steps)
    G_K = sigmoid(V, V_offset_K, sigma_K, g_K)
    I_K = G_K * (V - EL_K)
    G_Na = sigmoid(V, V_offset_Na, sigma_Na, g_Na)
    I_Na = G_Na * (V - EL_Na)


