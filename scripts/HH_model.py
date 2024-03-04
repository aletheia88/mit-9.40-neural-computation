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

def simulate_HH():
    pass

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


