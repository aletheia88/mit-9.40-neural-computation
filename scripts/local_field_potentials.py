# lecture 12 in class activities
# implementations to `notes/local_field_potentials.md`
# author: Alicia KY. Lu

from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt

def rectangle_kernel(
    input_size=100,
    kernel_size=20,
    random_seed=9041,
    plot=True
):
    np.random.seed(random_seed)
    # generated random stimulus
    a_n = np.arange(0, input_size) + 4 * np.random.randn(input_size)
    kernel = np.ones(kernel_size) / kernel_size

    # compute the convoluted output
    y_output = np.zeros(input_size + kernel_size - 1)
    time_vector = np.arange(input_size + kernel_size - 1)

    # modify kernel & signal sizes for easier calculations
    kernel_ = np.zeros(input_size + kernel_size - 1)
    kernel_[:kernel_size] = kernel

    stimulus_ = np.zeros(input_size + kernel_size - 1)
    stimulus_[:input_size] = a_n

    for n in range(input_size + kernel_size - 1):
        for m in range(kernel_size + 1):
            y_output[n] += stimulus_[n - m] * kernel_[m]

    if plot:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(a_n, 'b-', label='input signal')
        ax.plot(y_output, 'g', label='full output')
        ax.grid()
        ax.legend()
        ax2 = ax.twinx()
        K_shift = np.zeros(100)
        K_shift[40:60] = kernel
        ax2.plot(K_shift, color='k', label='kernel K(t)')
        ax2.set_ylabel('K(t)')
        ax2.legend()
        plt.show()

def general_stimulus(
    max_time,
    dt=0.05,
    plot=True
):
    tau = 10    # in ms
    time_vector = np.arange(0, max_time, dt)
    stimulus = (np.cos(np.pi * time_vector / 50) + np.sin(np.pi * time_vector / 20))
    kernel = np.exp(-time_vector / tau)
    V_numerical = np.zeros(time_vector.size)
    V = np.zeros(time_vector.size)

    for t in tqdm(range(time_vector.size - 1)):
        V[t + 1] = (V[t] + stimulus[t + 1] * dt) / (1 + dt / tau)

        for j in range(t):
            V_numerical[t + 1] += kernel[t - j] * stimulus[j] * dt

    V_convoluted = np.convolve(kernel, stimulus, mode='full') * dt

    if plot:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(time_vector, V, 'cyan', lw = 3, alpha=0.4)
        ax.plot(time_vector, V_numerical, 'k-.', label='numerical')
        ax.plot(
            time_vector,
            V_convoluted[:time_vector.size],
            color='orange',
            label='convolution'
        )
        plt.legend()
        plt.show()

general_stimulus(500)