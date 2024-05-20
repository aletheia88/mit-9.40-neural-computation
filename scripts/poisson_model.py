# lecture 13 in class activities
# implementations to `notes/poisson_model.md`
# author: Alicia KY. Lu

from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import math

def low_pass_filter(dataset_name="dataFee.mat"):

    raw_data = loadmat(f"data/{dataset_name}")["rawData"]
    N = np.size(raw_data)

    F_sampling = 30000  # sampling frequency in Hz
    dt = 1 / F_sampling
    time_vector = np.arange(0, N * dt, dt)

    filter_width = 0.01 # sec
    num_samples = int(np.round(F_sampling * filter_width))
    low_pass_filter = np.ones(num_samples) / num_samples

    _, ax1 = plt.subplots(figsize=(4, 4))
    ax1.stairs(low_pass_filter, color='b', lw=3)
    ax1.set_title("low-pass filter")
    ax1.set_xlabel("time steps")
    plt.show()

    _, (ax, ax2) = plt.subplots(1, 2, figsize=(4, 4))
    ax.plot(time_vector, raw_data)
    ax.set(
        xlim=(0, 2),
        xlabel="time (sec)",
        ylabel="voltage (μV)",
        title="raw signal"
    )

    low_pass_filtered_data = np.convolve(
        np.squeeze(raw_data),
        low_pass_filter,
        mode="same"
    )
    ax2.plot(time_vector, low_pass_filtered_data)
    ax2.set(
        xlim=(0, 2),
        ylim=(-800, 400),
        xlabel="time (sec)",
        ylabel="voltage (mV)",
        title="low-passed LFP"
    )
    plt.tight_layout()
    plt.show()


def high_pass_filter(dataset_name="dataFee.mat"):

    raw_data = loadmat(f"data/{dataset_name}")["rawData"]
    N = np.size(raw_data)

    F_sampling = 30000
    dt = 1/F_sampling

    time_vector = np.arange(0, N * dt, dt)

    filter_width = 0.01 # sec
    num_samples = int(np.round(F_sampling * filter_width)) # num. samples per filter width
    low_pass_filter = np.ones(num_samples) / num_samples
    high_pass_filter = - low_pass_filter

    # create unit impulse
    high_pass_filter[math.ceil(num_samples/2)] = (num_samples - 1) / num_samples
    _, ax = plt.subplots()
    ax.plot(high_pass_filter, "b", label="high pass filter")
    ax.plot(low_pass_filter, "k", lw=10, alpha=0.3, label="low pass filter")
    ax.legend()
    plt.show()

    high_pass_filtered_data = np.convolve(
        np.squeeze(raw_data),
        high_pass_filter,
        mode="same"
    )
    _, (ax1, ax2) = plt.subplots(1, 2, figsize=(4, 4))

    ax1.plot(time_vector, high_pass_filtered_data)
    ax1.set(
        xlim=(0, 2),
        ylim=(-800, 800),
        xlabel="time (sec)",
        ylabel="voltage (mV)",
        title="high-passed LFP"
    )
    ax2.plot(time_vector, np.squeeze(raw_data))
    ax2.set(
        xlim=(0, 2),
        ylim=(-800, 800),
        xlabel="time (sec)",
        ylabel="voltage (mV)",
        title="raw data"
    )
    plt.tight_layout()
    plt.show()

low_pass_filter()
high_pass_filter()
