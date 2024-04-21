from scipy import signal
import numpy as np
import  matplotlib.pyplot as plt
import copy

def power_spectral_density(plot=True):

    """
    Calculate the PSD of white noise

    i.e., a time series of numbers generated from a
    Guassian distribution with zero mean and variance equal to 1.
    """
    trials = [1000, 5000, 10000, 100000]
    if plot:
        fig, ax = plt.subplots(1, len(trials), figsize=(20,4))

    for count, num_trials in enumerate(trials):
        # FFT the Gaussian noise
        f_nu = np.fft.fft(np.random.randn(num_trials))
        power_spectral_density = np.abs(f_nu) ** 2 / num_trials
        print(f"PSD mean: {np.mean(power_spectral_density)}")
        print(f"PSD standard deviation: {np.std(power_spectral_density)}")

        if plot:
            ax[count].plot(power_spectral_density, 'b', ls=':')
            ax[count].set_title('size # trials = %d'%num_trials)
            ax[count].axhline(
                y = np.mean(power_spectral_density),
                color = 'r',
                ls ='-.',
                lw = 4,
                label='mean'
            )
            ax[0].legend()
            ax[0].set_ylabel(r'$\frac{1}{num_trials}$ psd $(\nu)$', fontsize=12)
    if plot:
        fig.tight_layout()
        plt.show()


def clean_vs_noisy_signal(
    max_time=10,
    num_time_steps=10000,
    f1=1,
    f2=2,
    plot=True
):
    dt = max_time / num_time_steps
    time_vector = np.arange(0, max_time, dt)

    clean_signal = np.cos(2 * np.pi * f1 * time_vector) + \
        0.5 * np.cos(2 * np.pi * f2 * time_vector)
    noisy_signal = clean_signal + 1.5 * np.random.randn(num_time_steps)

    ## FFT of noisy signal
    fft_signal = np.fft.fft(noisy_signal)
    real_fft_shift = np.fft.fftshift(fft_signal.real)
    imaginary_fft_shift = np.fft.fftshift(fft_signal.imag)

    ## shifted frequencies
    freq = np.fft.fftfreq(num_time_steps, dt)
    freq_shift = np.fft.fftshift(freq)

    ## power spectrum
    nyquist_freq = 1 / (2 * dt)
    positive_freq = np.arange(0, int(num_time_steps / 2)) / max_time

    if plot:
        fig1, ax = plt.subplots(1, 2, figsize=(12, 4))
        ax[0].plot(time_vector, clean_signal, "b")
        ax[1].plot(time_vector, noisy_signal, "b")
        ax[0].set_xlabel("time (sec)")
        ax[1].set_xlabel("time (sec)")
        ax[0].set_ylabel('clean signal', fontsize=10)
        ax[1].set_ylabel('noisy signal', fontsize=10)
        ax[0].grid(True)
        plt.show()

        fig2, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4))

        ### plot real and imaginary part of the noisy signal FT
        ax1.plot(
            np.arange(num_time_steps),
            (dt / max_time) * fft_signal.real,
            'k',
            label='real part'
        )
        ax1.plot(
            np.arange(num_time_steps),
            (dt / max_time) * fft_signal.imag,
            'r',
            label='imaginary part'
        )
        ax1.legend()
        ax1.set_xlabel('index', fontsize=10)
        ax1.set_ylabel('fft signal', fontsize=10)
        ax1.grid(True)

        ### plot shifted frequency
        ax2.plot(
            freq_shift,
            (dt / max_time) * real_fft_shift,
            'k',
            label='real part'
        )
        ax2.plot(
            freq_shift,
            (dt / max_time) * imaginary_fft_shift,
            'r',
            label='imaginary part'
        )
        ax2.legend()
        ax2.set_xlabel('shifted frequency', fontsize=10)
        ax2.set_ylabel('fft signal', fontsize=10)
        ax2.set_xlim(-10, 10)
        ax2.grid(True)

        ### plot the power spectrum
        ax3.plot(
            positive_freq,
            (dt / max_time) ** 2 * np.abs(fft_signal[:int(num_time_steps / 2)]) ** 2,
            'b'
        )
        ax3.set_xlabel('f (Hz)')
        ax3.set_ylabel('PSD')
        ax3.set_xlim(0, 10)
        fig2.tight_layout(pad=2)
        plt.show()


def correlation(num_time_steps=100000, plot=True):

    stimulus = np.random.randn(num_time_steps)
    autocorr_fftconv = signal.fftconvolve(
            stimulus,
            stimulus[::-1],
            mode='full'
    )
    autocorr_fft = signal.correlate(stimulus, stimulus, mode='full')
    power_spectral_density = np.abs(stimulus) ** 2 / len(stimulus)

    dt = 1
    freq = np.fft.fftfreq(num_time_steps, dt)
    freq_shift = np.fft.fftshift(freq)
    fft_of_autocorr = np.fft.fft(autocorr_fft)
    real_fft_of_autocorr = np.fft.fftshift(fft_of_autocorr.real)

    if plot:
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(12, 3))
        ax1.plot(stimulus)
        ax1.set_title('white noise')
        ax2.plot(
            np.arange(-len(stimulus) + 1, len(stimulus)),
            autocorr_fftconv,
            'b',
            lw = 6,
            label='fftconvolve',
            alpha=0.2
        )
        ax2.plot(
            np.arange(-len(stimulus) + 1, len(stimulus)),
            autocorr_fft,
            'k',
            label='signal.correlate'
        )
        ax2.set_title('autocorrelation')
        ax2.legend()
        ax3.plot(
            power_spectral_density,
            'k',
            label='PSD'
        )
        ax3.legend()
        ax4.plot(
            # freq_shift,
            (dt / num_time_steps) * real_fft_of_autocorr,
            label='FFT of autocorrelation'
        )
        ax4.legend()
        fig.tight_layout()
        plt.show()


def inverse_fourier_transform(
    max_time=10,
    num_time_steps=10000,
    f1=1,
    f2=2,
    plot=True
):
    dt = max_time / num_time_steps
    time_vector = np.arange(0, max_time, dt)

    clean_signal = np.cos(2 * np.pi * f1 * time_vector) + \
        0.5 * np.cos(2 * np.pi * f2 * time_vector)
    noisy_signal = clean_signal + 1.5 * np.random.randn(num_time_steps)

    fft_signal = np.fft.fft(noisy_signal)
    fft_cleaned_signal = copy.deepcopy(np.abs(fft_signal))
    fft_cleaned_signal[np.abs(fft_signal)< 450] = 0
    inv_fft_cleaned_signal = np.fft.ifft(fft_cleaned_signal)

    if plot:
        ### plot clean signal first
        plt.plot(
            np.fft.fftshift(np.fft.fftfreq(num_time_steps, dt)),
            (1 / num_time_steps) * np.roll(
                fft_cleaned_signal,
                int(num_time_steps/2 - 1)), 'b'
        )
        plt.xlabel('frequency (Hz)')
        plt.ylabel('fft')
        plt.xlim(-10, 10)
        plt.show()
        ### plot inverse FFT of the cleaned signal
        fig, ax = plt.subplots(figsize=(15, 4))
        ax.plot(
            time_vector[1:],
            inv_fft_cleaned_signal[1:],
            color ='b',
            label='cleaned data by FFT',
            zorder=2
        )
        ax.plot(
            time_vector,
            clean_signal,
            color = 'red',
            lw= 4,
            alpha=0.5,
            label='original clean signal',
            zorder =1
        )
        ax.plot(
            time_vector,
            noisy_signal,
            color = 'gray',
            lw= 1,
            alpha=0.3,
            label='noisy signal'
        )
        ax.legend()
        plt.show()


def band_pass_filter(
    num_samples=2**13,
    freq1=20,
    freq2=10,
    freq3=5,
    sigma=0.005,
    plot=True
):
    """
    Use FFT to suppress frequencies in the middle of the
    frequency band of the signal.
    """
    max_time = num_samples / 1000
    dt =  max_time / num_samples
    sampling_freq = 1. / dt

    time_vector = dt * np.arange(-num_samples / 2, num_samples / 2)
    freq_vector = sampling_freq / num_samples * \
            np.arange(-num_samples / 2, num_samples / 2)

    signal = np.cos(2 * np.pi * freq1 * time_vector) + \
        0.5 * np.cos(2 * np.pi * freq2 * time_vector) + \
        0.25 * np.cos(2 * np.pi * freq3 * time_vector)
    # shift zero point from center to first point in the array
    signal_shifted = np.roll(signal, int(num_samples / 2))

    # compute the FFT centered
    fft_signal_shifted = np.fft.fft(signal_shifted) / num_samples
    # shift the spectrum to put zero frequency at the middle of the array
    signal_centered = np.roll(fft_signal_shifted, int(num_samples / 2))
    spectrum = np.abs(signal_centered) ** 2

    if plot:
        ### plot the signal
        _, ax = plt.subplots(figsize=(20, 4))
        ax.plot(time_vector, signal, "b-", lw=3)
        ax.plot(time_vector, signal_shifted, "r-.", alpha=0.5)
        ax.set_xlabel("time (sec)")
        ax.set_xlim(-1, 1)
        ax.grid()
        plt.show()

        ### plot FFT of the signal
        _, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4))
        ax1.plot(freq_vector, np.real(signal_centered), linewidth=1)
        ax1.plot(freq_vector, np.imag(signal_centered), 'r', linewidth=1.5)
        ax1.set_xlabel('Frequency (Hz)', fontsize= 12)
        ax1.set_ylabel('FFT')
        ax1.set_xlim(-95, 95)
        ax1.grid()

        ### plot power spectrum
        ax2.plot(freq_vector, spectrum, linewidth=2)
        ax2.set_xlabel('Frequency (Hz)', fontsize=12)
        ax2.set_ylabel('Power', fontsize=12)
        ax2.set_xlim(-95, 95)
        ax2.grid()

        ### inverse FFT to get the original signal
        signal_restored = num_samples * np.fft.ifft(fft_signal_shifted)
        ax3.plot(time_vector, np.roll(signal_restored, int(num_samples / 2)))
        ax3.set_xlim(-1, 1)
        ax3.grid()
        plt.show()

power_spectral_density()
clean_vs_noisy_signal()
correlation()
inverse_fourier_transform()
band_pass_filter()

