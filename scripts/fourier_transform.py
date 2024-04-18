import numpy as np
import  matplotlib.pyplot as plt

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

if __name__ == "__main__":
    power_spectral_density()