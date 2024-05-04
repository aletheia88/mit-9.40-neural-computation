# lecture 1 in class activities
# implementations to `notes/random_diffusion.md`
# author: Alicia KY. Lu

import numpy as np
import matplotlib.pyplot as plt

def simulate_1D_diffusion(n_particles, total_time=2000, plot=True):

    time_vector = np.arange(0, total_time)

    if n_particles == 1:
        displacement = np.zeros(total_time)

        for step in range(1, total_time):
            delta = np.random.binomial(1, 0.5) * 2 - 1
            displacement[step] = displacement[step - 1] + delta

    elif n_particles > 1:
        ### alternatively ###
        deltas = 2 * np.random.binomial(1, 0.5, (n_particles, total_time-1)) - 1
        displacement = np.concatenate((
            np.zeros((n_particles, 1)),
            np.cumsum(deltas, axis = 1)),
            axis=1)

    if plot and n_particles == 1:
        _, ax = plt.subplots(figsize=(4,4))
        ax.plot(time_vector, displacement)
        ax.set(title='Simulating Particle: Trajectory')
        ax.set(xlabel='Timestep', ylabel='Position')
        plt.show()

    elif plot and n_particles > 1:
        _, ax = plt.subplots(figsize=(4,4))

        for i in range(n_particles):
            ax.plot(time_vector, displacement[i],linewidth=0.5)

        ax.set(title='Trajectory: %d particles'%n_particles)
        plt.show()

simulate_1D_diffusion(100)