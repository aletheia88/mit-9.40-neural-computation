from scipy import optimize
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


def binocular_rivalry(
    a_parameters,
    I_parameters,
    dt=None,
    time_max=None,
    plot1=True,
    plot2=True,
    plot3=True
):
    def F(x):
        """Sigmoid function of activation"""
        return 1 / (1 + np.exp(1 - x))

    def dF(x):
        """Derivative of the sigmoid function"""
        return np.exp(1 - x) * F(x)**2

    v = np.arange(0, 1.2, 0.01)
    # fix a vary I
    if len(a_parameters) == 1 and len(I_parameters) > 1:
        a = a_parameters[0]
        if plot1:
            fig, ax = plt.subplots()
            for I in I_parameters:
                # where `I - av` intersects v is the steady-state solution
                ax.plot(v, F(I - a * v), label='%.2f'%I)
                ax.plot(v, v, 'k--')
            ax.legend(title='I_app')
            ax.grid()
            ax.set_xlabel('v')
            plt.show()

    elif len(a_parameters) > 1 and len(I_parameters) == 1:
        I = I_parameters[0]
        if plot1:
            fig, ax = plt.subplots()
            for a in a_parameters:
                # where`I - av` intersects v is the steady-state solution
                ax.plot(v, F(I - a * v), label='%.2f'%a)
                ax.plot(v, v, 'k--')
            ax.legend(title='a')
            ax.grid()
            ax.set_xlabel('v')
            plt.show()

    def vect_func(v, t, I_app, a):
        """2D dynamical systems modeling the binocular rivalry"""
        v1, v2 = v
        return [-v1 + F(I_app - a * v2), -v2 + F(I_app - a * v1)]

    # initial conditions: 30 random paired values of v1 and v2
    X0 = np.random.uniform(0, 1, size=(30, 2))
    v_final =[]
    I_sweep = []
    time_vector = np.arange(0, time_max, dt)

    for x0 in (X0):
        for I_app in I_parameters:
            # initial condition of v = [v1, v2]
            v0 = [x0[0], x0[1]]
            sol = odeint(vect_func, v0, time_vector, args=(I_app, a))
            v1 = sol[:, 0]
            v2 = sol[:, 1]
            v_final.append(v1[-1])
            I_sweep.append(I_app)
    if plot2:
        _, ax = plt.subplots()
        ax.scatter(I_sweep, v_final, c='b', lw=0.5)
        ax.set_xlabel('I_app')
        ax.set_ylabel('v1')
        ax.set_title('a = %.2f'%a)
        ax.grid()
        plt.show()

    if plot3:
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 4))
        x0 = [0.2, 0.5]
        v_final = []
        for lw, I_app in enumerate(I_parameters):
            sol = odeint(vect_func, x0, time_vector, args=(I_app, a))
            v1 = sol[:, 0]
            v2 = sol[:, 1]
            # v1's time evolution
            ax1.plot(time_vector, v1,  lw =1+lw/2,  label='%.2f'%I_app)
            ax1.set_xlabel("time")
            ax1.set_ylabel("v1")
            ax1.legend(bbox_to_anchor=(-0.2, 1),title='I_app')
            v_final.append(v1[-1])
            # v2's time evolution
            ax2.plot(time_vector, v2, lw =1+lw/2)
            ax2.set_xlabel("time")
            ax2.set_ylabel("v2")
            # v1, v2 in relation to each other
            ax3.scatter(v1, v2, lw=0.5, alpha = 0.6)
            ax3.set_xlabel('v1')
            ax3.set_ylabel('v2')

        plt.show()
        fig.suptitle(
            r'$v_1(0), v_2(0) = (%.2f, %.2f),\quad a = %.2f$'%(x0[0], x0[1],
                a))


# binocular_rivalry(
#     a_parameters=[6],
#     I_parameters=[0, 1, 2, 4, 5, 6, 8, 10]
# )
# binocular_rivalry(
#     a_parameters=[0.1, 1, 2, 4, 6, 8, 10],
#     I_parameters=[4]
# )
# binocular_rivalry(
#     a_parameters=[5],
#     I_parameters=list(np.arange(0, 10, 0.5)),
#     dt=0.01,
#     time_max=40,
#     plot1=False,
#     plot2=True,
#     plot3=False,
# )
binocular_rivalry(
    a_parameters=[5],
    I_parameters=list(np.arange(0.5, 8, 1)),
    dt=0.01,
    time_max=40,
    plot1=False,
    plot2=False,
    plot3=True
)
