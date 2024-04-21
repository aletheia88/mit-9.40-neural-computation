## Firing rate models

Communication between single neurons is transmitted throguh firing rates. For a population of neurons, the average firing rate will be the macroscopic variable that we want to compute given a stimulus. A firing rate model provides an output average firing rates of interacting populations of neurons. It is considered a mean-field theory where correlations or fluctuations are neglected. It's often required that any presynaptic input is also uncorrelated. This means for the averaging to make sense, the intervals to average over should be much larger than presynaptic interspike intervals.

Consider $n$ interacting populations of neurons with averaging firing rates $r_{i=1,2,3,\cdots,n}$. Let $s_i$ be the fraction of active presynaptic channels to population $i$ as a result of incoming presynaptic firing. The dynamics is described as a coupled system:
$$\tau_i \frac{dr_i}{dt} = -r_i + f_i \left( \sum_{j=1}^n W_{ji} s_j \right)$$
$$\tau_{s_i} \frac{ds_i}{dt} = -s_i + F(r_i)$$
where $W_{ji}$ equal to the strength of connection between neurons in population $i$ and $j$. The weights $W_{ji}$ are positive for excitatory connections and negative for inhibitory connections.

The common choices of a firing rate function $f_i$ are:
$$\text{the linear-rectified unit (ReLu):} \quad f(I) = r_0 [I - I_{th}]_{+}^{\alpha}$$
$$\text{the sigmoid:} \quad \quad f(I) = \frac{r_0}{1 + \exp(-\frac{(I - I_0)}{\sigma})}$$

### time scale of synapses vs. firing rates
Changes in synaptic inputs occur on the order of 200 to 500 micro-seconds scale. The neuronal firings are of the order of 10-100 milliseconds. Therefore, $\tau_s \ll \tau_r$, so we replace the time-dependent synaptic dynamics by its steady-state. $$-s_i + F(r_i) = 0$$
where $F(r_i)$ is the synaptic gate function that is taken to be
$$F(r_i) = \frac{r}{r_{max}}$$

So, we will only simulate the following equation and solve the firing rate as a function of time:
$$\tau_i \frac{dr_i}{dt} = -r_i + f_i \left( \sum_{j=1}^n \frac{1}{r_{j}^{max}} W_{ji} r_j \right)$$

### binocular rivalry
The following model has been used to describe binocular rivalry,
$$\frac{d v_1}{d t} = - v_1 + F(I_{app} - a v_2 -  b u_1)$$
$$\frac{d v_2}{d t} = - v_2 + F(I_{app} - a v_1 -  b u_2)$$
$$\frac{d u_1}{d t} = - \frac{1}{\tau} (u_1 - v_1)$$
$$\frac{d u_2}{d t} = - \frac{1}{\tau} (u_2 - v_2)$$
where $v_1, v_2$ represent the right, left eye processing of the input; $u_1$, $u_2$ are the adaption parameters for the right and left eye; $\tau$ is the time constant to be fitted; the function $F$ models the firing rate input from the other cell and taken to be a sigmoid function, specifically
$$F(x) = \frac{1}{1 + \exp(1-x)}$$
In our case, we will ignore the adaption parameters $u_1$ and $u_2$ to make our system two-dimensional:
$$\frac{d v_1}{d t} = - v_1 + F(I_{app} - a v_2 )$$
$$\frac{d v_2}{d t} = - v_2 + F(I_{app} - a v_1)$$

See the simulation in `scripts/firing_rate_models.py`