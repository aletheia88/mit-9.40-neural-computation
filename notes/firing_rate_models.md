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

### noise induced transitions

Set $a_c = 5$ and $I_{app} = 3$. Let $\xi$ be a random source. Our 2D dynamical system now becomes $$\frac{d v_i}{dt} = - v_i + F(I - a v_j) + \xi$$

#### Euler-Maruyama (Euler) method
$$X_{n + 1} = X_n + f(X_n) \Delta t + \sqrt{\sigma^2 \Delta t} \ N(0, 1)$$

Using the Euler method, we update our system as follows:
$$v^{1}_{t+1} = v^{1,2}_{t} + dv^{1,2}_t$$
where the supscripts mean to represent that the value depends on both $v_1$ and $v_2$.
We write the update $d v^{1,2}_{t}$ as
$$d v^{1}_{t} = \frac{1}{\tau} \cdot \left( -v_{t}^{1} \ dt + F(I - a v_{t}^2) \ dt + \sqrt{\sigma^2 \Delta t} \ N(0, 1) \right)$$
where we group the random source as: $\sigma N(0, 1)\sqrt{\Delta t} \equiv x^i \sqrt{\Delta t}$
Substituting the temporal update yields
$$v_{t + 1}^i = v_{t}^i + \frac{1}{\tau} \ dt \left( - v_t^i + F(I - a v_t^j) \right) + x^i \sqrt{\Delta t}$$
This is implemented in `scripts/firing_rate_models.py`

### synaptic saturation with a train of Poisson spikes

Synaptic conductance can be modeled as a maximum conductance multiplied by a probability $P$ reflecting events at the presynaptic and postsynaptic terminal

$$g_s = \bar g_s P$$

We factorize the probability $P$ as the probability of release of a neurotransmitter and that of opening postsynaptic channels.

$$P = P_{post} P_{release}$$

The release probability decreases between spikes:
$$\tau_p \frac{d P_{release}}{dt} = P_{0, release} - P_{release}$$
where $P_{0, release}$ is the release probability in the absence of spikes (very low). We update $P_{release}$ after each spike as follows
$$P_{+, release} = P_{release} + \alpha (1 - P_{release}) \quad 0 \leq \alpha \leq 1$$

Let $T$ be the time in between the current and previous spike (let $P$ be $P_{release}$ for notation simplifity). Just before the new spike acts on vesicles, we have
$$P(T) = P_0 + (P + \alpha (1 - P) - P_0) e^{-T/\tau_P}$$