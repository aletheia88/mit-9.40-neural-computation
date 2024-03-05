The complete Hodgkin-Huxley (HH) model describes the interplay between the sodium and potassium flow in and out of the cell to create a spike event. The main ingredient is the equations that describe the time-dependence of the sodium and potassium channels. It is described by the following differential equations:
$$C_m \frac{dV}{dt} = -\bar{G_L}(V -E_L) - \bar{G_K} n^4 (V - E_K) - \bar{G_{Na}} m^3 h (V - E_{Na}) + I_e$$
The gating variables $m$, $n$, $h$ all satisfy
$$\tau_n (V) \frac{dn}{dt} = n_{\infty}(V) - n$$
Both potassium and sodium currents depend on their channel conductances, $$I_K = G_K(V)(V - E_K) \ \ \ I_{Na} = G_{Na}(V)(V - E_{Na})$$The conductances $G_K(V)$ and $G_{Na}(V)$ can be modeled as a sigmoid function.
$$S(V) = \frac{g_{max}}{ 1 + e^{-(V - V_{\text{offset}})/\sigma}}$$
The gates of the potassium and sodium channels change at certain time-dependent rate
$$\frac{dn}{dt} = -\beta_n (V) n \ + \ \alpha_n(V)(1-n)$$
$$\frac{dm}{dt} = -\beta_m (V) m \ + \ \alpha_m(V)(1-m)$$
where $\beta_{n, m} (V)$ is the rate of gate opening and $\alpha_{n,m}(V)$ is the rate of gate closing.

