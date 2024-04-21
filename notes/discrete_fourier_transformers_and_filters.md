## Discrete fourier transforms and filters (lecture 16)

### Time series and fourier transforms

#### Definition of the FT

The *Fourier transform* of a function $f: \mathcal{R} \to \mathcal{C}$ is defined to be
$$\mathcal{F}\{f\} = \hat{f}(k) = \int_{-\infty}^{\infty} f(x) e^{-ikx} dx$$
provided that the integral exists. To elaborate more on this criterion,

* if $f(x)$ is *absolutely integratable*, that is if $$\int_{-\infty}^{\infty} |f(x)| dx < \infty$$ then $\hat{f}(k)$ exists.
* if $f \to 0$ as $x \to \pm \infty$


#### The inverse of the FT
The inverse FT of a function $\hat{g}: \mathcal{R} \to \mathcal{C}$ is defined to be $$\mathcal{F}^{-1}\{\hat{g}\} = g(x) = \frac{1}{2\pi} \int_{-\infty}^{\infty} \hat{g}(k) e^{ikx} dx$$
provided that the integral exists.

#### Example

Compute the FT of $$f_2(x) = f_1(-x) = H(-x) e^x$$

**Solution:**

$$\mathcal{F}\{ H(-x)e^x \} = \int_{-\infty}^{\infty} H(-x)e^x  e^{-ikx} dx$$
$$= \int_{-\infty}^{0} H(x) \ e^{x(1-ik)} dx = \frac{e^{x(1-ik)}}{1 - ik} \Bigg|_{-\infty}^{0} = \frac{1}{1 - ik}$$

#### Properties of FT

##### reflection

If $\mathcal{F}\{f(x)\} = \hat{f}(k)$ exists, then $\mathcal{F\{f(-x)\}} = \hat{f}(-k)$

##### diation

If $\mathcal{F}\{f(x)\} = \hat{f}(k)$ exists, then for $a \in \mathbb{R}$, $$\mathcal{F\{f(ax)\}} = \frac{1}{|a|}\hat{f}\left( \frac{k}{a} \right)$$

##### linearity

The Fourier transform is a linear operator; that is for $a, b \in \mathbb{C}$, $$\mathcal{F\{af(x) + bf(x) \}} = a\hat{f}(k) + b\hat{f}(k)$$
assuming the FGT of $f, g$ both exist.

#### under convolution

The transformation of a convolution becomes simply the product of two functions without an integral!

$$\mathcal{F}[f * g] = \mathcal{F}[g] \mathcal{F}[f]$$

#### FT of Gaussian is Gaussian

Consider a Gaussian filter in time:
$$g(t) = \frac{1}{\sqrt{2 \pi \sigma^2}} e^{-t^2 / 2\sigma^2}$$
Its FT is
$$\mathcal{g(t)} = G(\nu) = \frac{1}{\sqrt{2 \pi \sigma^2}} \int_{-\infty}^{\infty} e^{i 2\pi \nu t} e^{- t^2/2\sigma^2} dt$$
or
$$G(\nu) = e^{-2 \pi^2 \nu^2 \sigma^2}\frac{1}{\sqrt{2 \pi \sigma^2}}\int_{-\infty}^{\infty} e^{- ( t - 2\pi i \sigma^2 \nu)^2/2\sigma^2} d t$$

The complex integral can be shown that it is equal to 1. Therefore, the important result is that the Fourier transform of a Gaussian is also a Gaussian
$$\boxed{G(\nu) = e^{-2 \pi^2 \sigma^2 \nu^2}}$$

#### the shifting theorem
**Theorem:** Suppose $f(x)$ is an absolutely integrable and $\mathcal{F}\{ f(x)\} = \hat{f}(k)$, then $\mathcal\{f(x + c)\} = e^{ikc} \hat{f}(k)$.

##### Proof
From the definition of the Fourier Transform,
$$\mathcal{F\{f(x + c)\}} = \int_{-\infty}^{\infty} f(x + c) e^{-ikx} dx$$
Let $z= x + c$, so $dz = dx$.
$$= \int_{z=-\infty}^{\infty} f(z)e^{-ik(z-c)} dz$$
$$= e^{-ikc}\int_{z=-\infty}^{\infty} f(z)e^{-ikz} dz$$
$$= e^{-ikc}\hat{f(z)}$$

### FT of the derivative

**Theorem:** Suppose $f(x)$ is a differentiable function, $f(x)$ and $\frac{df}{dx}$ are absolutely integrable, and $\lim_{x \to \pm \infty} f(x) = 0$, then $\mathcal{F} \left\{ \frac{df}{dx} \right\} = ik \hat{f}(k)$.

#### Proof
To solve $$\mathcal{F}\left\{ \frac{df}{dx} \right\} = \int_{-\infty}^{\infty} \frac{df}{dx} e^{-ikx} dx,$$
we use integration by parts; specifically, we let
$$u = e^{-ikx} \quad \quad dv = \frac{df}{dx} dx$$
which makes
$$du = (-ik) e^{-ikx} dx \quad \quad v = f(x)$$
Recall that the method of integration by parts states
$$\int_{x = -\infty}^{\infty} u dv = uv \Bigg|_{x = -\infty}^{\infty} - \int_{-\infty}^{\infty} v du$$

This means
$$\mathcal{F}\left\{ \frac{df}{dx} \right\} = \int_{-\infty}^{\infty} \frac{df}{dx} e^{-ikx} dx = f(x)e^{-ikx} \Bigg|_{x = -\infty}^{\infty} - \int_{-\infty}^{\infty} f(x)(-ik) e^{-ikx} dx = (ik) \int_{-\infty}^{\infty} f(x)e^{-ikx} dx = ik \hat{f}(k) $$
Note the boundary term $f(x)e^{-ikx} \Bigg|_{x = -\infty}^{\infty}$ disappears because $\lim_{x \to \pm \infty} f(x) = 0$.

FT turns differentiation into multiplication!

### Solving ODEs with the FT

#### Example:
Solve for $y(x)$ that satisfies
$$\textbf{DE:} \quad y' + y = H(x) e^{-2x}, -\infty < x < \infty$$
$$\textbf{BC:} \quad \lim_{x \to \pm \infty} y(x) = 0$$

**Solution:**

Let $\hat{y}(k) = \mathcal{F}\{ y(x) \}$ and Fourier transform the DE.
$$\mathcal{F}\{ y' + y \} = \mathcal{F}\{ H(x) e^{-2x}\}$$
$$\mathcal{F}\{ y' \} + \mathcal{F}\{ y \} = \frac{1}{2 + ik}$$
$$ik \hat{y}(k) + \hat{y}(k) = \frac{1}{2 + ik}$$
$$(1 + ik) \hat{y}(k) = \frac{1}{2 + ik}$$
$$\hat{y}(k) = \frac{1}{(1 + ik)(2 + ik)}$$

Now, we take the inverse FT:
$$y(x) = \mathcal{F}^{-1} \left\{ \frac{1}{(1 + ik)(2 + ik)} \right\} = \mathcal{F}^{-1} \left\{ \frac{1}{(1 + ik)} - \frac{1}{(2 + ik)} \right\}$$
$$= \mathcal{F}^{-1}\left\{ \frac{1}{(1 + ik)} \right\} -  \mathcal{F}^{-1}\left\{ \frac{1}{(2 + ik)} \right\} = H(x)e^{-x} - H(x)e^{-2x}$$

### $\delta$-function
A $\delta$-function mathematically models an impulse the transfer of a finite amount of momentum in an infinitesimal amount of time. It is not a function, but a distribution. It can be described as a limit but it is usually defined by its action on other functions.

#### the sampling property
$$f(x_0) = \int_{x=-\infty}^{\infty} f(x) \delta(x - x_0) dx$$
Using the sampling property,
$$\mathcal{F}\{\delta(x - x_0)\} = \int_{x=-\infty}^{\infty} \delta(x - x_0) e^{-ikx} dx = e^{-ikx_0}$$
At the origin ($x_0 = 0$):
$$\mathcal{F}\{ \delta(x) \} = \int_{-\infty}^{\infty} \delta(x) e^{-ikx} dx = 1$$

### The discrete (fast) Fourier Transform (FFT)
Define a function $g(t)$ on a finite set $[a, b]$, then the FT is $$G(\nu) = \int_a^b e^{i 2\pi \nu t} g(t) dt$$
We then center $g(t)$ around 0 and assume that $g(t)$ is zero outside of the bound $[-1/2, 1/2]$. The FT becomes:
$$G(\nu) = \int_{-1/2}^{1/2} e^{i 2\pi \nu t} g(t) dt$$
The discrete form of this is
$$G(\nu) = \Delta t \sum_{n = -N/2}^{N/2-1} e^{i 2\pi \nu n \Delta t} g(n\Delta t)$$
where $N = \frac{1}{\Delta t}$ and $g_n = g(n \Delta t)$. $N$ is called the sampling rate.

So far, we have assumed $\nu$ to be continous. To discretize $\nu$, instead of assuming that $g(t)$ is zero outside $[-1/2, 1/2]$, we assume that it has a period of $T=1$. Periodicity of $1$ implies that $g(t)$ has the same value at the current time point of samping and the next sampling point; in other words,
$$g_n = g_{n + N}$$

#### Nyquist frequency

The Nyquist frequency is the maximum frequency present in the approximation to the signal and any higher frequencies in the original signal will be missing. It is defined as
$$f_{\text{Nyquist}} = \frac{1}{2 \Delta t}$$

#### Sampling theorem
The *sampling theorem* states that any signal with frequencies less than $f_{Nyquist}$ is completely determined by its samples $g_n$.

To get an informative sampling of a signal, the sampling rate $\frac{1}{\Delta t}$ should be higher than twice the Nyquist frequency. For example, if the signal is a sinusoidal wave of frequency $f$, the sampling rate should be at least twice of $f$ or 2 sampled points per period $T = 1/f$.

### Power Spectral Density (PSD)

The total power in a signal can be caluclated either in the time domain or frequency domain. It is given by

$$\text{Power} = \int_{-\infty}^\infty | g(t)|^2 dt = \int_{-\infty}^\infty | G(\nu) |^2 d \nu$$

It can be easily derived as follows:

$$\int_{-\infty}^\infty | G(\nu) |^2 d \nu  = \int_{-\infty}^\infty  G(\nu) G^*(\nu) d \nu = \int d\nu \int e^{i 2 \pi \nu t} g(t) dt  \int e^{-i 2 \pi \nu t'} g^*(t') dt' =  \int dt dt' g(t) g^*(t') \delta(t -t')= \int dt | g(t) |^2$$

*(Note that the last step makes use of the inverse FT of the delta function)*

Consider a continous white noise signal $x(t)$ that satisfies
$$\langle x(t) \rangle = 0 \quad \langle x(t)x(\tau)\rangle = \delta(t - \tau)$$
The PSD of $x(t)$ is defined as $\langle | X(f) |^2 \rangle$

### Correlation

The correaltion between two signals $f(t)$ and $g(t)$ is given by
$$\text{Corr}(f, g)(t) = \int_{-\infty}^{\infty} f(t + \tau) g(\tau) d\tau$$
where $t$ represents the lag between $f(t)$ and $g(t)$. **Autocorrelation** is when $f = g$.

#### Wiener-Khinchin Theorem
The Fourier transform of the autocorrelation for a real signal $g(t)$ is the PSD of the signal.

$$\text{Corr}(g, g)(t) = \mathcal{F} \left[ \int_{-\infty}^{\infty} g(t + \tau)g(\tau) d\tau \right] = \int_{-\infty}^{\infty} dt \ e^{i2\pi \nu t} \int_{-\infty}^{\infty} g(t + \tau)g(\tau) d\tau$$
$$=\int_{-\infty}^{\infty} dt \ e^{i 2\pi \nu (t + \tau)} g(t + \tau) \int_{-\infty}^{\infty} e^{-i 2\pi \nu t} g(\tau) d\tau = G(\nu)G(-\nu) = |G(\nu)|^2$$

### White noise and Spike-triggered Average (STA)

In a linear system, the response $y(t)$ to the stimulus $x(t)$ is $$y(t) = \int K(t - t') x(t') dt$$

If $y(t)$ represents the measured spikes, then $y(t) = \sum_{s} \delta(t - t_s)$.

Let $t_s$ be the time of spike occuring and $K$ be the temporal receptive field (RF) of a neuron:

$$\text{STA}(\tau) = \frac{1}{N_s} \sum_s x(t_s - \tau) = \frac{1}{N_s} \sum_s \int dt \ x(t - \tau) \delta(t - t_s)$$
$$=\frac{1}{N_s} \int dt \ x(t - \tau) \sum_s \delta(t - t_s) = \frac{1}{N_s} \int dt \ y(t)x(t - \tau)$$
$$ = \frac{1}{N_s} \int dt \int dt' x(t - \tau) K(t - t') x(t')$$

Since the stimulus $x(t)$ is taken to be white noise, its average follows the property that $\langle x(t) x(t') \rangle = \delta(t - t')$. Substituting this to the above, we find that
$$\langle \text{STA} (\tau) \rangle = \frac{1}{N_s} \int dt \int dt' \ K(t -t') \delta(t - \tau - t') = \dots = \frac{T}{N} K(\tau)$$
where $T$ represents the number of sampled points in a discrete signal. STA and RF differ only by a constant.

### Filters in frequency space: inverse Fourier transform

