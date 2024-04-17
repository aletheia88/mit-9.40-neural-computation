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

#### the shifting theorem
**Theorem:** Suppose $f(x)$ is an absolutely integrable and $\mathcal{F}\{ f(x)\} = \hat{f}(k)$, then $\mathcal\{f(x + c)\} = e^{ikc} \hat{f}(k)$.

##### Proof
From the definition of the Fourier Transform,
$$\mathcal{F\{f(x + c)\}} = \int_{-\infty}^{\infty} f(x + c) e^{-ikx} dx$$
Let $z= x + c$, so $dz = dx$.
$$= \int_{z=-\infty}^{\infty} f(z)e^{-ik(z-c)} dz$$
$$= e^{-ikc}\int_{z=-\infty}^{\infty} f(z)e^{-ikz} dz$$
$$= e^{-ikc}\hat{f(z)}$$

### FT of the Derivative

**Theorem:** Suppose $f(x)$ is a differentiable function, $f(x)$ and $\frac{df}{dx}$ are absolutely integrable, and $\lim_{x \to \pm \infty} f(x) = 0$, then $\mathcal{F}\left\{ \frac{df}{dx} \right\} = ik \hat{f}(k)$.

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
