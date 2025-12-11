---
layout: post
title: "Practical Implementation of Langevin Dynamics"
description: "An intuitive introduction to stochastic differential equations, starting from the Wiener process and progressing to Itô's lemma, the Fokker-Planck equation, and Langevin dynamics."
tags: Math MD
categories: Academic
date: 2025-12-10
toc:
  beginning: true
---

This post is a continuation of our previous discussion, where we demonstrated why Langevin Dynamics correctly samples from the canonical ensemble. Here, we move from theory to practice and explore how to actually implement Langevin Dynamics (often called a Langevin Integrator) in simulation.

Prerequisites: [Langevin Dynamics: From Wiener Process to Stochastic Differential Equations](/blog/2025/langevin-dynamics1/)

---

## Problem Setting

Recall the Langevin equation of motion:

$$
\left\{
\begin{aligned}
    dx_t &= v_t dt \\
    dv_t &= \frac{F(x_t)}{m}dt - \gamma v_t dt + \sqrt{\frac{2\gamma k_BT}{m}}dW_t
\end{aligned}
\right. \tag{1}
$$

A naive approach to solving this numerically would be to simply discretize the equation (the Euler-Maruyama method):

$$
v_{n+1} = v_n + \frac{F_n}{m}\Delta t - \gamma v_n \Delta t + \sqrt{\frac{2\gamma k_BT}{m}}\sqrt{\Delta t}\xi_n
$$

where $\xi_n$ is a random number drawn from a standard normal distribution.

While conceptually simple, this approach fails for molecular systems. It introduces significant integration errors that distort the sampling of the Boltzmann distribution, effectively simulating the wrong temperature or energy landscape. To fix this, we must avoid treating the stochastic term like a standard differential. The solution is to decompose the dynamics into distinct components and integrate them sequentially.

## The Solution: Operator Splitting

To solve this robustly, we break the complex Langevin equation into three distinct, analytically solvable parts. We imagine the system evolving under three separate "operators":

1. **A (Drift):** The system moves at constant velocity (position update).
    * $dx = v dt$
    * $dv = 0$

2. **B (Kick):** The velocities change due to conservative forces (velocity update).
    * $dx = 0$
    * $dv = \frac{F(x)}{m} dt$

3. **O (Ornstein-Uhlenbeck):** The velocities change due to the thermostat (friction + noise).
    * $dx = 0$
    * $dv = -\gamma v dt + \sqrt{\frac{2\gamma k_B T}{m}} dW$

While the first two parts are deterministic ODEs that are trivial to integrate, the third part is a Stochastic Differential Equation (SDE). We will now show how to solve this exactly.

### Mathematical Interlude: The Itô Isometry

Before solving the SDE, we must prove a useful result known as the **Itô Isometry**. This will allow us to calculate the variance of the stochastic noise. The theorem states:

$$
    \mathbb{E}\left[\left(\int_0^T X_t dW_t\right)^2\right] = \mathbb{E}\left[\int_0^T X_t^2 dt\right] \tag{2}
$$

**Proof:**

We define a partition over the time interval $[0,T]$ such that $0=t_0<t_1<\cdots < t_N=T$. The Itô integral is defined as the limit of the sum:

$$
    I = \int_0^T X_t dW_t \approx \sum_{i=0}^{N-1}X_{t_i}(W_{t_{i+1}}-W_{t_i})
$$

Let $\Delta W_i = W_{t_{i+1}}-W_{t_i}$. The square of the summation is:

$$
\begin{aligned}
    I^2 &= \left(\sum_{i=0}^{N-1}X_{t_i}\Delta W_i\right)\left(\sum_{j=0}^{N-1}X_{t_j}\Delta W_j\right)\\
    &=\sum_{i=0}^{N-1}X^2_{t_i}(\Delta W_i)^2 + 2\sum_{0\leq i<j\leq N-1}X_{t_i}X_{t_j}\Delta W_i\Delta W_j
\end{aligned}
$$

where the first summation contains the **diagonal terms** ($i=j$) and the second summation contains the **cross terms** ($i \neq j$).

Next, we take the expectation of the entire expression:

$$
    \mathbb{E}[I^2] = \mathbb{E}\left[\sum_{i=0}^{N-1}X^2_{t_i}(\Delta W_i)^2\right] + 2\mathbb{E}\left[\sum_{0\leq i<j\leq N-1}X_{t_i}X_{t_j}\Delta W_i\Delta W_j\right]
$$

For the **second expectation** (the cross terms), we condition on the filtration $\mathcal{F}_{t_j}$ (the information available at time $t_j$). Since $i < j$, the terms before $t_j$ are all known at time $t_j$. We also recall the independent increment property of the Wiener process:

$$
\begin{aligned}
\mathbb{E}\left[\sum_{0\leq i<j\leq N-1}X_{t_i}X_{t_j}\Delta W_i\Delta W_j\right] &= \sum_{0\leq i<j\leq N-1}\mathbb{E}\left[X_{t_i}X_{t_j}\Delta W_i\Delta W_j\right]\\
&= \sum_{0\leq i<j\leq N-1}\mathbb{E}\left[\mathbb{E}\left[X_{t_i}X_{t_j}\Delta W_i\Delta W_j \mid \mathcal{F}_{t_j}\right]\right]\\
&= \sum_{0\leq i<j\leq N-1}\mathbb{E}\left[X_{t_i}X_{t_j}\Delta W_i \cdot \underbrace{\mathbb{E}\left[\Delta W_j \mid \mathcal{F}_{t_j}\right]}_{0}\right]\\
&= 0
\end{aligned}
$$

For the **first expectation** (the diagonal terms), we condition on the filtration $\mathcal{F}_{t_i}$. We recall that the squared increment of a Wiener process has an expectation equal to the time step, $\mathbb{E}[(\Delta W_i)^2] = \Delta t$:

$$
\begin{aligned}
\mathbb{E}\left[\sum_{i=0}^{N-1}X^2_{t_i}(\Delta W_i)^2\right] &= \sum_{i=0}^{N-1}\mathbb{E}\left[X^2_{t_i}(\Delta W_i)^2\right]\\
&= \sum_{i=0}^{N-1}\mathbb{E}\left[\mathbb{E}\left[X^2_{t_i}(\Delta W_i)^2 \mid \mathcal{F}_{t_i}\right]\right]\\
&= \sum_{i=0}^{N-1}\mathbb{E}\left[X^2_{t_i}\underbrace{\mathbb{E}\left[(\Delta W_i)^2 \mid \mathcal{F}_{t_i}\right]}_{\Delta t}\right]\\
&= \mathbb{E}\left[\sum_{i=0}^{N-1}X^2_{t_i}\Delta t\right]
\end{aligned}
$$

The term inside the expectation is a Riemann sum. As the partition becomes finer ($N \to \infty$), it converges to the integral:

$$
    \mathbb{E}\left[\left(\int_0^T X_t dW_t\right)^2\right] = \mathbb{E}\left[\int_0^T X_t^2 dt\right]
$$

### Solution of the Ornstein-Uhlenbeck Process

Now we return to the "O" step. The SDE for the Ornstein-Uhlenbeck process is:

$$
    dv_t=-\gamma v_tdt+\sigma dW_t
$$

where $\sigma=\sqrt{\frac{2\gamma k_BT}{m}}$.

Similar to a linear ODE, we can use the method of integrating factors to solve this linear SDE.

$$
\begin{aligned}
&\quad\ \ dv_t=-\gamma v_tdt+\sigma dW_t\\
&\Rightarrow dv_t+\gamma v_tdt=\sigma dW_t\\
&\Rightarrow \mathrm{e}^{\gamma t}(dv_t+\gamma v_tdt)=\mathrm{e}^{\gamma t}\sigma dW_t\\
&\Rightarrow d(\mathrm{e}^{\gamma t}v_t)=\mathrm{e}^{\gamma t}\sigma dW_t\\
&\Rightarrow \int_{t}^{t+\Delta t}d(\mathrm{e}^{\gamma s}v_s)=\int_{t}^{t+\Delta t}\mathrm{e}^{\gamma s}\sigma dW_s\\
&\Rightarrow \mathrm{e}^{\gamma(t+\Delta t)}v(t+\Delta t)-\mathrm{e}^{\gamma t}v(t)=\sigma\int_{t}^{t+\Delta t}\mathrm{e}^{\gamma s} dW_s\\
&\Rightarrow v(t+\Delta t)=\underbrace{\mathrm{e}^{-\gamma \Delta t}v(t)}_{\mathrm{Deterministic\ Decay}}+\underbrace{\sigma\mathrm{e}^{-\gamma(t+\Delta t)}\int_{t}^{t+\Delta t}\mathrm{e}^{\gamma s} dW_s}_{\mathrm{Stochastic\ Integral}}\\
\end{aligned}
$$

The second term is a sum of Gaussian variables, which is itself a Gaussian. For simulation, we only need to determine its Mean and Variance.

Let $I = \sigma\mathrm{e}^{-\gamma(t+\Delta t)}\int_{t}^{t+\Delta t}\mathrm{e}^{\gamma s} dW_s$.

* **Mean:** $\mathbb{E}[I]=0$ (since $\mathbb{E}[dW_s]=0$).
* **Variance:** We use the Itô Isometry we proved earlier ($\mathbb{E}[(\int fdW)^2]=\mathbb{E}[\int f^2dt]$):

$$
\begin{aligned}
    \mathrm{Var}(I)&=\mathbb{E}\left[\left(\sigma\mathrm{e}^{-\gamma(t+\Delta t)}\int_{t}^{t+\Delta t}\mathrm{e}^{\gamma s} dW_s\right)^2\right]\\
    &=\sigma^2\mathrm{e}^{-2\gamma(t+\Delta t)}\int_t^{t+\Delta t}\mathrm{e}^{2\gamma s}ds\\
    &=\sigma^2\mathrm{e}^{-2\gamma(t+\Delta t)}\left[\frac{1}{2\gamma}\mathrm{e}^{2\gamma s}\right]_t^{t+\Delta t}\\
    &=\frac{\sigma^2}{2\gamma}(1-\mathrm{e}^{-2\gamma\Delta t})\\
    &=\frac{2\gamma k_BT/m}{2\gamma}(1-\mathrm{e}^{-2\gamma\Delta t})\\
    &=\frac{k_BT}{m}(1-\mathrm{e}^{-2\gamma\Delta t})
\end{aligned}
$$

Putting it all together, the **exact** update rule for the O-step is:

$$
    v(t+\Delta t) = \underbrace{e^{-\gamma \Delta t}}_{c_1} v(t) + \underbrace{\sqrt{\frac{k_B T}{m}(1 - e^{-2\gamma \Delta t})}}_{c_2} \xi
$$

where $\xi \sim \mathcal{N}(0,1)$.

## The BAOAB Algorithm

This derivation leads us to the **BAOAB** scheme (Kick-Drift-Noise-Drift-Kick). This algorithm is currently considered one of the most accurate schemes for sampling the Boltzmann distribution.

**Constants:**

* $c_1 = e^{-\gamma \Delta t}$
* $c_2 = \sqrt{\frac{k_B T}{m}(1 - c_1^2)}$

**The Loop (for each time step $\Delta t$):**

1. **B (Half Kick):**

    $$v \leftarrow v + \frac{\Delta t}{2} \frac{F(x)}{m}$$

2. **A (Half Drift):**

    $$x \leftarrow x + \frac{\Delta t}{2} v$$

3. **O (Full Noise/Friction):**

    $$v \leftarrow c_1 v + c_2 \xi \quad (\xi \sim \mathcal{N}(0,1))$$

4. **A (Half Drift):**

    $$x \leftarrow x + \frac{\Delta t}{2} v$$

5. **B (Half Kick):**

    $$\text{Recalculate forces } F(x)$$

    $$v \leftarrow v + \frac{\Delta t}{2} \frac{F(x)}{m}$$

## Conclusion

In this post, we finally bridged the gap between abstract math and running actual simulations. We saw that you can't just throw standard calculus at Langevin dynamics because thermal noise is too "jagged" for simple numerical methods to handle. By using **Operator Splitting** and the **BAOAB** scheme, we were able to isolate that messy stochastic part and solve it exactly using the **Ornstein-Uhlenbeck** process. This guarantees that our thermostat isn't just adding random interference, but is actually driving the system to the correct **Boltzmann distribution**, providing the rigorous math that makes tools like OpenMM work under the hood.
