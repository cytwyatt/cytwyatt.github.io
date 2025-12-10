---
layout: post
title: "Langevin Dynamics: From Wiener Process to Stochastic Differential Equations"
description: "An intuitive introduction to stochastic differential equations, starting from the Wiener process and progressing to Itô's lemma, the Fokker-Planck equation, and Langevin dynamics."
tags: Math MD
categories: Academic
date: 2025-12-09
---

Welcome to my personal notes on **Stochastic Differential Equations (SDEs)**. We will begin with one of the simplest stochastic processes, the **Wiener process (or Brownian motion)**, and progress to **Itô's lemma** and the **Fokker-Planck equation**. By the end of this post, we will discuss how **Langevin dynamics** can be used to sample the correct **Boltzmann distribution** as a practical application of the Fokker-Planck equation.

Please note that this is not intended to be a rigorous mathematical treatise. Instead, my goal is to build an **intuitive understanding** of these theorems. Enjoy reading!

Prerequisites: Probability Theory, Molecular Dynamics

---

## Introduction

Welcome to my personal notes on **Stochastic Differential Equations (SDEs)**. We will begin with one of the simplest stochastic processes, the **Wiener process (or Brownian motion)**, and progress to **Itô's lemma** and the **Fokker-Planck equation**. By the end of this post, we will discuss how **Langevin dynamics** can be used to sample the correct **Boltzmann distribution** as a practical application of the Fokker-Planck equation.

Please note that this is not intended to be a rigorous mathematical treatise. Instead, my goal is to build an **intuitive understanding** of these theorems. Enjoy reading!

Prerequisites: Probability Theory, Molecular Dynamics

## Wiener Process

In probability theory, a **stochastic process** (or random process) is a collection of random variables indexed by time, denoted as $\{X(t)\}_{t\in T}$. Let's examine the example of a **1D random walk without drift**.

Starting from $X(0)=0$, we define the transition probabilities as:

$$
P(X(t+1)=x+1 \mid X(t)=x)=\frac{1}{2}
$$

$$
P(X(t+1)=x-1 \mid X(t)=x)=\frac{1}{2}
$$

This means that given the current position, the process moves either one unit length to the right or one unit length to the left in the next step. We can formalize this process by introducing independent and identically distributed (i.i.d.) random variables $Y_t$, defined as:

$$
Y_t \in \{-1, 1\} \quad \text{with} \quad P(Y_t=1) = P(Y_t=-1) = \frac{1}{2}
$$

Using this, we can express the position $X(T)$ as a cumulative sum of the steps $Y_i$:

$$
\begin{aligned}
    X(T) &= X(T-1) + Y_T \\
         &= X(T-2) + Y_T + Y_{T-1} \\
         &=\cdots \\
         &= \sum_{i=1}^{T}Y_i
\end{aligned}
$$

We observe that the individual steps have the following properties:

$$
    \mathbb{E}[Y_i] = 0 \quad \text{and} \quad \mathrm{Var}[Y_i] = 1
$$

Consequently, the expectation and variance of the cumulative sum are:

$$
    \mathbb{E}\left[\sum_{i=1}^T Y_i\right] = 0 \quad \text{and} \quad \mathrm{Var}\left[\sum_{i=1}^T Y_i\right] = T
$$

According to the **Central Limit Theorem**, as $T$ becomes large, the cumulative sum converges in distribution to a Gaussian. Therefore, $X(T)$ is approximately normally distributed with mean $0$ and variance $T$:

$$
    X(T) \sim \mathcal{N}(0, T)
$$

The **Wiener process** (or Brownian motion) can be viewed as the continuous limit of a random walk where the time step and spatial step size approach zero.

Consider a total time interval $T$ divided into $N$ equal sub-intervals of length $h$, such that:

$$
    T=N h
$$

To ensure the process converges to a finite variance, we define the spatial step size for each transition as $\sigma \sqrt{h}$. The transition probabilities become:

$$
P(X(t+h)=x+\sigma\sqrt{h} \mid X(t)=x)=\frac{1}{2}
$$

$$
P(X(t+h)=x-\sigma\sqrt{h} \mid X(t)=x)=\frac{1}{2}
$$

We can retain our previous definition of the i.i.d. random variables $Y_i \in \{-1, 1\}$. The position $X(T)$ can then be expressed as:

$$
    \begin{aligned}
        X(T) = X(Nh) &= X((N-1)h) + \sigma \sqrt{h}Y_N \\
            &= \sigma \sqrt{h}\sum_{i=1}^N Y_i
    \end{aligned}
$$

Using the properties of expectation and variance ($\text{Var}(cY) = c^2\text{Var}(Y)$), and applying the Central Limit Theorem for large $N$:

$$
    \begin{aligned}
        X(T) &\sim \sigma \sqrt{h} \, \mathcal{N}(0, N) \\
             &= \mathcal{N}(0, \sigma^2 h N) \\
             &= \mathcal{N}(0, \sigma^2 T)
    \end{aligned}
$$

By setting $\sigma=1$, we can formally define the **standard Wiener process**, denoted as $W_t$, with the following properties:

* **Initialization:** $W_0=0$

* **Normality:** $W_t$ is normally distributed with mean $0$ and variance $t$:

    $$
    W_t \sim \mathcal{N}(0, t)
    $$

* **Independent Increments:** For any disjoint time intervals $(t_1, t_2)$ and $(s_1, s_2)$, the increments $W_{t_2}-W_{t_1}$ and $W_{s_2}-W_{s_1}$ are independent random variables.

* **Stationary Increments:** The probability distribution of the increment $W_{t+s}-W_{s}$ depends only on the time difference $t$. It follows the same distribution as $W_t$:

    $$
    W_{t+s}-W_{s} \sim \mathcal{N}(0, t)
    $$

An interesting and defining property of the Wiener process is its covariance function, $\mathrm{cov}(W_{t_1}, W_{t_2})$. Since the expectation $\mathbb{E}[W_t] = 0$ for all $t$, the covariance definition simplifies to the expectation of the product:

$$
    \mathrm{cov}(W_{t_1}, W_{t_2}) = \mathbb{E}[(W_{t_1}-\mathbb{E}[W_{t_1}])\cdot (W_{t_2}-\mathbb{E}[W_{t_2}])] = \mathbb{E}[W_{t_1}\cdot W_{t_2}]
$$

To evaluate this, let's assume without loss of generality that $t_1 < t_2$. We can use a standard trick of rewriting $W_{t_2}$ as the sum of $W_{t_1}$ and the increment $(W_{t_2} - W_{t_1})$:

$$
    \begin{aligned}
        \mathbb{E}[W_{t_1}\cdot W_{t_2}] &= \mathbb{E}[W_{t_1}\cdot((W_{t_2}-W_{t_1})+W_{t_1})] \\
        &= \mathbb{E}[W_{t_1}\cdot(W_{t_2}-W_{t_1})] + \mathbb{E}[W_{t_1}^2]
    \end{aligned}
$$

Using the property of **Independent Increments**, we know that the future increment $(W_{t_2}-W_{t_1})$ is independent of the current state $W_{t_1}$. Therefore, the expectation of their product splits:

$$
    \begin{aligned}
        \mathbb{E}[W_{t_1}\cdot(W_{t_2}-W_{t_1})] &= \mathbb{E}[W_{t_1}] \cdot \mathbb{E}[W_{t_2}-W_{t_1}] \\
        &= 0 \cdot 0 \\
        &= 0
    \end{aligned}
$$

This leaves us with only the variance term. Since $W_{t_1} \sim \mathcal{N}(0, t_1)$, we have $\mathbb{E}[W_{t_1}^2] = t_1$.

Thus, combining both cases ($t_1 < t_2$ and $t_2 < t_1$), the general formula for the covariance is:

$$
    \mathbb{E}[W_{t_1}\cdot W_{t_2}] = \min(t_1, t_2)
$$

Having established the Wiener process, we can now define $dW_t$ as the increment of the Wiener process over a tiny time step $dt$:

$$
    dW_t = W_{t+dt}-W_{t}\sim \mathcal{N}(0, dt)
$$

This leads to the most important intuition in stochastic calculus: since the variance is $dt$, the standard deviation—representing the typical "size" of the step—scales as $\sqrt{dt}$.

$$
    dW_t \propto \sqrt{dt}
$$

Consequently, the squared increment scales linearly with time:

$$
    (dW_t)^2 \propto dt
$$

This is a strange property that fundamentally distinguishes stochastic calculus from standard calculus. In classical deterministic calculus, for a smooth function $f(x)$, the Taylor expansion is:

$$
    df(x) = f'(x)dx + o(dx)
$$

If we square this increment, we get terms of order $(dx)^2$:

$$
    (df(x))^2 \propto (dx)^2 = o(dx)
$$

In the limit where $dx \to 0$, $(dx)^2$ becomes negligible compared to $dx$, allowing us to discard second-order terms. However, for a Wiener process, **$(dW_t)^2$ is of the same scale as $dt$**. It is *not* negligible.

As a consequence, we must retain second-order terms in our expansions. This leads to the specific multiplication rules of Stochastic Calculus:

* **The Expectation is Zero:**
    $$
    \mathbb{E}[dW_t]=0
    $$

* **The Variance (The Non-Zero Quadratic Variation):**
    $$
    \mathbb{E}[(dW_t)^2]=dt
    $$
    For the purpose of calculation within integrals, we treat this deterministically:
    $$
    (dW_t)^2 = dt
    $$

* **Cross Terms Vanish:**
    Mixed terms involving higher powers of $dt$ become negligible:
    $$
    dW_t \cdot dt = 0
    $$

Note that though we can define the infinitesimal increment, the Wiener process is nowhere differentiable, because

$$
    \frac{dW}{dt}\approx \frac{W(t+h)-W(t)}{h}\approx \frac{\sqrt{h}}{h}=\frac{1}{\sqrt{h}}
$$

which is diverged when $h\to0$.

## Stochastic Differential Equations

### Definition and Interpretation

We define a general **Stochastic Differential Equation (SDE)** for a variable $X_t$ as:

$$
    dX_t=\mu(X_t,t)dt + \sigma(X_t,t)dW_t \tag{1}
$$

This equation describes how the system will evolve over an infinitesimal time step $dt$. It consists of two parts:

* Drift Term $\mu dt$: This is the deterministic part of the motion, which represents the standard laws of physics (like force) that drive the system.
* Diffusion Term $\sigma dW_t$: This is the stochastic part. It represents the random noise or fluctuations interacting with the system.

Note that, since $\frac{dW}{dt}$ is undefined, the SDE is actually shorthand notation for an integral equation:

$$
    X_t=X_0+\int_0^t\mu(X_s,s)ds+\int_0^t\sigma(X_s,s)dW_s \tag{2}
$$

The first integral (Riemann) is standard, but the result of the second integral (stochastic) depends entirely on the specific definition of the integration rule we choose. 

### Itô vs Stratonovich Integrals

To demonstrate this, let's analyze the integral of the Wiener process against itself, $\int_0^T W_t dW_t$, using a discrete Riemann sum approximation:

$$
    I = \int_0^T W_t dW_t \approx \sum_{i=0}^{N-1} W_{\tau_i}(W_{t_{i+1}}-W_{t_i})
$$

In standard calculus, the result would simply be $\frac{1}{2}W_T^2$. However, in stochastic calculus, the result depends on how we choose the evaluation point $\tau_i$ within each interval. Suppose we define $W_{\tau_i}$ as a weighted average of the two endpoints, controlled by a parameter $\alpha \in [0,1]$:

$$
    W_{\tau_i} = (1-\alpha)W_{t_i} + \alpha W_{t_{i+1}}
$$

Substituting this into our sum, we get:

$$
\begin{aligned}
    I &= \sum_{i=0}^{N-1} \left[(1-\alpha)W_{t_i} + \alpha W_{t_{i+1}}\right](W_{t_{i+1}}-W_{t_i}) \\
      &= \sum_{i=0}^{N-1} \left[(1-\alpha)W_{t_i} + \alpha (W_{t_i}+\Delta W_i)\right]\Delta W_i \\
      &= \sum_{i=0}^{N-1} (W_{t_i} + \alpha \Delta W_i)\Delta W_i \\
      &= \sum_{i=0}^{N-1} W_{t_i}\Delta W_i + \alpha \sum_{i=0}^{N-1}(\Delta W_i)^2 \\
      &= I_1 + I_2
\end{aligned}
$$

We can analyze the limits of these two terms separately as $N \to \infty$:

**1. The First Term ($I_1$)**
Using the algebraic identity $a(b-a) = \frac{1}{2}(b^2-a^2) - \frac{1}{2}(b-a)^2$, we can expand this sum:

$$
\begin{aligned}
    I_1 &= \sum_{i=0}^{N-1} W_{t_i}\Delta W_i \\
        &= \sum_{i=0}^{N-1}\frac{1}{2}(W_{t_{i+1}}^2-W_{t_i}^2) - \sum_{i=0}^{N-1}\frac{1}{2}(W_{t_{i+1}}-W_{t_i})^2 \\
        &= \frac{1}{2}W_T^2 - \frac{1}{2}T \quad (\text{via Telescoping Sum and Quadratic Variation})
\end{aligned}
$$

**2. The Second Term ($I_2$)**
This term depends on the quadratic variation of the process:

$$
    I_2 = \alpha \sum_{i=0}^{N-1}(\Delta W_i)^2 \to \alpha T
$$

**The General Result**
Combining these, the value of the stochastic integral is:

$$
    I = \frac{1}{2}W_T^2 + \left(\alpha-\frac{1}{2}\right)T
$$

This leads to two distinct interpretations:

* **Itô Integral ($\alpha=0$):** We evaluate at the left endpoint (past information).
    $$
    I = \frac{1}{2}W_T^2 - \frac{1}{2}T
    $$
    There is an additional correction term ($-\frac{1}{2}T$). Because of the existence of this term, standard calculus rules fail, and we must modify the chain rule, which leads to **Itô's Lemma**.

* **Stratonovich Integral ($\alpha=\frac{1}{2}$):** We evaluate at the midpoint.
    $$
    I = \frac{1}{2}W_T^2 + \left(\frac{1}{2}-\frac{1}{2}\right)T = \frac{1}{2}W_T^2
    $$
    The result matches standard calculus.

For the remainder of these notes, we will primarily use the **Itô integral**, as it respects causality.

## Itô's Lemma

### Derivation

Given a process $X_t$ governed by a standard SDE:

$$
    dX_t = \mu(X_t,t)dt + \sigma(X_t,t)dW_t
$$

Suppose we are interested in analyzing a function $f(X_t,t)$ that depends on the state $X_t$ and time $t$. To find the differential $df(X_t,t)$, we must determine how $f$ evolves as $X_t$ changes.

Just like in standard calculus, we perform a Taylor expansion. However, a crucial difference arises here. Because $(dX_t)^2$ contains terms of order $dt$ (unlike in deterministic calculus where $(dx)^2 \approx 0$), we **must** expand up to the second order to capture all terms affecting the drift:

$$
    df = \frac{\partial f}{\partial t}dt + \frac{\partial f}{\partial x}dX_t + \frac{1}{2}\frac{\partial^2f}{\partial x^2}(dX_t)^2 + \cdots
$$

To simplify this, let's examine the squared term $(dX_t)^2$:

$$
    \begin{aligned}
        (dX_t)^2 &= [\mu dt + \sigma dW_t]^2 \\
                 &= \mu^2(dt)^2 + 2\mu\sigma(dt)(dW_t) + \sigma^2(dW_t)^2
    \end{aligned}
$$

Applying the **Itô multiplication rules** (where terms higher than order $dt$ vanish, and $(dW_t)^2 \to dt$):

$$
    \begin{aligned}
        (dX_t)^2 &= 0 + 0 + \sigma^2(dt) \\
                 &= \sigma^2 dt
    \end{aligned}
$$

Substituting this result back into the Taylor expansion, we arrive at the fundamental formula of Stochastic Calculus, **Itô's Lemma**:

>$$
    df(X_t,t) = \left( \frac{\partial f}{\partial t} + \mu\frac{\partial f}{\partial x} + \frac{1}{2}\sigma^2\frac{\partial^2f}{\partial x^2} \right)dt + \sigma \frac{\partial f}{\partial x}dW_t \tag{3}
$$

The term $\frac{1}{2}\sigma^2\frac{\partial^2 f}{\partial x^2}$ is known as the **Itô correction**. It represents the drift induced by the convexity of the function $f$ combined with the variance of the process.

### Converting Stratonovich to Itô

We can also convert a Stratonovich SDE to an Itô SDE. Let's denote the **Stratonovich SDE** as:

$$
    dX_t = \mu(X_t,t)dt + \sigma(X_t,t)\circ dW_t \tag{4}
$$

The difference arises from how we approximate the stochastic integral term $\int \sigma \circ dW$. Recall that the Stratonovich integral evaluates $\sigma$ at the midpoint $\bar{X}=\frac{X_t+X_{t+dt}}{2}$. We can approximate $\sigma (\bar{X})$ using a first-order Taylor expansion around $X_t$:

$$
    \sigma(\bar{X}) \approx \sigma(X_t) + \sigma'(X_t)(\bar{X}-X_t) \approx \sigma(X_t) + \sigma'(X_t)\frac{1}{2}dX_t
$$

Substituting this approximation back into the noise term, and recalling that $dX_t \approx \sigma(X_t)dW_t$ for the leading order terms:

$$
\begin{aligned}
    \sigma(X_t,t)\circ dW_t = \sigma(\bar{X})dW_t &\approx \left[\sigma(X_t) + \frac{1}{2}\sigma'(X_t)(\mu dt + \sigma(X_t)dW_t)\right]dW_t \\
    &= \sigma(X_t)dW_t + \underbrace{\frac{1}{2}\mu\sigma'(dt)(dW_t)}_{\approx 0} + \frac{1}{2}\sigma(X_t)\sigma'(X_t)\underbrace{(dW_t)^2}_{dt} \\
    &= \sigma(X_t)dW_t + \frac{1}{2}\sigma(X_t)\sigma'(X_t)dt
\end{aligned}
$$

Thus, we have the relationship:

$$
    \sigma\circ dW_t = \sigma dW_t + \frac{1}{2}\sigma\sigma'dt
$$

Consequently, the original Stratonovich SDE

$$
    dX_t = \mu(X_t,t)dt + \sigma(X_t,t)\circ dW_t
$$

is equivalent to the **Itô SDE**:

$$
    dX_t = \left(\mu(X_t,t) + \frac{1}{2}\sigma(X_t,t)\frac{\partial\sigma}{\partial x}\right)dt + \sigma(X_t,t)dW_t \tag{5}
$$

## Fokker-Planck Equation

### Derivation

As a direct application of Itô's lemma, we can derive a deterministic partial differential equation that describes the time evolution of the probability density, starting from the SDE:

$$
    dX_t = \mu(X_t,t)dt + \sigma(X_t,t)dW_t
$$

While the trajectory of a single particle $X_t$ is random and jagged, the evolution of the **probability distribution** $p(x,t)$ of finding the particle at a specific location is deterministic and smooth. To derive this relationship, we consider the time evolution of an arbitrary smooth test function $f(X_t)$, which doesn't explicitly depend on time.

By Itô's lemma, we can write

$$
\begin{aligned}
    df(X_t) &= \frac{\partial f}{\partial t}dt + \frac{\partial f}{\partial x}dX_t + \frac{1}{2}\frac{\partial^2 f}{\partial x^2}(dX_t)^2 \\
    &= \left(\mu(x,t) \frac{\partial f}{\partial x} + \frac{1}{2}\sigma^2(x,t) \frac{\partial^2 f}{\partial x^2} \right)dt + \sigma \frac{\partial f}{\partial x} dW_t
\end{aligned}
$$

To find the equation for the probability density $p(x,t)$, we take the expectation $\mathbb{E}[\cdot]$ of the entire equation. Since the expectation of the Wiener increment is zero ($\mathbb{E}[dW_t]=0$), the stochastic term vanishes:

$$
    \frac{d}{dt}\mathbb{E}[f(X_t)] = \mathbb{E}\left[ \mu(X_t,t)\frac{\partial f}{\partial x} + \frac{1}{2}\sigma^2(X_t,t)\frac{\partial^2 f}{\partial x^2} \right]
$$

We can rewrite these expectations as integrals over the state space. By definition, $\mathbb{E}[g(X)] = \int g(x)p(x,t)dx$. Therefore:

$$
    \int_{-\infty}^{\infty} f(x) \frac{\partial p(x,t)}{\partial t} dx = \int_{-\infty}^{\infty} \left[ \mu(x,t)\frac{\partial f}{\partial x} + \frac{1}{2}\sigma^2(x,t)\frac{\partial^2 f}{\partial x^2} \right] p(x,t) dx
$$

The goal is to isolate $f(x)$. We use **Integration by Parts** to move the spatial derivatives from the test function $f$ onto the density $p$. We assume that the probability density and its flux vanish at the boundaries ($x \to \pm \infty$).

**1. The Drift Term (One integration by parts):**

$$
    \int_{-\infty}^{\infty} \left(\mu p\right) \frac{\partial f}{\partial x} dx = \underbrace{\left[ \mu p f \right]_{-\infty}^{\infty}}_{0} - \int_{-\infty}^{\infty} f(x) \frac{\partial}{\partial x}(\mu p) dx
$$

**2. The Diffusion Term (Two integrations by parts):**

$$
\begin{aligned}
    \int_{-\infty}^{\infty} \left(\frac{1}{2}\sigma^2 p\right) \frac{\partial^2 f}{\partial x^2} dx &= \underbrace{\left[ \dots \right]_{-\infty}^{\infty}}_{0} - \int_{-\infty}^{\infty} \frac{\partial}{\partial x}\left(\frac{1}{2}\sigma^2 p\right) \frac{\partial f}{\partial x} dx \\
    &= \underbrace{\left[ \dots \right]_{-\infty}^{\infty}}_{0} + \int_{-\infty}^{\infty} f(x) \frac{\partial^2}{\partial x^2}\left(\frac{1}{2}\sigma^2 p\right) dx
\end{aligned}
$$

Substituting these back into the main equation and grouping everything under one integral:

$$
    \int_{-\infty}^{\infty} f(x) \left\{ \frac{\partial p}{\partial t} + \frac{\partial}{\partial x}[\mu(x,t) p] - \frac{\partial^2}{\partial x^2}\left[\frac{1}{2}\sigma^2(x,t) p\right] \right\} dx = 0
$$

Since this integral must equal zero for **any** arbitrary smooth test function $f(x)$, the term inside the brackets must be identically zero.

This yields the **Fokker-Planck Equation** (also known as the Kolmogorov Forward Equation):

>$$
    \frac{\partial p(x,t)}{\partial t} = -\frac{\partial}{\partial x}\left[ \mu(x,t) p(x,t) \right] + \frac{\partial^2}{\partial x^2}\left[ \frac{1}{2}\sigma^2(x,t) p(x,t) \right] \tag{6}
$$

This PDE describes the deterministic time evolution of the probability distribution of a particle undergoing Brownian motion.

In physics applications, it is standard to define the **Diffusion Coefficient** $D(x,t)$ as half of the variance rate:

$$
    D(x,t) = \frac{1}{2}\sigma^2(x,t)
$$

Substituting this into our result, the Fokker-Planck equation becomes:

>$$
    \frac{\partial p(x,t)}{\partial t} = -\frac{\partial}{\partial x}\left[ \mu(x,t) p(x,t) \right] + \frac{\partial^2}{\partial x^2}\left[ D(x,t) p(x,t) \right] \tag{7}
$$

This equation reveals a deep physical insight: it is exactly a **Continuity Equation**. In hydrodynamics and electromagnetism, a continuity equation expresses the conservation of a quantity (mass, charge, or in this case, probability). It states that the rate of change of a density $\rho$ is equal to the negative divergence of its corresponding flux $\mathbf{J}$:

$$
    \frac{\partial \rho}{\partial t} = -\nabla \cdot \mathbf{J} \tag{8}
$$

In the 1D case of the Fokker-Planck equation, we can identify the **Probability Flux**, denoted as $J(x,t)$, by rewriting the equation as:

$$
    \frac{\partial p}{\partial t} + \frac{\partial}{\partial x} J(x,t) = 0
$$

By comparing terms, we define the flux $J(x,t)$ as the sum of two distinct physical transport mechanisms:

$$
    J(x,t) = \underbrace{\mu(x,t) p(x,t)}_{\text{Drift / Convection}} - \underbrace{\frac{\partial}{\partial x}\left[ D(x,t) p(x,t) \right]}_{\text{Diffusion}} \tag{9}
$$

1.  **Deterministic Convection (Drift):** The term $\mu p$ represents the probability mass being carried along by the deterministic velocity field $\mu$.
2.  **Diffusion:** The term $-\frac{\partial}{\partial x}(Dp)$ represents the spreading of probability due to random fluctuations (noise).

### Stationary Distribution

One of the most important questions we can ask about a stochastic process is: **"What happens after a long time?"**

If we let the system evolve for $t \to \infty$, the probability distribution often settles into a fixed shape that no longer changes with time. This is called the **Stationary Distribution**, denoted as $p_{st}(x)$.

Mathematically, "stationary" means the time derivative is zero:

$$
    \frac{\partial p(x,t)}{\partial t} = 0
$$

Substituting this into the Fokker-Planck equation, we get:

$$
    0 = -\frac{\partial}{\partial x} J(x)
$$

This implies that the probability flux $J(x)$ must be constant everywhere in space. For a physical system confined to a finite region (where the probability vanishes at $\pm \infty$), this constant flux must be **zero**.

$$
    J(x) = \mu(x)p_{st}(x) - \frac{d}{dx}[D(x)p_{st}(x)] = 0
$$

This is the condition of **Detailed Balance**: the deterministic drift current exactly cancels the diffusive spreading current at every point in space.

#### Solving for the Distribution
This is now a first-order Ordinary Differential Equation (ODE). We can rearrange it to solve for $p_{st}(x)$:

$$
    \frac{d}{dx} [D(x) p_{st}(x)] = \mu(x) p_{st}(x)
$$

To solve this, let's define a helper function $f(x) = D(x)p_{st}(x)$. The equation becomes $f'(x) = \frac{\mu(x)}{D(x)}f(x)$, or via separation of variables:

$$
    \frac{f'(x)}{f(x)} = \frac{\mu(x)}{D(x)}
$$

Integrating both sides and exponentiating gives the general solution:

$$
    p_{st}(x) = \frac{\mathcal{N}}{D(x)} \exp\left( \int_{x_0}^x \frac{\mu(y)}{D(y)} dy \right) \tag{10}
$$

*(Where $\mathcal{N}$ is a normalization constant).*

## Langevin Dynamics

### Motivation and Equation of Motion

Armed with this theoretical framework, we can finally examine **Langevin Dynamics**. This method is a cornerstone of modern Molecular Dynamics (MD) simulations and is widely used to sample the correct Boltzmann distribution in the canonical ($NVT$) ensemble.

To understand the motivation, consider the standard Newton's equation of motion:

$$
    m d\mathbf{v}_t = -\nabla U(\mathbf{x}_t) dt
$$

If we integrate this equation directly, the total energy of the system is conserved. Consequently, we sample the **Microcanonical ($NVE$) ensemble** (fixed particles, volume, and energy). However, this is often undesirable because, in experimental reality, it is much easier to control **temperature** than total energy. Therefore, we usually aim to simulate the Canonical ($NVT$) ensemble.

Langevin dynamics achieves this by mathematically coupling the system to a "heat bath." The standard equation of motion is augmented with two additional terms: a friction term (representing viscous drag from the solvent) and a noise term (representing random collisions with solvent molecules).

$$
    m d\mathbf{v}_t = \underbrace{-\nabla U(\mathbf{x}_t) dt}_{\text{Conservative Force}} - \underbrace{m\gamma \mathbf{v}_t dt}_{\text{Friction}} + \underbrace{\sqrt{2m\gamma k_B T} d\mathbf{W}_t}_{\text{Random Noise}} \tag{11}
$$

### Equipartition Theorem for Free Particle

First, let us demonstrate how the addition of these two terms enables the equation to sample from a fixed-temperature ensemble.

From the perspective of statistical mechanics, temperature is a macroscopic metric that measures the average kinetic energy of the microscopic particles in the system. This relationship is quantified by the **Equipartition Theorem**, which states that each quadratic degree of freedom contributes $\frac{1}{2}k_BT$ to the average energy.

In one dimension, this relationship is expressed as:

$$
    \frac{1}{2}m\langle v^2\rangle = \frac{1}{2}k_BT \tag{12}
$$

To demonstrate this, let us consider a one-dimensional **free particle** (where the potential $U(x)=0$). The velocity evolves according to the Langevin SDE:

$$
    dv_t = -\gamma v_t dt + \sqrt{\frac{2\gamma k_BT}{m}}dW_t
$$

We are interested in the average kinetic energy, which depends on $v^2$. Let us define $f(v)=v^2$ and find its differential using **Itô's lemma**:

$$
    df = d(v^2) = f'(v)dv + \frac{1}{2}f''(v)(dv)^2
$$

Knowing that $f'(v)=2v$ and $f''(v)=2$, and using the Itô rule $(dv_t)^2 = \frac{2\gamma k_B T}{m}dt$ (since we only keep terms of order $dt$), we substitute these into the expansion:

$$
\begin{aligned}
    d(v^2) &= 2v\left(-\gamma v_t dt + \sqrt{\frac{2\gamma k_BT}{m}}dW_t\right) + \frac{1}{2}(2)\left( \frac{2\gamma k_BT}{m}dt \right) \\
    &= \left(-2\gamma v^2 + \frac{2\gamma k_BT}{m}\right)dt + 2v\sqrt{\frac{2\gamma k_BT}{m}}dW_t
\end{aligned}
$$

Next, we take the expectation $\mathbb{E}[\cdot]$ of both sides. Since the expectation of the Wiener increment is zero ($\mathbb{E}[dW_t]=0$), the stochastic term vanishes:

$$
    \frac{d}{dt}\mathbb{E}[v^2] = -2\gamma\mathbb{E}[v^2] + \frac{2\gamma k_BT}{m}
$$

We are interested in the **equilibrium limit** ($t\to\infty$), where the system reaches a steady state and the time derivative becomes zero ($\frac{d}{dt}=0$):

$$
    0 = -2\gamma\langle v^2\rangle + \frac{2\gamma k_BT}{m}
$$

Solving for the mean squared velocity $\langle v^2 \rangle$:

$$
    \langle v^2\rangle = \frac{k_BT}{m}
$$

Multiplying both sides by $\frac{1}{2}m$, we recover the exact result predicted by the **Equipartition Theorem**:

$$
    \frac{1}{2}m\langle v^2\rangle = \frac{1}{2}k_BT
$$

### Sampling the Boltzmann Distribution

Finally, we demonstrate that in the presence of a conservative force, Langevin dynamics correctly samples the target **Boltzmann distribution** $p_{eq}(x,v)$, defined as:

$$
    p_{eq}(x,v) = \mathcal{Z}^{-1} \exp \left(-\frac{\frac{1}{2}mv^2+U(x)}{k_BT}\right)
$$

where $\mathcal{Z}$ is the partition function.

Note that when a potential $U(x)$ is present, the Langevin dynamics are described by a system of coupled Stochastic Differential Equations (SDEs):

$$
\left\{
\begin{aligned}
    dx_t &= v_t dt \\
    dv_t &= -\frac{U'(x_t)}{m}dt - \gamma v_t dt + \sqrt{\frac{2\gamma k_BT}{m}}dW_t
\end{aligned}
\right.
$$

To analyze the evolution of the probability density $p(x,v,t)$, we write the corresponding Fokker-Planck equation (often called the **Klein-Kramers equation**) in the full 2D phase space:

$$
\begin{aligned}
    \frac{\partial p}{\partial t} &= -\frac{\partial }{\partial x}(vp) - \frac{\partial}{\partial v}\left[\left(-\frac{U'(x)}{m}-\gamma v\right)p\right] + \frac{\gamma k_BT}{m}\frac{\partial^2p}{\partial v^2} \\
    &= \underbrace{-\frac{\partial}{\partial x}(vp) + \frac{\partial}{\partial v}\left(\frac{U'(x)}{m}p\right)}_{\text{Hamiltonian Flow}} + \underbrace{\gamma\frac{\partial}{\partial v}(vp) + \frac{\gamma k_BT}{m}\frac{\partial^2p}{\partial v^2}}_{\text{Thermostat / Dissipation}}
\end{aligned}
$$

### Proof of Stationary Distribution

To prove that $p_{eq}(x,v)$ is the stationary solution, we must show that substituting $p_{eq}$ into the equation makes the time derivative zero ($\frac{\partial p}{\partial t} = 0$).

First, let us calculate the partial derivatives of the proposed Boltzmann distribution:

1.  **Spatial Derivative:**

    $$\frac{\partial p_{eq}}{\partial x} = p_{eq} \cdot \frac{\partial}{\partial x}\left(-\frac{U(x)}{k_B T}\right) = -\frac{U'(x)}{k_B T} p_{eq}$$

2.  **Velocity Derivative:**

    $$\frac{\partial p_{eq}}{\partial v} = p_{eq} \cdot \frac{\partial}{\partial v}\left(-\frac{mv^2}{2k_B T}\right) = -\frac{mv}{k_B T} p_{eq}$$

Now, we substitute these into the two distinct parts of the Klein-Kramers equation.

**Part 1: The Hamiltonian Terms (Reversible Dynamics)**
Since $x$ and $v$ are independent phase space coordinates, $v$ can be pulled out of the spatial derivative and $U'(x)$ out of the velocity derivative:

$$
\begin{aligned}
    \text{Term}_1 &= -v\frac{\partial p_{eq}}{\partial x} + \frac{U'(x)}{m}\frac{\partial p_{eq}}{\partial v} \\
    &= -v\left(-\frac{U'(x)}{k_B T} p_{eq}\right) + \frac{U'(x)}{m}\left(-\frac{mv}{k_B T} p_{eq}\right) \\
    &= \left(\frac{v U'(x)}{k_B T} - \frac{v U'(x)}{k_B T}\right)p_{eq} \\
    &= 0
\end{aligned}
$$

This confirms that the conservative forces obey Liouville's Theorem and do not alter the distribution on their own.

**Part 2: The Thermostat Terms (Dissipation and Fluctuation)**
Let us look at the friction and noise terms:

$$
    \text{Term}_2 = \gamma\frac{\partial}{\partial v}(vp_{eq}) + \frac{\gamma k_BT}{m}\frac{\partial}{\partial v}\left(\frac{\partial p_{eq}}{\partial v}\right)
$$

Substitute $\frac{\partial p_{eq}}{\partial v} = -\frac{mv}{k_B T} p_{eq}$ into the diffusion term:

$$
\begin{aligned}
    \text{Term}_2 &= \gamma\frac{\partial}{\partial v}(vp_{eq}) + \frac{\gamma k_BT}{m}\frac{\partial}{\partial v}\left(-\frac{mv}{k_B T} p_{eq}\right) \\
    &= \gamma\frac{\partial}{\partial v}(vp_{eq}) - \gamma\frac{\partial}{\partial v}(v p_{eq}) \\
    &= 0
\end{aligned}
$$

The friction (damping) and diffusion (noise) terms cancel each other out exactly.

**Conclusion**

Since both parts of the equation sum to zero independently, the total time derivative is zero:

$$
    \frac{\partial p_{eq}}{\partial t} = 0
$$

This proves that the Langevin dynamics equation naturally preserves the canonical Boltzmann distribution as its stationary state.

## Conclusion

In this post, we bridged the gap between stochastic calculus and practical molecular dynamics.

We began by deriving **Itô's Lemma**, revealing how the "jagged" nature of Brownian motion introduces a necessary drift correction that standard calculus misses. From there, we established the **Fokker-Planck equation**, which translates chaotic individual trajectories into a deterministic evolution of probability density.

Finally, we applied this framework to **Langevin Dynamics**. By analyzing the interplay between friction and thermal noise, we proved mathematically that this specific SDE naturally converges to the **Boltzmann distribution**.

This rigorous foundation explains why integrators like `integrator = sd` in GROMACS or `LangevinIntegrator` in OpenMM work. The random noise isn't just interference; it is a carefully tuned mechanism that drives the simulation toward correct thermodynamic equilibrium.
