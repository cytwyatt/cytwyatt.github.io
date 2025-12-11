---
layout: post
title: "Langevin Sampling and Score Matching"
description: "We often treat the 'normalization constant' in deep learning as a mathematical headache. But by looking at it through the lens of physics, we find a solution. This post breaks down how Langevin dynamics allows us to bypass the intractable integral and generate data using the 'score function'—physics' answer to generative modeling."
tags: Math MD DeepLearning
categories: Academic
date: 2025-12-11
---

This post demonstrates another important application of Langevin dynamics, particularly within the context of diffusion models.

The core idea is that even if we do not know the exact Boltzmann distribution, we can still use Langevin dynamics to sample from it by utilizing the potential (known as the force field in molecular dynamics).

We see a striking parallel in deep learning. Given training data, the true data distribution is unknown and difficult to model explicitly. This is due to the existence of the normalization factor—a high-dimensional integral known as the *partition function* in statistical mechanics.

However, using the analogy of molecular dynamics, we can instead model the "force," which in statistics is known as the **score function**. This allows us to bypass modeling the true distribution explicitly. With knowledge of the score function alone, we can draw correct samples from the true distribution, forming the basis of modern generative models.

Prerequisites: [Langevin Dynamics: From Wiener Process to Stochastic Differential Equations](/blog/2025/langevin-dynamics1/)

---

## Introduction

In generative modeling, we begin with a dataset $\{x_i\}_{i=1}^n$ originating from an underlying distribution $p(x)$. The goal is to train a neural network $p_\theta(x)$, parameterized by $\theta$, to approximate this true distribution. If the network successfully learns to mimic $p(x)$, we can then draw samples from $p_\theta(x)$ to produce new, artificially generated data.

While this approach sounds promising, effectively modeling the true distribution directly is often impossible. To understand why, consider the standard Maximum Likelihood Estimation (MLE) objective, where we aim to maximize the log-likelihood of the observed data with respect to $\theta$:

$$
    \hat{\theta} = \mathrm{argmax}_{\theta}\sum_{i=1}^n \log p_\theta(x_i)
$$

The difficulty arises when we define the probability density $p_\theta(x)$. To ensure it is a valid probability distribution, it must sum (or integrate) to one. Therefore, the model is typically expressed as an unnormalized energy function $f_\theta(x)$ divided by a normalization constant:

$$
    p_\theta(x) = \frac{e^{-f_\theta(x)}}{Z_\theta}
$$

Here, $Z_\theta$ is the partition function, which ensures the distribution is normalized. If we substitute this back into our MLE objective, the log-likelihood becomes:

$$
    \log p_\theta(x) = -f_\theta(x) - \log Z_\theta
$$

This reveals the core problem: to maximize the likelihood, we must compute $\log Z_\theta$. However, $Z_\theta$ is defined as the integral over the entire data space:

$$
    Z_\theta = \int e^{-f_\theta(x)} dx
$$

This integral extends over a high-dimensional input space, presenting a computational challenge identical to calculating the partition function in statistical mechanics. Since we cannot simply sum over every possible configuration to find the normalization constant, direct likelihood maximization becomes intractable.

Fortunately, we can bypass this hurdle by taking the gradient of the log-likelihood with respect to $x$. Because the partition function $Z_\theta$ does not depend on $x$, it vanishes during differentiation:

$$
    \nabla_x \log p_\theta(x) = -\nabla_x f_{\theta}(x)
$$

This quantity—the gradient of the log-likelihood—is known as the **score function**. It is directly analogous to the *force* in molecular dynamics. As we established earlier, simply knowing this "force" allows us to sample from the true distribution using Langevin dynamics, without ever needing to compute the intractable normalization constant.

## Langevin Sampling

Once we have the score function, we can simulate the following Stochastic Differential Equation (SDE):

$$
    dx_t = \frac{1}{2}\nabla_x\log p_\theta(x)dt + dW_t
$$

This is effectively the **overdamped** version of the Langevin equation, where the inertia term ($m dv_t$) is negligible and approaches zero. To demonstrate that this SDE indeed samples from the true distribution in the long-term limit ($t \to \infty$), we analyze the evolution of the probability density $p_t(x)$ using the **Fokker-Planck equation**.

$$
    \frac{\partial p_t(x)}{\partial t} = -\nabla \cdot \left( \frac{1}{2} (\nabla \log p_\theta(x)) p_t(x) \right) + \frac{1}{2}\Delta p_t(x)
$$

We are interested in the stationary distribution, where the probability density stops changing over time ($\frac{\partial p_t}{\partial t} = 0$). Let us test if the true distribution $p_\theta(x)$ is the solution. Substituting $p_t(x) = p_\theta(x)$ into the right-hand side:

$$
    \begin{aligned}
    \text{RHS} &= -\nabla \cdot \left( \frac{1}{2} \frac{\nabla p_\theta(x)}{p_\theta(x)} p_\theta(x) \right) + \frac{1}{2}\nabla^2 p_\theta(x) \\
               &= -\frac{1}{2} \nabla \cdot (\nabla p_\theta(x)) + \frac{1}{2}\nabla^2 p_\theta(x) \\
               &= -\frac{1}{2} \nabla^2 p_\theta(x) + \frac{1}{2}\nabla^2 p_\theta(x) \\
               &= 0
    \end{aligned}
$$

Since the time derivative is zero, $p_\theta(x)$ is indeed the stationary distribution. This proves that by following the "force" (the score function) and adding noise, the particle will eventually settle into the correct data distribution.

## Score Matching

We have now established that instead of struggling to model the probability distribution directly, we can focus on modeling the **score function**. By doing so, we can employ Langevin dynamics to sample from the true underlying distribution. The technique used to learn this score function is known as **Score Matching**.

Intuitively, we want to train a neural network $s_\theta(x)$ to approximate the true score $\nabla_x \log p_{data}(x)$. A natural objective would be to minimize the squared difference between the model's score and the true score, averaged over our data. This quantity is known as the **Fisher Divergence**:

$$
    \mathcal{L}(\theta) = \mathbb{E}_{x \sim p_{data}} \left[ \frac{1}{2} \| s_\theta(x) - \nabla_x \log p_{data}(x) \|^2_2 \right]
$$

However, we immediately encounter a paradox: to train the model, we need the "ground truth" target $\nabla_x \log p_{data}(x)$, but this is exactly the unknown quantity we are trying to estimate! Fortunately, there is a mathematical workaround. We can simplify the loss function using integration by parts:

We begin by expanding the squared Euclidean norm inside the expectation. This allows us to break the objective function into three distinct components:

$$
\begin{aligned}
    \mathcal{L}(\theta) &= \mathbb{E}_{x \sim p_{\mathrm{data}}} \left[ \frac{1}{2} \| s_\theta(x) - \nabla_x \log p_{\mathrm{data}}(x) \|^2_2 \right]\\
    &= \underbrace{\frac{1}{2}\mathbb{E}_{x \sim p_{\mathrm{data}}} \left[\|s_\theta(x)\|^2_2\right]}_{\text{(1)}} - \underbrace{\mathbb{E}_{x \sim p_{\mathrm{data}}} \left[s_\theta(x)^\top \nabla_x\log p_{\mathrm{data}}(x)\right]}_{\text{(2)}} + \underbrace{\frac{1}{2}\mathbb{E}_{x \sim p_{\mathrm{data}}}\left[\|\nabla_x\log p_{\mathrm{data}}(x)\|^2_2\right]}_{\text{(3)}}
\end{aligned}
$$

The first term (1) can be easily approximated via Monte Carlo sampling over a finite batch of data points. The third term (3) depends only on the true data distribution and not on our model parameters $\theta$; therefore, it acts as a constant offset and can be ignored during optimization.

The challenge lies in the second term (2). To make this tractable, we apply **integration by parts**. For clarity, let us demonstrate this in one dimension (though the logic generalizes to high dimensions using the divergence theorem):

$$
\begin{aligned}
    \mathbb{E}_{x \sim p_{\mathrm{data}}} \left[s_\theta(x)\nabla_x\log p_{\mathrm{data}}(x)\right] &= \int_{-\infty}^{\infty} s_\theta(x) \nabla_x \log p_{\mathrm{data}}(x) p_{\mathrm{data}}(x) dx \\
    &= \int_{-\infty}^{\infty} s_\theta(x) \frac{\nabla_x p_{\mathrm{data}}(x)}{p_{\mathrm{data}}(x)} p_{\mathrm{data}}(x) dx \\
    &= \int_{-\infty}^{\infty} s_\theta(x) \nabla_x p_{\mathrm{data}}(x) dx
\end{aligned}
$$

Now, applying integration by parts $\int u dv = uv - \int v du$:

$$
    = \left[ s_\theta(x) p_{\mathrm{data}}(x) \right]_{-\infty}^{\infty} - \int_{-\infty}^{\infty} \nabla_x s_\theta(x) p_{\mathrm{data}}(x) dx
$$

We assume that the data density vanishes at the boundaries of the domain (i.e., $p_{\mathrm{data}}(x) \to 0$ as $|x| \to \infty$). Consequently, the boundary term becomes zero. We are left with:

$$
\begin{aligned}
    \mathbb{E}_{x \sim p_{\mathrm{data}}} \left[s_\theta(x)\nabla_x\log p_{\mathrm{data}}(x)\right] &= - \int_{-\infty}^{\infty} \nabla_x s_\theta(x) p_{\mathrm{data}}(x) dx \\
    &= - \mathbb{E}_{x \sim p_{\mathrm{data}}} \left[ \nabla_x s_\theta(x) \right]
\end{aligned}
$$

Substituting this result back into our original loss function, we arrive at the final **Score Matching Objective**:

$$
    \mathcal{L}(\theta) = \mathbb{E}_{x \sim p_{\mathrm{data}}} \left[ \frac{1}{2} \|s_\theta(x)\|^2 + \nabla_x s_\theta(x) \right]
$$

In the general high-dimensional case, the term $\nabla_x s_\theta(x)$ becomes the **trace of the Jacobian matrix** of the score function, denoted as $\text{tr}(\nabla_x s_\theta(x))$.

## Conclusion

In summary, this post illustrates how we can circumvent the intractable partition function in generative modeling by shifting our focus from the probability density to the **score function**—analogous to the force field in molecular dynamics. By minimizing the score matching objective derived above, we can learn this field directly from data without ever computing the normalization constant. This enables us to use Langevin dynamics to guide random noise into high-fidelity data samples, essentially bridging the gap between statistical mechanics and the modern diffusion models that drive generative AI today.
