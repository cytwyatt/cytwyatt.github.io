---
layout: post
title: "Collective Variables and Free Energy Profile"
date: 2025-11-03
description: "Understanding collective variables (CVs) and how to estimate free energy profiles from molecular simulations using kernel density estimation."
tags: Simulation Sampling MD
categories: Academic
related_posts: true
toc:
  beginning: true
---

This post is largely inspired by the materials in [3rd i-CoMSE Workshop](https://github.com/icomse/3rd_workshop_advanced_sampling). 

Prerequisite: Basic Molecular Dynamics, Thermodynamics, Probability Theory.

---

## Collective Variables

### Definition and Motivation

Consider a system of $N$ particles, $\lbrace \mathbf{x}_i \rbrace _{i=1}^N$, simulated in the NVT ensemble. The configuration space of this system has a dimensionality of $3N$, which makes it far too complex to fully monitor or visualize. To address this, we can define a carefully chosen $d$-dimensional vector

$$
\boldsymbol{\theta} = \big(\theta_1(\mathbf{x}), \theta_2(\mathbf{x}), \dots, \theta_d(\mathbf{x})\big),
$$

known as the **collective variables (CVs)**. These serve as a reduced, lower-dimensional representation obtained by projecting the full $3N$-dimensional configuration space onto a subspace that captures the essential physics of the system.

An effective CV should capture the slowest degrees of freedom and distinguish between the relevant metastable states. For instance, in protein folding studies, the **root-mean-square deviation (RMSD)** from the native structure often serves as a good CV because it reflects large-scale conformational changes.

However, identifying suitable CVs is often non-trivial, especially for systems where the slowest collective motions are not obvious. In recent years, **machine learning–based approaches** have shown promise in automatically discovering meaningful CVs. For practical examples, you can explore the [mlcolvar](https://github.com/luigibonati/mlcolvar) package, which implements some data-driven methods.

### Free Energy

We often care about the free energy of a system, as it serves as a measure of its stability. In the NVT ensemble, the Helmholtz free energy is defined as

$$
F = -\frac{1}{\beta} \ln Z = -\frac{1}{\beta} \ln \int \mathrm{e}^{-\beta U(\mathbf{x})} \, \mathrm{d}\mathbf{x}
$$

By selecting a collective variable (CV) $ \theta(\mathbf{x}) $, we effectively group together all configurations that share the same value of $ \theta $. This corresponds to the **marginal probability** along $ \theta $, which can be expressed using a Dirac delta function:

$$
p(\theta) = \frac{1}{Z} \int \mathrm{e}^{-\beta U(\mathbf{x})} \, \delta(\theta - \theta(\mathbf{x})) \, \mathrm{d}\mathbf{x} = \frac{Z_\theta}{Z}
$$

We can then define the free energy along this coordinate, often referred to as the **free energy profile** or **potential of mean force (PMF)** (these two terms are generally interchangeable) as

$$
F(\theta) = -\frac{1}{\beta} \ln Z_{\theta} = -\frac{1}{\beta} \ln p(\theta) + \text{const.} = -\frac{1}{\beta} \ln \int \mathrm{e}^{-\beta U(\mathbf{x})} \, \delta(\theta - \theta(\mathbf{x})) \, \mathrm{d}\mathbf{x} + \text{const.}
$$

The constant term is unimportant since we can choose an arbitrary reference level for the free energy. Thus, the **free energy profile** can be written as

$$
F(\theta) = -\frac{1}{\beta} \ln p(\theta) = -\frac{1}{\beta} \ln \int \mathrm{e}^{-\beta U(\mathbf{x})} \, \delta(\theta - \theta(\mathbf{x})) \, \mathrm{d}\mathbf{x} = -\frac{1}{\beta} \ln \langle \delta(\theta - \theta(\mathbf{x})) \rangle
$$

## Estimating Free Energy from Simulation

From the previous derivation, we know that estimating the free energy profile reduces to estimating the probability distribution of a given CV $p(\theta)$. A common and effective method for this is **Kernel Density Estimation (KDE)**.  

If we have sampled a set of configurations and computed their corresponding CV values $\lbrace \theta(\mathbf{x}_i) \rbrace _{i=1}^M$
, then the estimated distribution is given by

$$
p(\theta) = \frac{1}{M} \sum_{i=1}^M K_h \big(\theta - \theta(\mathbf{x}_i)\big),
$$

where $K$ is the **kernel function**, and $h$ is the **bandwidth** (a smoothing parameter). A larger bandwidth results in a smoother estimate, while a smaller one makes the distribution appear more peaked, approaching a delta function.  

There are various choices of kernels; the most common one is the **Gaussian kernel**:

$$
K_h(x) = \frac{1}{\sqrt{2\pi}h} \, \mathrm{e}^{-\frac{x^2}{2h^2}}.
$$

If the configurations are sampled from a **biased distribution**, as discussed in the [previous post](https://cytwyatt.github.io/blog/2025/enhanced-sampling1/), we can correct for this bias using a **weighted KDE**:

$$
p(\theta) = \frac{1}{M} \sum_{i=1}^M w_i \, K_h(\theta - \theta_i), \quad
w_i = \frac{\mathrm{e}^{\beta U_b(\theta(\mathbf{x}_i))}}{\sum\limits_{i=1}^M \mathrm{e}^{\beta U_b(\theta(\mathbf{x}_i))}}.
$$

Here, $U_b(\theta)$ is the applied bias potential, and the weights $w_i$ reweight the biased samples back to the unbiased ensemble.

## Conclusion

In this post, I briefly discussed the concepts of **collective variables (CVs)** and **free energy profiles**, which are central to enhanced sampling methods. In the next post, we will explore how to use **OPES (On-the-fly Probability Enhanced Sampling)** to construct the free energy profile in practice.

---

**Enhanced Sampling Series:**

- [Part 1: Reweighting Technique in Molecular Simulations](/blog/2025/enhanced-sampling1/)
- **Part 2: Collective Variables and Free Energy Profile** (current)
- [Part 3: OPES and Implementation in Toy Model](/blog/2025/enhanced-sampling3/)
