---
layout: post
title: "Reweighting Technique in Molecular Simulations: From Umbrella Sampling to MBAR"
date: 2025-10-01
description: Enhanced Sampling 1
tags: Simulation Math Sampling MD
categories: Academic
related_posts: false
---

This post is largely inspired by the materials in [3rd i-CoMSE Workshop](https://github.com/icomse/3rd_workshop_advanced_sampling). 

Prerequisite: Basic Molecular Dynamics, Thermodynamics, Probability Theory.

---

#### Boltzmann Distribution
Starting from the second law of thermodynamics, we can derive that, at a given temperature, the probability of finding a specific microscopic configuration $x_i$ follows the **Boltzmann distribution** (sometimes called the Gibbs distribution)

$$
    p_i \propto e^{-\beta E_i} = \frac{1}{Z} \mathrm{e}^{-\beta E_i}
$$

where $\beta = 1/(k_B T)$ is inversely proportional to the temperature of the system, and have the unit of inverse energy (it is useful to remember that at room temperature, $k_B T \approx 2.5 \,\mathrm{kJ/mol}$).

The normalization factor $Z$, called the **partition function**, is defined as

$$
    Z = \sum_{i} \mathrm{e}^{-\beta E_i}
$$

There are many interpretations of the partition function. One useful way to think about it is as a **weighted count of all microscopic configurations**. To see this, we can insert a magic one in the expression (we will use this trick many times in the post)

$$
    Z = \sum_{i} 1 \cdot \mathrm{e}^{-\beta E_i}
$$

If the exponential factor were always $1$, then $Z$ would simply count the number of accessible configurations. But the exponential factor, also called the **Boltzmann factor** or **Boltzmann weight**, assigns more importance to lower-energy states. Thus, $Z$ is not just a count, but a weighted sum that reflects how nature prefers lower-energy configurations at finite temperature.

From this formula we also see a key feature of the Boltzmann distribution: all microstates with the same energy are equally probable. In other words, the probability depends only on the energy, not on any other label of the configuration.

#### Ensemble
In statistical mechanics, an ensemble is a collection of hypothetical copies of a system that all share the same macroscopic properties (such as temperature, pressure, volume, and chemical potential) but may differ in their microscopic details (such as the positions and velocities of individual particles). The most common ensemble is the so-called canonical ensemble (also $NVT$ ensemble), where the number of particles $N$, volume $V$, and temperature $T$ are fixed. The probability distribution is given by the Boltzmann distribution.

In a molecular system, the microscopic configuration is determined by the continuous positions $\mathbf{x}$ and momenta $\mathbf{p}$ (note that if the system has $N$ particles, the dimension of both $\mathbf{x}$ and $\mathbf{p}$ is $3N$). Thus we can rewrite the Boltzmann distribution as

$$
    p(\mathbf{x}, \mathbf{p}) = 
    \frac{\mathrm{e}^{-\beta \left(U(\mathbf{x}) + K(\mathbf{p})\right)}}
         {\iint \mathrm{e}^{-\beta \left(U(\mathbf{x}) + K(\mathbf{p})\right)} \, \mathrm{d}\mathbf{x} \, \mathrm{d}\mathbf{p}}
$$

Most of the time, we are only interested in the configurational distribution, so we can integrate out the momenta. This gives

$$
    p(\mathbf{x}) = 
    \frac{\int \mathrm{e}^{-\beta \left(U(\mathbf{x}) + K(\mathbf{p})\right)} \, \mathrm{d}\mathbf{p}}
         {\iint \mathrm{e}^{-\beta \left(U(\mathbf{x}) + K(\mathbf{p})\right)} \, \mathrm{d}\mathbf{x} \, \mathrm{d}\mathbf{p}}
    = \frac{\mathrm{e}^{-\beta U(\mathbf{x})}}
           {\int \mathrm{e}^{-\beta U(\mathbf{x})} \, \mathrm{d}\mathbf{x}}
    = \frac{1}{Z} \, \mathrm{e}^{-\beta U(\mathbf{x})}
$$

where $Z$ is the configurational partition function. Note that $Z$ is a muitl-dimensional intergral ($3N$) which is generally impossible to evaluate.

#### Expectation

Often, we are interested in how to measure a physical observable $O$ in a given ensemble. From probability theory, we can calculate its expectation value as

$$
    \langle O \rangle_p = \int O(\mathbf{x}) \, p(\mathbf{x}) \, \mathrm{d}\mathbf{x} 
    = \frac{\int O(\mathbf{x}) \, \mathrm{e}^{-\beta U(\mathbf{x})} \, \mathrm{d}\mathbf{x}}
           {\int \mathrm{e}^{-\beta U(\mathbf{x})} \, \mathrm{d}\mathbf{x}}
$$

However, for most systems beyond the trivial cases, evaluating this integral exactly is impossible, and we must seek approximations. Suppose that we can somehow sample from the Boltzmann distribution (we will return to sampling methods later). In that case, we can use **Monte Carlo integration** to approximate the expectation value (I will use equal sign in the rest of the post, but keep on mind this is an approximation and the error decreases as the increase of samples. Meanwhile the convergence is only guaranteed for ergodic system, which is not always true):

$$
    \langle O \rangle_p = \int O(\mathbf{x}) \, p(\mathbf{x}) \, \mathrm{d}\mathbf{x} 
    \approx \frac{1}{N} \sum_{n=1}^N O(\mathbf{x}_n),\quad\mathbf{x}_n\sim p(\mathbf{x})
$$

where each $\mathbf{x}_n$ is sampled from $p(\mathbf{x})$. It is important to note that there is no explicit Boltzmann factor in the Monte Carlo estimator, because it is already implicitly included in the sampling procedure.

Everything looks good, but the problem is that we cannot efficiently sample from an arbitrary distribution, especially for the complex, high-dimensional Boltzmann distribution. This is where the method of **importance sampling** comes in. It allows us to estimate observables under any distribution by using an easier-to-sample proposal distribution.

Suppose we want to measure an observable $O$ under distribution $p(\mathbf{x})$, and we can find another proposal distribution $q(\mathbf{x})$. Then

$$
    \langle O \rangle_p = \int O(\mathbf{x}) \, p(\mathbf{x}) \, \mathrm{d}\mathbf{x} 
    = \int O(\mathbf{x}) \, p(\mathbf{x}) \, \frac{q(\mathbf{x})}{q(\mathbf{x})} \, \mathrm{d}\mathbf{x}
    = \int O(\mathbf{x}) \, \frac{p(\mathbf{x})}{q(\mathbf{x})} \, q(\mathbf{x}) \, \mathrm{d}\mathbf{x}
    = \frac{1}{N} \sum_{n=1}^N O(\mathbf{x}_n) \, \frac{p(\mathbf{x}_n)}{q(\mathbf{x}_n)}, \quad \mathbf{x}_n \sim q(\mathbf{x})
$$

This leads to an important identity

$$
    \langle O \rangle _p
    = \frac{1}{N} \sum_{n=1}^N O(\mathbf{x}_n), \quad \mathbf{x}_n \sim p(\mathbf{x})
    \;\;\;\; = \frac{1}{N} \sum_{n=1}^N O(\mathbf{x}_n) \, \frac{p(\mathbf{x}_n)}{q(\mathbf{x}_n)}, \quad \mathbf{x}_n \sim q(\mathbf{x})
$$

or in a more compact notation

$$
    \langle O \rangle _p = \left\langle O \frac{p}{q}\right\rangle _q
$$

This means that if we want to measure an observable $O$ under distribution $p$, instead of sampling directly from $p$, we can sample from another distribution $q$. The expectation value of $O$ is then obtained as the average of $O(\mathbf{x}_n)$ weighted by the factor $\tfrac{p(\mathbf{x}_n)}{q(\mathbf{x}_n)}$. 

This is a very powerful technique, and it is valid for any observable $O$ and for any pair of distributions $p$ and $q$. The effectiveness of the method depends on how well $q$ approximates $p$. If the ratio $\tfrac{p(\mathbf{x}_n)}{q(\mathbf{x}_n)}$ fluctuates too much, the estimator will become very noisy. Of course, finding a suitable $q$ that both approximates $p$ well and is easy to sample from is not trivial — but before worrying about that, let’s first look at some examples of how to use it.

#### Removing Bias

In many enhanced sampling methods for example the **umbrella sampling**, we apply a bias potential to help the system overcome energy barriers and explore the full configurational space. The problem is how to remove the bias and recover the true expectation value of an observable. This is where importance sampling comes in.

Suppose we want to measure an observable $O$ under the distribution $p(\mathbf{x})$, which follows the Boltzmann distribution

$$
    p(\mathbf{x}) = \frac{\mathrm{e}^{-\beta U(\mathbf{x})}}{Z_p},
$$

where $Z_p = \int \mathrm{e}^{-\beta U(\mathbf{x})} \mathrm{d}\mathbf{x}$ is the partition function. But in simulation we apply a bias $U_b(\mathbf{x})$ to the original system, which makes us sample from another Boltzmann distribution $q(\mathbf{x})$:

$$
    q(\mathbf{x}) = \frac{\mathrm{e}^{-\beta (U(\mathbf{x}) + U_b(\mathbf{x}))}}{Z_q},
$$

where $Z_q = \int \mathrm{e}^{-\beta (U(\mathbf{x})+U_b(\mathbf{x}))} \mathrm{d}\mathbf{x}$. According to importance sampling, the expectation value of $O$ under the original distribution can be estimated from biased sampling:

$$
    \langle O \rangle_p = \frac{1}{N} \sum_{n=1}^N O(\mathbf{x}_n) \frac{p(\mathbf{x}_n)}{q(\mathbf{x}_n)}
    = \frac{1}{N} \sum_{n=1}^N O(\mathbf{x}_n) \, \mathrm{e}^{\beta U_b(\mathbf{x}_n)} \, \frac{Z_q}{Z_p}.
$$

The problem now is how to evaluate the partition function ratio $\tfrac{Z_q}{Z_p}$. The trick is to set $O = 1$, in which case the equation becomes

$$
    1 = \frac{1}{N} \sum_{n=1}^N \mathrm{e}^{\beta U_b(\mathbf{x}_n)} \, \frac{Z_q}{Z_p},
$$

so that

$$
    \frac{Z_p}{Z_q} = \frac{1}{N} \sum_{n=1}^N \mathrm{e}^{\beta U_b(\mathbf{x}_n)}.
$$

Recall the thermodynamic relation between free energy and partition function:

$$
    A = -k_B T \ln Z.
$$

The above identity can be written as

$$
    \mathrm{e}^{-\beta (A_p - A_q)} = \langle \mathrm{e}^{-\beta (U_p - U_q)} \rangle_q,
$$

which is the **Zwanzig equation**. It provides a way to estimate the free energy difference between two states and is the foundation of Free Energy Perturbation (FEP).

Putting this back, the expectation of $O$ under the unbiased distribution becomes

$$
    \langle O \rangle_p =
    \frac{1}{N} \sum_{n=1}^N \frac{O(\mathbf{x}_n) \, \mathrm{e}^{\beta U_b(\mathbf{x}_n)}}{\frac{1}{N}\sum\limits_{n=1}^N \mathrm{e}^{\beta U_b(\mathbf{x}_n)}}
    = \sum\limits_{n=1}^N \frac{O(\mathbf{x}_n) \, \mathrm{e}^{\beta U_b(\mathbf{x}_n)}}{\sum\limits_{n=1}^N \mathrm{e}^{\beta U_b(\mathbf{x}_n)}}
$$

This is nothing but a **weighted average** of the sampled values:

$$
    \langle O \rangle_p = \sum_{n=1}^N O(\mathbf{x}_n) \, w_n,
$$

where the weights are normalized:

$$
    w_n = \frac{\mathrm{e}^{\beta U_b(\mathbf{x}_n)}}{\sum\limits_{n=1}^N \mathrm{e}^{\beta U_b(\mathbf{x}_n)}}\quad \sum_{n=1}^N w_n = 1.
$$

It can be also written as a ratio of two expectation:

$$
    \langle O \rangle_p
    = \sum\limits_{n=1}^N \frac{O(\mathbf{x}_n) \, \mathrm{e}^{\beta U_b(\mathbf{x}_n)}}{\sum\limits_{n=1}^N \mathrm{e}^{\beta U_b(\mathbf{x}_n)}}
    = \frac{\langle O \mathrm{e}^{\beta U_b} \rangle _q}{\langle \mathrm{e}^{\beta U_b} \rangle _q}
$$

Always remember the $\mathbf{x}_n$ is sampled from $q$, otherwise we don't need the weight.

#### Removing Multiple Bias

In many cases, a single biased sampling is not enough. For example, in **umbrella sampling**, we pin the system at multiple positions by applying different bias potentials centered at those positions. Each simulation then explores a restricted region of configurational space. To reconstruct the unbiased distribution, we need to combine data from all these biased simulations and remove the effects of their respective bias potentials.

Suppose we have sampled from $k$ Boltzmann distributions $p_i$, each corresponding to $N_i$ samples, with the total number of samples being $N = \sum_{i=1}^k N_i$. If we combine all the samples together, the effective distribution that describes them is the **mixture distribution** $p_m$, defined as

$$
    p_m(\mathbf{x}) = \sum_{i=1}^k \frac{N_i}{N}\, p_i(\mathbf{x}).
$$

This reflects the fact that, when pooling all the data, a fraction $\tfrac{N_i}{N}$ of the total samples originates from distribution $p_i$.  

Now, how can we recover the expectation value of an observable $O$ under the unbiased distribution $p$? By analogy with the single-bias case, we have

$$
    \langle O \rangle_p \approx \frac{1}{N} \sum_{n=1}^N O(\mathbf{x}_n)\, \frac{p(\mathbf{x}_n)}{p_m(\mathbf{x}_n)}
    = \frac{1}{N} \sum_{n=1}^N O(\mathbf{x}_n)\,
      \frac{\mathrm{e}^{-\beta U(\mathbf{x}_n)} \frac{1}{Z}}
           {\sum\limits_{i=1}^k \frac{N_i}{N}\, \mathrm{e}^{-\beta U_i(\mathbf{x}_n)} \frac{1}{Z_i}}.
$$

While this looks complicated, it reduces to a weighted average of the sampled values:

$$
    \langle O \rangle_p = \sum_{n=1}^N w_n\, O(\mathbf{x}_n),
$$

with weights

$$
    w_n = \frac{\mathrm{e}^{-\beta U(\mathbf{x}_n)} \tfrac{1}{Z}}
    {\sum\limits_{i=1}^k N_i \, \mathrm{e}^{-\beta U_i(\mathbf{x}_n)} \tfrac{1}{Z_i}}.
$$

The difficulty here is that the partition functions $Z_i$ are unknown. To handle this, we introduce the unnormalized weights

$$
    \tilde{w}_n = \frac{\mathrm{e}^{-\beta U(\mathbf{x}_n)}}{\sum\limits_{i=1}^k N_i\, \mathrm{e}^{-\beta U_i(\mathbf{x}_n)} \tfrac{1}{Z_i}},
    \qquad
    w_n = \frac{\tilde{w}_n}{\sum\limits_{m=1}^N \tilde{w}_m}.
$$

The unknown partition functions $Z_i$ can themselves be determined by applying importance sampling to each biased ensemble. Setting $O=1$ and $p=p_i$ gives

$$
    1 = \frac{1}{N} \sum\limits_{n=1}^N \frac{p_i(\mathbf{x}_n)}{p_m(\mathbf{x}_n)}
    = \frac{1}{N} \sum\limits_{n=1}^N \frac{p_i(\mathbf{x}_n)}{\sum\limits_{j=1}^k \tfrac{N_j}{N}\, p_j(\mathbf{x}_n)},
$$

where $p_i(\mathbf{x}) = \tfrac{\mathrm{e}^{-\beta U_i(\mathbf{x})}}{Z_i}$. Using the thermodynamic relation $A_i = -k_B T \ln Z_i$, we can rewrite

$$
    p_i(\mathbf{x}) = \exp\!\big[\,\beta (A_i - U_i(\mathbf{x}))\,\big].
$$

Substituting back, we obtain a set of $k$ nonlinear equations for the free energies $A_i$:

$$
    \mathrm{e}^{-\beta A_i}
    = \frac{1}{N} \sum\limits_{n=1}^N \frac{\mathrm{e}^{-\beta U_i(\mathbf{x}_n)}}
    {\sum\limits_{j=1}^k \tfrac{N_j}{N}\, \exp\!\big[\beta(A_j - U_j(\mathbf{x}_n))\big]}.
$$

These equations can be solved self-consistently to determine $\{A_i\}$. In practice, this procedure is implemented in efficient algorithms such as [pymbar](https://github.com/choderalab/pymbar).  

Finally, the unbiased expectation value is recovered as

$$
    \langle O \rangle_p = \sum_{n=1}^N w_n\, O(\mathbf{x}_n),
$$

with normalized weights

$$
    \tilde{w}_n = \frac{\mathrm{e}^{-\beta U(\mathbf{x}_n)}}
    {\sum\limits_{i=1}^k N_i\, \exp\!\big[\beta(A_i - U_i(\mathbf{x}_n))\big]},
    \qquad
    w_n = \frac{\tilde{w}_n}{\sum\limits_{m=1}^N \tilde{w}_m}.
$$

This is the central result of **MBAR (Multistate Bennett Acceptance Ratio)**, which provides the statistically optimal, minimum-variance estimate of observables and free energies when combining multiple biased simulations.

#### Conclusion

In this post, we explored the theoretical background of reweighting techniques and importance sampling in molecular simulations. These methods provide the foundation for correcting biased ensembles and recovering unbiased thermodynamic properties, even when direct sampling is impractical. 

In the next post, I will present concrete examples and applications to show how these ideas are implemented in practice.
