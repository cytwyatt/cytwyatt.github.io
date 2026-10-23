"""Compare Euler and Verlet energy for the unit harmonic oscillator."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def verlet(h, duration):
    """Unit mass and frequency, starting from q = 1 and p = 0."""
    time = h * np.arange(int(np.floor(duration / h)) + 1)
    q = np.empty_like(time)
    p = np.empty_like(time)
    q[0], p[0] = 1.0, 0.0
    for n in range(len(time) - 1):
        half_p = p[n] - 0.5 * h * q[n]
        q[n + 1] = q[n] + h * half_p
        p[n + 1] = half_p - 0.5 * h * q[n + 1]
    return time, q, p


def euler(h, duration):
    """Forward Euler: both updates use the old position and momentum."""
    time = h * np.arange(int(np.floor(duration / h)) + 1)
    q = np.empty_like(time)
    p = np.empty_like(time)
    q[0], p[0] = 1.0, 0.0
    for n in range(len(time) - 1):
        q[n + 1] = q[n] + h * p[n]
        p[n + 1] = p[n] - h * q[n]
    return time, q, p


def main():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.titleweight": "bold",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.labelcolor": "#334155",
            "text.color": "#182b45",
            "xtick.color": "#475569",
            "ytick.color": "#475569",
            "grid.color": "#e2e8f0",
            "grid.linewidth": 0.7,
            "legend.frameon": False,
            "savefig.facecolor": "white",
        }
    )
    blue, green, orange, gray = "#2563a6", "#008574", "#d05b36", "#64748b"
    period = 2 * np.pi
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.7), layout="constrained")
    fig.suptitle(
        "Harmonic oscillator: accumulating error versus bounded oscillations",
        fontsize=16,
        weight="bold",
    )
    for ax in axes:
        ax.grid(alpha=0.8)
        ax.set_axisbelow(True)

    ax = axes[0]
    h = 0.2
    t, q, p = euler(h, 6 * period)
    energy = q**2 + p**2  # H / H_0, with H_0 = 1/2.
    predicted = (1 + h**2) ** np.arange(len(t))
    np.testing.assert_allclose(energy, predicted, rtol=2e-13)
    ax.semilogy(t / period, energy, color=orange, lw=2, label="Forward Euler")
    t, q, p = verlet(h, 6 * period)
    ax.semilogy(t / period, q**2 + p**2, color=blue, lw=2, label="Velocity Verlet")
    ax.set(
        title=r"A. Same time step: $h=0.2$",
        xlabel="Time / physical period",
        ylabel=r"Physical energy $H_n/H_0$ (log scale)",
        xlim=(0, 6),
        ylim=(0.8, 3000),
    )
    ax.legend(loc="upper left", fontsize=10)

    ax = axes[1]
    for h, color in [(0.4, blue), (0.2, green)]:
        t, q, p = verlet(h, 6 * period)
        energy = q**2 + p**2
        invariant = p**2 + (1 - h**2 / 4) * q**2  # I_h / H_0.
        np.testing.assert_allclose(invariant, 1 - h**2 / 4, rtol=0, atol=2e-13)
        theta = 2 * np.arcsin(h / 2)
        analytic = 1 - h**2 / 4 * np.sin(np.arange(len(t)) * theta) ** 2
        np.testing.assert_allclose(energy, analytic, rtol=0, atol=2e-13)
        ax.plot(t / period, energy, color=color, lw=1.7, label=rf"$H_n/H_0$, $h={h}$")
        ax.plot(t / period, invariant, color=color, ls="--", lw=1.4, label=rf"$I_h/H_0$, $h={h}$")
    ax.axhline(1, color=gray, ls=":", lw=1)
    ax.set(
        title="B. Verlet: physical and modified energy",
        xlabel="Time / physical period",
        ylabel="Energy / initial physical energy",
        xlim=(0, 6),
        ylim=(0.953, 1.026),
        yticks=[0.96, 0.98, 1.0],
    )
    ax.legend(loc="upper left", ncol=2, fontsize=9, columnspacing=1.5, handlelength=2.0)

    output = Path(__file__).with_name("harmonic-oscillator.png")
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
