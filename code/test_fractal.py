import numpy as np
import matplotlib.pyplot as plt
from numba import njit, prange

# -----------------------------
# Physical parameters
# -----------------------------
m1 = 1.0
m2 = 1.0
l1 = 1.0
l2 = 1.0
g = 9.82

dt = 0.01
steps = 3000

# -----------------------------
# Derivatives (Numba JIT)
# -----------------------------
@njit(fastmath=True)
def derivatives(state):
    th1, w1, th2, w2 = state
    d = th1 - th2

    c = np.cos(d)
    s = np.sin(d)

    den1 = (m1 + m2) * l1 - m2 * l1 * c * c
    den2 = (l2 / l1) * den1

    dw1 = (
        m2 * l1 * w1 * w1 * s * c
        + m2 * g * np.sin(th2) * c
        + m2 * l2 * w2 * w2 * s
        - (m1 + m2) * g * np.sin(th1)
    ) / den1

    dw2 = (
        -m2 * l2 * w2 * w2 * s * c
        + (m1 + m2) * g * np.sin(th1) * c
        - (m1 + m2) * l1 * w1 * w1 * s
        - (m1 + m2) * g * np.sin(th2)
    ) / den2

    return np.array((w1, dw1, w2, dw2))


# -----------------------------
# RK4 step
# -----------------------------
@njit(fastmath=True)
def rk4_step(state):
    k1 = derivatives(state)
    k2 = derivatives(state + 0.5 * dt * k1)
    k3 = derivatives(state + 0.5 * dt * k2)
    k4 = derivatives(state + dt * k3)
    return state + dt / 6.0 * (k1 + 2*k2 + 2*k3 + k4)


# -----------------------------
# Fractal computation
# -----------------------------
@njit(parallel=True, fastmath=True)
def compute_fractal(theta1_vals, theta2_vals):
    N = len(theta1_vals)
    fractal = np.zeros((N, N))

    for i in prange(N):
        for j in range(N):

            state = np.array((
                theta1_vals[i],
                0.0,
                theta2_vals[j],
                0.0
            ))

            for _ in range(steps):
                state = rk4_step(state)

            fractal[j, i] = np.sign(state[3])

    return fractal


# -----------------------------
# Grid
# -----------------------------
N = 600
theta1_vals = np.linspace(-np.pi, np.pi, N)
theta2_vals = np.linspace(-np.pi, np.pi, N)

# -----------------------------
# Run (first call compiles)
# -----------------------------
fractal = compute_fractal(theta1_vals, theta2_vals)

# -----------------------------
# Plot
# -----------------------------
plt.figure(figsize=(8, 8))
plt.imshow(
    fractal,
    extent=[-np.pi, np.pi, -np.pi, np.pi],
    origin="lower",
    cmap="inferno"
)
plt.xlabel(r"$\theta_1$")
plt.ylabel(r"$\theta_2$")
plt.title("Double Pendulum Fractal (Numba accelerated)")
plt.colorbar(label="sign($\omega_2$)")
plt.show()
