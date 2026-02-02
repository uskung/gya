import numpy as np

# --- parameters ---
the1 = float(eval(input("Ange startvinkeln för theta_1: ")))
the2 = float(eval(input("Ange startvinkeln för theta_2: ")))
ome1 = 0.0
ome2 = 0.0

h = 0.00005
t_tot = float(input("Hur många sekunder vill du simulera pendeln? "))

m_1 = m_2 = 1.0
l_1 = l_2 = 1.0
g = 9.82

# --- initial states ---
state = np.array([the1, ome1, the2, ome2])
delta0 = 1e-8
state_p = state + np.array([delta0, 0, 0, 0])

# --- dynamics ---
def derivative(state):
    the1, ome1, the2, ome2 = state
    dtheta = the1 - the2

    alpha = (m_1 + m_2)*l_1
    beta = m_2*l_2*np.cos(dtheta)
    gamma = m_2*l_1*np.cos(dtheta)
    delta = m_2*l_2

    epsilon = -m_2*l_2*ome2**2*np.sin(dtheta) - (m_1+m_2)*g*np.sin(the1)
    zeta = m_2*l_1*ome1**2*np.sin(dtheta) - m_2*g*np.sin(the2)

    dome1 = (delta*epsilon - beta*zeta)/(alpha*delta - beta*gamma)
    dome2 = (alpha*zeta - gamma*epsilon)/(alpha*delta - beta*gamma)

    return np.array([ome1, dome1, ome2, dome2])

def rk4_step(state):
    K1 = derivative(state)
    K2 = derivative(state + h/2 * K1)
    K3 = derivative(state + h/2 * K2)
    K4 = derivative(state + h * K3)
    return state + h/6 * (K1 + 2*K2 + 2*K3 + K4)

def phase_space_distance(a, b):
    d = a - b
    d[0] = (d[0] + np.pi) % (2*np.pi) - np.pi
    d[2] = (d[2] + np.pi) % (2*np.pi) - np.pi
    return np.linalg.norm(d)

# --- Lyapunov calculation ---
lyap_sum = 0.0
t = 0.0
renorm_steps = 20
step = 0

# optional: ignore transient
t_transient = 5.0

while t < t_tot:
    state = rk4_step(state)
    state_p = rk4_step(state_p)

    if step % renorm_steps == 0:
        d = phase_space_distance(state_p, state)

        if d > 0 and t > t_transient:
            lyap_sum += np.log(d / delta0)

        # renormalize
        state_p = state + delta0 * (state_p - state) / d

    t += h
    step += 1

lyapunov = lyap_sum / (t_tot - t_transient)
print("Largest Lyapunov exponent:", lyapunov)
