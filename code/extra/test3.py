import matplotlib.pyplot as plt 
import numpy as np
from matplotlib.animation import FuncAnimation

## code based on euler.py
## but reconverted to runge-kutta method as to increase accuracy

the1 = float(eval(input("Ange startvinkeln för theta_1: ")))
the2 = float(eval(input("Ange startvinkeln för theta_2: ")))
ome1 = 0
ome2 = 0
h = 0.00005
t_tot= float(input("Hur många sekunder vill du simulera pendeln? "))
t0 = 0
t_frames = np.linspace(0,40,501)

m_1 = 1
m_2 = 1
l_1 = 1
l_2 = 1
g = 9.82

the1_list = [the1]
the2_list = [the2]
x1pos = []
y1pos = []
x2pos = []
y2pos = []

def derivatives(state):
    the1, ome1, the2, ome2 = state
    dtheta = the1 - the2
    alpha = (m_1 + m_2) * l_1
    beta = m_2 * l_2 * np.cos(dtheta)
    gamma = m_2 * l_1 * np.cos(dtheta)
    delta = m_2 * l_2
    epsilon = -m_2 * l_2 * ome2**2 * np.sin(dtheta) - (m_1 + m_2) * g * np.sin(the1)
    zeta = m_2 * l_2 * ome1**2 * np.sin(dtheta) - m_2 * g * np.sin(the2)
    denom = alpha * delta - beta * gamma
    domega1 = (delta * epsilon - beta * zeta) / denom
    domega2 = (alpha * zeta - gamma * epsilon) / denom
    return np.array([ome1, domega1, ome2, domega2])

state = np.array([the1, ome1, the2, ome2])
t0 = 0

while t0 < t_tot + h:
    # Store positions
    x1 = l_1 * np.sin(state[0])
    y1 = -l_1 * np.cos(state[0])
    x2 = x1 + l_2 * np.sin(state[2])
    y2 = y1 - l_2 * np.cos(state[2])
    x1pos.append(x1)
    y1pos.append(y1)
    x2pos.append(x2)
    y2pos.append(y2)
    the1_list.append(state[0])
    the2_list.append(state[2])

    # Runge-Kutta 4th order step
    k1 = derivatives(state)
    k2 = derivatives(state + 0.5 * h * k1)
    k3 = derivatives(state + 0.5 * h * k2)
    k4 = derivatives(state + h * k3)
    state = state + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

    t0 += h

fig, axis = plt.subplots()
animated_l_1 = axis.plot([],[], color='blue')[0]
animated_l_2 = axis.plot([],[], color='blue')[0]
animated_m1 = axis.plot([],[], 'o', markersize=15, color='red')[0]
animated_m2 = axis.plot([],[],'o', markersize=15,color='red')[0]
animated_path_m2 = axis.plot([],[], color='red')[0]

axis.set_xlim([-2.5,2.5])
axis.set_ylim([-2.5,2.5])
axis.set_title('Animering av dubbelpendel')

plt.grid()

frames=round((t_tot/25)*10**3)
animation_const = len(x2pos)/frames

def update_data(frame):    
    animated_l_1.set_data([0,x1pos[round(frame*animation_const)]], [0, y1pos[round(frame*animation_const)]])
    animated_l_2.set_data([x1pos[round(frame*animation_const)], x2pos[round(frame*animation_const)]], [y1pos[round(frame*animation_const)], y2pos[round(frame*animation_const)]])

    animated_m1.set_data([x1pos[round(frame*animation_const)]],[y1pos[round(frame*animation_const)]])
    animated_m2.set_data([x2pos[round(frame*animation_const)]],[y2pos[round(frame*animation_const)]])

    animated_path_m2.set_data(x2pos[:round(frame*animation_const)], y2pos[:round(frame*animation_const)])
    return animated_l_1, animated_l_2,  animated_m1, animated_m2, animated_path_m2


animation = FuncAnimation(
    fig=fig,
    func=update_data,
    frames=frames,
    interval=25,
)     

plt.show()