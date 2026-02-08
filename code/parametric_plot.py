import matplotlib.pyplot as plt 
import numpy as np
from matplotlib.animation import FuncAnimation

plt.rc('font', size = 11, family='serif')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern')

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
# x1pos = []
# y1pos = []
# x2pos = []
# y2pos = []

state = np.array([the1, ome1, the2, ome2])

def derivative(state):
    the1, ome1, the2, ome2 = state
    dtheta = the1 - the2
    alpha = (m_1 + m_2)*l_1
    beta = m_2*l_2*np.cos(dtheta)
    gamma = m_2*l_1*np.cos(dtheta)
    delta = m_2*l_2
    epsilon = -1*m_2*l_2 * ome2**2 * np.sin(dtheta) - (m_1 + m_2)*g*np.sin(the1)
    zeta = m_2*l_2 * ome1**2 * np.sin(dtheta) - m_2*g*np.sin(the2)

    domega1 = (delta*epsilon - beta*zeta)/(alpha*delta - beta*gamma) # calculates new domega1
    domega2 = (alpha*zeta - gamma*epsilon)/(alpha*delta - beta*gamma) # calculates new domega2

    return np.array([ome1, domega1, ome2, domega2])

while t0 < t_tot + h:
    #calculate coordinates of pendulum
    ## calculates coordinates of mass 1
    x1 = l_1*np.sin(state[0])
    # x1pos.append(x1)
    # y1 = -1 * l_1 * np.cos(state[0])
    # y1pos.append(y1)

    # ## calculates coordinates of mass 2
    # x2 = l_1*np.sin(state[0]) + l_2*np.sin(state[2])
    # x2pos.append(x2)
    # y2 = -1 * l_1*np.cos(state[0]) - l_2*np.cos(state[2])
    # y2pos.append(y2)

    K1 = derivative(state)
    K2 = derivative(state + h/2 * K1)
    K3 = derivative(state + h/2 * K2)
    K4 = derivative(state + h * K3)

    state = state + h/6 * (K1 + 2*K2 + 2*K3 + K4)

    the1_list.append(state[0])
    the2_list.append(state[2])

    t0 += h

# fig, axis = plt.subplots()
# animated_theta = axis.plot([], [],'o', markersize=10, color='blue')[0]
# animated_path = axis.plot([], [], color='red')[0]

# axis.set_xlim([-3,3])
# axis.set_ylim([-3,3])
# axis.set_title('Animering av dubbelpendel - Parametric plot')

# plt.grid()

# frames=round((t_tot/25)*10**3)
# animation_const = len(the1_list)/frames
# path_splice_limit = 500

# def update_data(frame):
#     animated_theta.set_data([the1_list[round(frame*animation_const)]], [the2_list[round(frame*animation_const)]])
#     #animated_path.set_data(the1_list[max(0,round(frame*animation_const)-path_splice_limit):round(frame*animation_const)], the2_list[max(0,round(frame*animation_const)-path_splice_limit):round(frame*animation_const)])
#     animated_path.set_data(the1_list[:round(frame*animation_const)], the2_list[:round(frame*animation_const)])
#     return (animated_theta,)

# animation = FuncAnimation(
#     fig=fig,
#     func=update_data,
#     frames=frames,
#     interval=25,
# )

xaxis_limit = abs(np.array(the1_list)).max() + 0.01 * abs(np.array(the1_list)).max()
yaxis_limit = abs(np.array(the2_list)).max() + 0.01 * abs(np.array(the2_list)).max()

ax = plt.gca()
ax.set_xlim([-xaxis_limit, xaxis_limit])
ax.set_ylim([-yaxis_limit, yaxis_limit])
#ax.set_title(f'Parametric plot - the1={the1}, the2={the2}, t={t_tot}s')
#ax.set_aspect('equal', adjustable='box')
plt.plot(the1_list,the2_list, color='red')
plt.xlabel(' $\\theta_1$ [rad]')
plt.ylabel('$\\theta_2$ [rad]')


plt.grid(True)

plt.savefig(f'parametric_plots/parametric_plot_the1={the1}_the2={the2}_t={t_tot}s.png', dpi=300)

plt.show()