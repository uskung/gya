import matplotlib.pyplot as plt 
import numpy as np
from matplotlib.animation import FuncAnimation

plt.rc('font', size = 11, family='serif')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern')

the1 = float(eval(input("Ange startvinkeln för theta_1: ")))
the2 = float(eval(input("Ange startvinkeln för theta_2: ")))
ome1 = 0
ome2 = 0
h = 0.00005
t_tot= float(input("Hur många sekunder vill du simulera pendeln? "))
t0 = 0

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
    x1pos.append(x1)
    y1 = -1 * l_1 * np.cos(state[0])
    y1pos.append(y1)

    ## calculates coordinates of mass 2
    x2 = l_1*np.sin(state[0]) + l_2*np.sin(state[2])
    x2pos.append(x2)
    y2 = -1 * l_1*np.cos(state[0]) - l_2*np.cos(state[2])
    y2pos.append(y2)

    K1 = derivative(state)
    K2 = derivative(state + h/2 * K1)
    K3 = derivative(state + h/2 * K2)
    K4 = derivative(state + h * K3)

    state = state + h/6 * (K1 + 2*K2 + 2*K3 + K4)

    the1_list.append(state[0])
    the2_list.append(state[2])

    t0 += h

ax = plt.gca()
ax.set_xlim([-2.5,2.5])
ax.set_ylim([-2.5,2.5])
#ax.set_title('Animering av dubbelpendel - RK4 - t=0s')

plt.grid()

plt.plot(x2pos, y2pos, color='red', ls=(0,(1,10)), label='Path of mass 2')
plt.plot([0,x1pos[-1]], [0,y1pos[-1]], color='blue')
plt.plot([x1pos[-1], x2pos[-1]], [y1pos[-1], y2pos[-1]], color='blue')
plt.plot(x1pos[-1], y1pos[-1], 'o', markersize=15, color='red')
plt.plot(x2pos[-1], y2pos[-1], 'o', markersize=15, color='red')

plt.xlabel(' $x$-position [m]')
plt.ylabel('$y$-position [m]')

ax = plt.gca()
ax.set_xlim([-2.5, 2.5])
ax.set_ylim([-2.5, 2.5])
ax.set_aspect('equal', adjustable='box')
# plt.plot(xpos,ypos)

plt.savefig(f'runge-kutta_plots/runge-kutta_plot_the1={the1}_the2={the2}_at_{t_tot}s.png', dpi=300)

plt.show()


