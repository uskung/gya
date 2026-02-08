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
h = 0.00005 ## step length in seconds
t_tot= float(input("Hur många sekunder vill du simulera pendeln? "))
t0 = 0
procent_difference = float(input("Hur många procent av vinklarna vill du förändra i andra pendeln? "))/100

m_1 = 1
m_2 = 1
l_1 = 1
l_2 = 1
g = 9.82

the1_list_1 = [the1]
the2_list_1 = [the2]
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

x1 = l_1*np.sin(state[0])
x1 = l_1*np.sin(state[0])
x1pos.append(x1)
y1 = -1 * l_1 * np.cos(state[0])
y1pos.append(y1)

## calculates coordinates of mass 2
x2 = l_1*np.sin(state[0]) + l_2*np.sin(state[2])
x2pos.append(x2)
y2 = -1 * l_1*np.cos(state[0]) - l_2*np.cos(state[2])
y2pos.append(y2)


## first pendulum
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

    the1_list_1.append(state[0])
    the2_list_1.append(state[2])

    t0 += h

## second pendulum
the1_list_2 = [the1 + procent_difference*the1]
the2_list_2 = [the2 + procent_difference*the2]
X1pos = []
Y1pos = []
X2pos = []
Y2pos = []

t0 = 0
state = np.array([the1 + procent_difference*the1, ome1 + procent_difference*ome1, the2 + procent_difference*the2, ome2 + procent_difference*ome2])
while t0 < t_tot + h:
    X1 = l_1*np.sin(state[0])
    X1pos.append(X1)
    Y1 = -1 * l_1 * np.cos(state[0])
    Y1pos.append(Y1)

    ## calculates coordinates of mass 2
    X2 = l_1*np.sin(state[0]) + l_2*np.sin(state[2])
    X2pos.append(X2)
    Y2 = -1 * l_1*np.cos(state[0]) - l_2*np.cos(state[2])
    Y2pos.append(Y2)
    ## calculates everything for ome1

    K1 = derivative(state)
    K2 = derivative(state + h/2 * K1)
    K3 = derivative(state + h/2 * K2)
    K4 = derivative(state + h * K3)

    state = state + h/6 * (K1 + 2*K2 + 2*K3 + K4)

    the1_list_2.append(state[0])
    the2_list_2.append(state[2])

    t0 += h

steps_last2s = int(2/h) ## calculates how many steps are in the last 5 seconds of the simulation, used to plot only the last 2 seconds of the pendulum movement

# Draw path for last 2 seconds for mass 5 for pendulum 1
plt.plot(x2pos[-steps_last2s:], y2pos[-steps_last2s:], color='blue', ls = ':', lw = '1')

# Draw path for last 2 seconds for mass 2 for pendulum 2
plt.plot(X2pos[-steps_last2s:], Y2pos[-steps_last2s:], color='red', ls = ':', lw = '1')

### pendulum 1
plt.plot(x1pos[-1], y1pos[-1], color='blue', markersize=8, marker = 'o')
plt.plot(x2pos[-1], y2pos[-1], color='blue', markersize=8, marker = 'o',label='Pendel 1')
plt.plot([0,x1pos[-1]], [0,y1pos[-1]], color='blue')
plt.plot([x1pos[-1], x2pos[-1]], [y1pos[-1], y2pos[-1]], color='blue')

### pendulum 2, has deviation
plt.plot(X1pos[-1], Y1pos[-1], color='red', markersize=8, marker = 'o')
plt.plot(X2pos[-1], Y2pos[-1], color='red', markersize=8, marker = 'o', label=f'Pendel 2')
plt.plot([0,X1pos[-1]], [0,Y1pos[-1]], color='red')
plt.plot([X1pos[-1], X2pos[-1]], [Y1pos[-1], Y2pos[-1]], color='red')

plt.legend()
ax = plt.gca()
ax.set_xlim([-2.5, 2.5])
ax.set_ylim([-2.5, 2.5])
ax.set_aspect('equal', adjustable='box')
plt.grid()
plt.xlabel('$x$-position (m)')
plt.ylabel('$y$-position (m)') 

#plt.savefig(f'multiple_pendulums_plots/multiple_pendulums_plot_at_t={t_tot}s_the1={the1}_the2={the2}_procent_diff={procent_difference}.png', dpi=300)

plt.show()