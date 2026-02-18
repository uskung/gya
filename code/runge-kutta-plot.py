### This code generates a still image of a double pendulum and its corresponding path of the second mass.
# The double pendulum is simulated with the RK4 method. 

import matplotlib.pyplot as plt 
import numpy as np

################
# Initial conditions
m_1 = 1 # mass 1[kg]
m_2 = 1 # mass 2[kg]
l_1 = 1 # length of rod 1 [m]
l_2 = 1 # length of rod 2 [m]
g = 9.82 # gravitational acceleration [m/s^2]

ome1 = 0 # angular velocity of the1 [rad/s]
ome2 = 0 # angular velocity of the2 [rad/s]
################ 

##############
# User input for starting angles
the1 = float(eval(input("Ange startvinkeln for theta_1: "))) # angle 1 [rad]
the2 = float(eval(input("Ange startvinkeln for theta_2: "))) # angle 2 [rad]

# User input for simulation time
t_tot= float(input("Hur manga sekunder vill du simulera pendeln? "))
##############

##############
# Conditions for simulation
h = 0.00005 # time step [s]
t0 = 0 # initial time [s]
###############

################
# Define lists that will be used in simulation to store data
the1_list = [the1] # list with the1 angles
the2_list = [the2] # list with the2 angles
x1pos = [] # list with x1 positions
y1pos = [] # list with y1 positions
x2pos = [] # list with x2 positions
y2pos = [] # list with y2 positions
################

# Define state vector, which describes the pendulums state
state = np.array([the1, ome1, the2, ome2])

### Define function that calculates next domega1 and domega2, returns an array which is the time derivative of the 'state' array 
def derivative(state):
    ### calculates coefficients for RK4-method
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

    return np.array([ome1, domega1, ome2, domega2]) ## returns array which is the time derivative of the 'state' array

while t0 < t_tot + h:
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

    ## calculates K1 through K4 according to RK4 method
    K1 = derivative(state)
    K2 = derivative(state + h/2 * K1)
    K3 = derivative(state + h/2 * K2)
    K4 = derivative(state + h * K3)

    ## define next state
    state = state + h/6 * (K1 + 2*K2 + 2*K3 + K4) 

    the1_list.append(state[0])
    the2_list.append(state[2])

    ## for each iteration, add time step to current time
    t0 += h 

###############
# Makes font to Metafont (same as in LaTeX)
plt.rc('font', size = 11, family='serif')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern')
################

################
# Plots axis, axis-labels and grid
ax = plt.gca()
ax.set_xlim([-2.5,2.5])
ax.set_ylim([-2.5,2.5])
ax.set_aspect('equal', adjustable='box')

plt.xlabel(' $x$-position [m]')
plt.ylabel('$y$-position [m]')


plt.grid()
################

################
# Plots the pendulum
plt.plot(x1pos[-1], y1pos[-1], 'o', markersize=15, color='red') ## plots mass 1
plt.plot(x2pos[-1], y2pos[-1], 'o', markersize=15, color='red') ## plots mass 2
plt.plot([0,x1pos[-1]], [0,y1pos[-1]], color='blue') # plots rod 1
plt.plot([x1pos[-1], x2pos[-1]], [y1pos[-1], y2pos[-1]], color='blue') # plots rod 2
plt.plot(x2pos, y2pos, color='red', ls=':', label='Fardvag av massa 2') # plots path of mass 2
################

# Saves png to directory, commented away currently as to make it easer to run independently
#plt.savefig(f'runge-kutta_plots/runge-kutta_plot_the1={the1}_the2={the2}_at_{t_tot}s.png', dpi=300)

# Plots figure
plt.show()