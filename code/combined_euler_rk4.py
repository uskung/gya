import matplotlib.pyplot as plt 
import numpy as np
from matplotlib.animation import FuncAnimation

the1_gbl = float(eval(input("Ange startvinkeln för theta_1: ")))
the2_gbl = float(eval(input("Ange startvinkeln för theta_2: ")))
ome1_gbl = 0
ome2_gbl = 0
h = 0.00005
t_tot= float(input("Hur många sekunder vill du simulera pendeln? "))
t0 = 0

m_1 = 1
m_2 = 1
l_1 = 1
l_2 = 1
g = 9.82

## eulers method
def euler_step():
    t0 = 0
    x1pos_euler = []
    y1pos_euler = []
    x2pos_euler = []
    y2pos_euler = []

    the1 = the1_gbl
    the2 = the2_gbl
    ome1 = ome1_gbl
    ome2 = ome2_gbl

    while t0 < t_tot:
        #########################
        #calculates coordinates of pendulum of mass 1
        x1 = l_1*np.sin(the1)
        x1pos_euler.append(x1)
        y1 = -1 * l_1 * np.cos(the1)
        y1pos_euler.append(y1)

        ## calculates coordinates of mass 2
        x2 = l_1*np.sin(the1) + l_2*np.sin(the2)
        x2pos_euler.append(x2)
        y2 = -1 * l_1*np.cos(the1) - l_2*np.cos(the2)
        y2pos_euler.append(y2)
        #################

        #########################
        #calculate coefficients
        dtheta = the1 - the2
        alpha = (m_1 + m_2)*l_1
        beta = m_2*l_2*np.cos(dtheta)
        gamma = m_2*l_1*np.cos(dtheta)
        delta = m_2*l_2
        epsilon = -1*m_2*l_2 * ome2**2 * np.sin(dtheta) - (m_1 + m_2)*g*np.sin(the1)
        zeta = m_2*l_2 * ome1**2 * np.sin(dtheta) - m_2*g*np.sin(the2)
        ############################

        #########################
        # calculates next step
        domega1 = (delta*epsilon - beta*zeta)/(alpha*delta - beta*gamma) # calculates new domega1
        domega2 = (alpha*zeta - gamma*epsilon)/(alpha*delta - beta*gamma) # calculates new domega2

        ome1 = ome1 + h*domega1 # calculates new omega1
        ome2 = ome2 + h*domega2 # calculates new omega2
        the1 = the1 + h*ome1 # calculates new theta1
        the2 = the2 + h*ome2 # calculates new theta2

        t0 += h # calculates new time
        ########################
    t0 = 0
    return [x1pos_euler, y1pos_euler, x2pos_euler, y2pos_euler]

def rk4_step():
    t0 = 0
    x1pos_rk4 = []
    y1pos_rk4 = []
    x2pos_rk4 = []
    y2pos_rk4 = []

    the1 = the1_gbl
    the2 = the2_gbl
    ome1 = ome1_gbl
    ome2 = ome2_gbl

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
    
    while t0 < t_tot:
        #########################
        #calculates coordinates for mass 1
        x1 = l_1*np.sin(state[0])
        x1pos_rk4.append(x1)
        y1 = -1 * l_1 * np.cos(state[0])
        y1pos_rk4.append(y1)

        ## calculates coordinates for mass 2
        x2 = l_1*np.sin(state[0]) + l_2*np.sin(state[2])
        x2pos_rk4.append(x2)
        y2 = -1 * l_1*np.cos(state[0]) - l_2*np.cos(state[2])
        y2pos_rk4.append(y2)
        ########################

        K1 = derivative(state)
        K2 = derivative(state + h/2 * K1)
        K3 = derivative(state + h/2 * K2)
        K4 = derivative(state + h * K3)

        state = state + h/6 * (K1 + 2*K2 + 2*K3 + K4) ### defines next state vector

        t0 += h # calculates new time
        ########################
    t0=0
    return [x1pos_rk4, y1pos_rk4, x2pos_rk4, y2pos_rk4]

################################
# Get all positions from rk4 method
rk4_results = rk4_step()
x1pos_rk4, y1pos_rk4, x2pos_rk4, y2pos_rk4 = rk4_results

# Get all positions from euler method
euler_results = euler_step()
x1pos_euler, y1pos_euler, x2pos_euler, y2pos_euler = euler_results
###############################
###############################
# Calculate how many steps correspond to the last 2 seconds
steps_last2s = int(2 / h)

# Draw path for last 2 seconds for mass 2 for rk4 method
plt.plot(x2pos_rk4[-steps_last2s:], y2pos_rk4[-steps_last2s:], color='red', ls = ':', lw = '1',label='RK4 färdväg (senaste 2s)')

# Draw path for last 2 seconds for mass 2 for euler method
plt.plot(x2pos_euler[-steps_last2s:], y2pos_euler[-steps_last2s:], color='blue', ls = ':', lw = '1', label='Euler färdväg (senaste 2s)')
################################

########################
#RK4 pendulum

#plots pendulum at final position for rk4 method
plt.plot([0,x1pos_rk4[-1]], [0, y1pos_rk4[-1]], color='red')
plt.plot([x1pos_rk4[-1], x2pos_rk4[-1]], [y1pos_rk4[-1], y2pos_rk4[-1]], color='red')
plt.plot(x1pos_rk4[-1], y1pos_rk4[-1], marker='o', markersize=8, color='red')
plt.plot(x2pos_rk4[-1], y2pos_rk4[-1], marker='o', markersize=8, label='RK4', color='red')
##########################

##########################
## Euler pendulum

#plots pendulum at final position for euler method
plt.plot([0,x1pos_euler[-1]], [0, y1pos_euler[-1]], color='blue') ## draws line for l_1
plt.plot([x1pos_euler[-1], x2pos_euler[-1]], [y1pos_euler[-1], y2pos_euler[-1]], color='blue') ## draws line for l_2
plt.plot(x1pos_euler[-1], y1pos_euler[-1], marker='o', markersize=8, color='blue') ## draws mass 1
plt.plot(x2pos_euler[-1], y2pos_euler[-1], marker='o', markersize=8, label='Euler', color='blue') ## draws mass 2
###########################

############################
### plots legend, title and labels
plt.legend()
plt.title(f'Jämförelse av RK4 och Euler - t={t_tot}s')
plt.xlabel(' x-axel [m]')
plt.ylabel('y-axel [m]')
############################
## sets limits and grid
ax = plt.gca()
ax.set_xlim([-2.5, 2.5])
ax.set_ylim([-2.5, 2.5])
plt.grid()
############################

### saves plot to directory
plt.savefig(f'combined_euler_rk4_plots/combined_euler_rk4_plot_at_t={t_tot}s_the1={the1_gbl}_the2={the2_gbl}.png', dpi=300)

## shows plot
plt.show()

