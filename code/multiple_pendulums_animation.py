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

fig, axis = plt.subplots()
animated1_l_1 = axis.plot([],[], color='blue')[0]
animated1_l_2 = axis.plot([],[], color='blue')[0]
animated1_m1 = axis.plot([],[], 'o', markersize=15, color='red')[0]
animated1_m2 = axis.plot([],[],'o', markersize=15,color='red')[0]
animated1_path_m2 = axis.plot([],[], color='red')[0]

animated2_l_1 = axis.plot([], [], color='black')[0]
animated2_l_2 = axis.plot([],[], color='black')[0]
animated2_m1 = axis.plot([],[], 'o', markersize=15, color='green')[0]
animated2_m2 = axis.plot([],[],'o', markersize=15,color='green')[0]
animated2_path_m2 = axis.plot([],[], color='black')[0]

axis.set_xlim([-2.5,2.5])
axis.set_ylim([-2.5,2.5])
axis.set_title('Animering av två dubbelpendlar - RK4 NY')

plt.grid()

frames=round((t_tot/25)*10**3)
animation_const = len(x2pos)/frames

path_splice_limit=200

def update_data(frame):    
    # first pendulum
    animated1_l_1.set_data([0,x1pos[round(frame*animation_const)]], [0, y1pos[round(frame*animation_const)]])
    animated1_l_2.set_data([x1pos[round(frame*animation_const)], x2pos[round(frame*animation_const)]], [y1pos[round(frame*animation_const)], y2pos[round(frame*animation_const)]])

    animated1_m1.set_data([x1pos[round(frame*animation_const)]],[y1pos[round(frame*animation_const)]])
    animated1_m2.set_data([x2pos[round(frame*animation_const)]],[y2pos[round(frame*animation_const)]])

    animated1_path_m2.set_data(x2pos[:round(frame*animation_const):path_splice_limit], y2pos[:round(frame*animation_const):path_splice_limit])

    
    # second pendulum
    animated2_l_1.set_data([0,X1pos[round(frame*animation_const)]], [0, Y1pos[round(frame*animation_const)]])
    animated2_l_2.set_data([X1pos[round(frame*animation_const)], X2pos[round(frame*animation_const)]], [Y1pos[round(frame*animation_const)], Y2pos[round(frame*animation_const)]])

    animated2_m1.set_data([X1pos[round(frame*animation_const)]],[Y1pos[round(frame*animation_const)]])
    animated2_m2.set_data([X2pos[round(frame*animation_const)]],[Y2pos[round(frame*animation_const)]])

    animated2_path_m2.set_data(X2pos[:round(frame*animation_const):path_splice_limit], Y2pos[:round(frame*animation_const):path_splice_limit])
    
    
    return animated1_l_1, animated1_l_2,  animated1_m1, animated1_m2, animated1_path_m2, animated2_l_1, animated2_l_2, animated2_m1, animated2_m2, animated2_path_m2


animation = FuncAnimation(
    fig=fig,
    func=update_data,
    frames=frames,
    interval=25,
) 

# ax = plt.gca()
# ax.set_xlim([-2.5, 2.5])
# ax.set_ylim([-2.5, 2.5])

# plt.plot(xpos,ypos)
plt.show()
