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
    K1 = derivative(state)
    K2 = derivative(state + h/2 * K1)
    K3 = derivative(state + h/2 * K2)
    K4 = derivative(state + h * K3)

    state = state + h/6 * (K1 + 2*K2 + 2*K3 + K4)

    the1_list.append(state[0])
    the2_list.append(state[2])

    t0 += h

fig, axis = plt.subplots()
# animated_l_1 = axis.plot([],[], color='blue')[0]
# animated_l_2 = axis.plot([],[], color='blue')[0]
# animated_m1 = axis.plot([],[], 'o', markersize=15, color='red')[0]
# animated_m2 = axis.plot([],[],'o', markersize=15,color='red')[0]
# animated_path_m2 = axis.plot([],[], color='red')[0]
# animated_the1 = axis.plot([], [],'o', markersize=5, color='blue')[0]
# animated_the1_path = axis.plot([],[], color='blue')[0]
# animated_the2 = axis.plot([], [], 'o', markersize=5, color='red')[0]
# animated_the2_path = axis.plot([],[], color='red')[0]

the1_max = abs(np.array(the1_list)).max() + 0.01 * abs(np.array(the1_list)).max()
the2_max = abs(np.array(the2_list)).max() + 0.01 * abs(np.array(the2_list)).max()
if the1_max > the2_max:
    axis.set_ylim([-the1_max, the1_max])
else:
    axis.set_ylim([-the2_max, the2_max])

axis.set_xlim([0,t_tot])

axis.set_title(f'angles_as_functions_of_time_plot, the1={the1}_the2={the2}_t={t_tot}s')

plt.grid()

frames=round((t_tot/25)*10**3)
animation_const = len(the1_list)/frames

path_splice_limit=200
t_list = np.linspace(0, t_tot, len(the1_list))

# def update_data(frame):    
#     # animated_l_1.set_data([0,x1pos[round(frame*animation_const)]], [0, y1pos[round(frame*animation_const)]])
#     # animated_l_2.set_data([x1pos[round(frame*animation_const)], x2pos[round(frame*animation_const)]], [y1pos[round(frame*animation_const)], y2pos[round(frame*animation_const)]])

#     # animated_m1.set_data([x1pos[round(frame*animation_const)]],[y1pos[round(frame*animation_const)]])
#     # animated_m2.set_data([x2pos[round(frame*animation_const)]],[y2pos[round(frame*animation_const)]])

#     # animated_path_m2.set_data(x2pos[:round(frame*animation_const):path_splice_limit], y2pos[:round(frame*animation_const):path_splice_limit])
#     # return animated_l_1, animated_l_2,  animated_m1, animated_m2, animated_path_m2
#     animated_the1.set_data([t_list[round(frame*animation_const)]], [the1_list[round(frame*animation_const)]])
#     animated_the1_path.set_data([t_list[:round(frame*animation_const):path_splice_limit]], [the1_list[:round(frame*animation_const):path_splice_limit]])
#     animated_the2.set_data([t_list[round(frame*animation_const)]], [the2_list[round(frame*animation_const)]])
#     animated_the2_path.set_data([t_list[:round(frame*animation_const):path_splice_limit]], [the2_list[:round(frame*animation_const):path_splice_limit]])
#     return animated_the1, animated_the2, animated_the1_path, animated_the2_path


# animation = FuncAnimation(
#     fig=fig,
#     func=update_data,
#     frames=frames,
#     interval=25,
# ) 

# ax = plt.gca()
#axis.set_xlim([0, 2.5])
# axis.set_ylim([-2.5, 2.5])

# plt.plot(xpos,ypos)
plt.plot(t_list, the1_list, label='Theta 1')
plt.plot(t_list, the2_list, label='Theta 2')

plt.savefig(f'angles_as_functions_of_time_plot, the1={the1}_the2={the2}_t={t_tot}s.png', dpi=300)

plt.show()
