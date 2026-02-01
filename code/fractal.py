import matplotlib.pyplot as plt 
import numpy as np
from matplotlib.animation import FuncAnimation

## code based on euler.py
## but reconverted to runge-kutta method as to increase accuracy

h = float(input("Ange step size för RK4: ")) ## t.ex h=0.001
t_tot = float(input("Ange tiden för varje simulering: "))
ome1 = 0
ome2 = 0

the1_list = []
the2_list = []

stable_list = []
unstable_list = []

the1_splice = np.linspace(-1*np.pi, np.pi,10)
the2_splice = np.linspace(-1*np.pi, np.pi, 10)

for i in range(len(the1_splice)):
    for k in range(len(the2_splice)):  
        the1_list.append(the1_splice[i])
        the2_list.append(the2_splice[k])

print("Done splice of theta values")

unstable_the1_list = []
unstable_the2_list = []
stable_the1_list = []
stable_the2_list = []

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

m_1 = 1
m_2 = 1
l_1 = 1
l_2 = 1
g = 9.82

for i in range(len(the1_list)):
    the1 = the1_list[i]
    the2 = the2_list[i]

    t0 = 0

    state = np.array([the1, ome1, the2, ome2])

    stable_marker = True
    while t0 < t_tot + h:
        K1_next = derivative(state)
        K2_next = derivative(state + h/2 * K1_next)
        K3_next = derivative(state + h/2 * K2_next)
        K4_next = derivative(state + h * K3_next)

        state_next = state + h/6 * (K1_next + 2*K2_next + 2*K3_next + K4_next)
        
        
        #calculate coordinates of pendulum
        ## calculates coordinates of mass 1
        x1 = l_1*np.sin(state[0])
        #x1pos.append(x1)
        y1 = -1 * l_1 * np.cos(state[0])
        #y1pos.append(y1)

        x1_next = l_1*np.sin(state_next[0])
        #next_x1pos.append(x1_next)
        y1_next = -1 * l_1 * np.cos(state_next[0])
        #next_y1pos.append(y1_next)

        state = state_next

        t0 += h

        if (x1 < 0 and y1 > 0) and (x1_next > 0 and y1_next > 0):
            # unstable_the1_list.append([the1])
            # unstable_the2_list.append([the2])
            stable_marker = False
            break
        if (x1 > 0 and y1 > 0) and (x1_next < 0 and y1_next > 0):
            # unstable_list.append([the1,the2])
            stable_marker = False
            break

    if stable_marker:
        stable_the1_list.append(the1)
        stable_the2_list.append(the2)



# print(f"stable_list: {stable_list}")
# print(f"unstable_list: {unstable_list}")

ax = plt.gca()
ax.set_xlim([-np.pi, np.pi])
ax.set_ylim([-np.pi, np.pi])
plt.plot(stable_the1_list,stable_the2_list, ".")
plt.savefig('figure.png')
plt.show()



# fig, axis = plt.subplots()
# animated_l_1 = axis.plot([],[], color='blue')[0]
# animated_l_2 = axis.plot([],[], color='blue')[0]
# animated_m1 = axis.plot([],[], 'o', markersize=15, color='red')[0]
# animated_m2 = axis.plot([],[],'o', markersize=15,color='red')[0]
# animated_path_m2 = axis.plot([],[], color='red')[0]

# axis.set_xlim([-2.5,2.5])
# axis.set_ylim([-2.5,2.5])
# axis.set_title('Animering av dubbelpendel - RK4 NY')

# plt.grid()

# frames=round((t_tot/25)*10**3)
# animation_const = len(x2pos)/frames

# path_splice_limit=10

# def update_data(frame):    
#     animated_l_1.set_data([0,x1pos[round(frame*animation_const)]], [0, y1pos[round(frame*animation_const)]])
#     animated_l_2.set_data([x1pos[round(frame*animation_const)], x2pos[round(frame*animation_const)]], [y1pos[round(frame*animation_const)], y2pos[round(frame*animation_const)]])

#     animated_m1.set_data([x1pos[round(frame*animation_const)]],[y1pos[round(frame*animation_const)]])
#     animated_m2.set_data([x2pos[round(frame*animation_const)]],[y2pos[round(frame*animation_const)]])

#     animated_path_m2.set_data(x2pos[:round(frame*animation_const):path_splice_limit], y2pos[:round(frame*animation_const):path_splice_limit])
#     return animated_l_1, animated_l_2,  animated_m1, animated_m2, animated_path_m2

# animation = FuncAnimation(
#     fig=fig,
#     func=update_data,
#     frames=frames,
#     interval=25,
# ) 

# # ax = plt.gca()
# # ax.set_xlim([-2.5, 2.5])
# # ax.set_ylim([-2.5, 2.5])

# # plt.plot(xpos,ypos)
# plt.show()