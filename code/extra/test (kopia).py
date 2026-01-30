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

X_state = np.array[ome1, ome2, domega1, domega2]

while t0 < t_tot + h:
    def epsilon(ome2):
        return -1*m_2*l_2 * ome2**2 * np.sin(dtheta) - (m_1 + m_2)*g*np.sin(the1)
    def zeta(ome1):
        return m_2*l_2 * ome1**2 * np.sin(dtheta) - m_2*g*np.sin(the2)
    def domega1(epsilon, zeta):
        return (delta*epsilon - beta*zeta)/(alpha*delta - beta*gamma)
    def domega2(epsilon, zeta):
        return (alpha*zeta - gamma*epsilon)/(alpha*delta - beta*gamma)
    
    dtheta = the1 - the2
    alpha = (m_1 + m_2)*l_1
    beta = m_2*l_2*np.cos(dtheta)
    gamma = m_2*l_1*np.cos(dtheta)
    delta = m_2*l_2

    

    K1 = domega1(epsilon(ome2), zeta(ome1))

    ome1 = ome1 + h/6 * (K1 + 2*K2 + 2*K3 + K4)
