import matplotlib.pyplot as plt 
import numpy as np
from matplotlib.animation import FuncAnimation

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

while t0 < t_tot + h:
  #calculate coordinates of pendulum
  ## calculates coordinates of mass 1
  x1 = l_1*np.sin(the1)
  x1pos.append(x1)
  y1 = -1 * l_1 * np.cos(the1)
  y1pos.append(y1)

  ## calculates coordinates of mass 2
  x2 = l_1*np.sin(the1) + l_2*np.sin(the2)
  x2pos.append(x2)
  y2 = -1 * l_1*np.cos(the1) - l_2*np.cos(the2)
  y2pos.append(y2)

  #calculate coefficients
  dtheta = the1 - the2
  alpha = (m_1 + m_2)*l_1
  beta = m_2*l_2*np.cos(dtheta)
  gamma = m_2*l_1*np.cos(dtheta)
  delta = m_2*l_2
  epsilon = -1*m_2*l_2 * ome2**2 * np.sin(dtheta) - (m_1 + m_2)*g*np.sin(the1)
  zeta = m_2*l_2 * ome1**2 * np.sin(dtheta) - m_2*g*np.sin(the2)

  domega1 = (delta*epsilon - beta*zeta)/(alpha*delta - beta*gamma) # calculates new domega1
  domega2 = (alpha*zeta - gamma*epsilon)/(alpha*delta - beta*gamma) # calculates new domega2

  ome1 = ome1 + h*domega1 # calculates new omega1
  ome2 = ome2 + h*domega2 # calculates new omega2

  the1 = the1 + h*ome1 # calculates new theta1
  the1_list.append(the1)
  the2 = the2 + h*ome2 # calculates new theta2
  the2_list.append(the2)

  t0 += h

plt.grid()

frames=round((t_tot/25)*10**3)
animation_const = len(x2pos)/frames

path_limit = 10

ax = plt.gca()
ax.set_xlim([-2.5, 2.5])
ax.set_ylim([-2.5, 2.5])

plt.plot(x2pos,y2pos)
plt.show()
