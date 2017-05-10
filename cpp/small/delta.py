import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

a = 0.1
N = 1000
Nq = N*10

dR = a
dq = 0.0005*(2*np.pi/a)

###################################################################################

iR_vec = np.arange(0,N)
R_vec = dR*iR_vec

# iq_vec = np.arange(0,3*N)
# q_vec = (2.0*np.pi/a/N)*iq_vec

q_vec = np.linspace(-np.pi/a,np.pi/a,Nq)
dq = q_vec[1]-q_vec[0]

J = np.zeros(q_vec.shape,dtype=type(1j))
for R in R_vec:
	J = J + np.exp(1j*q_vec*R)
J = J/N


plt.ion() # enables interactive mode for ploting. No plt.show() is needed!

# ########################################

fig = plt.figure()
ax = fig.gca()
ax.plot([q_vec[0], q_vec[int(Nq/2)-1], q_vec[Nq-1]], [0, 0, 0], color='red', linestyle='none', linewidth=3, marker='o', markersize=10)
ax.plot(q_vec, np.absolute(J), linestyle='solid', linewidth=3, marker='', markersize=10)

# print('sum =', np.absolute(np.sum(J))*dq)
print('sum =', np.sum(np.conj(J)*J)*dq/(2*np.pi/(N*a)))

input("Press Enter to exit...")