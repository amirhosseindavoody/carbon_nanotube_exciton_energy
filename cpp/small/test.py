import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

a = 0.01
N = 100

u_vec = np.arange(0,N)
R_vec = a*u_vec
k_vec = (2*np.pi/a/N)*np.arange(0,N)

K, Kp = np.meshgrid(k_vec,k_vec)
# R, Rp = np.meshgrid(R_vec,R_vec)

V = np.zeros(K.shape,dtype=type(1j))

for R in R_vec:
	for Rp in R_vec:
		I = 1/np.sqrt((a/10)**2+(R-Rp)**2)
		V = V + np.exp(1j*R*K)*np.exp(-1j*Rp*Kp)*I

V = V/N

plt.ion() # enables interactive mode for ploting. No plt.show() is needed!
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(K, Kp, np.absolute(V), cmap=cm.coolwarm, linewidth=0, antialiased=True)
fig.colorbar(surf, shrink=0.5, aspect=5)

input("Press Enter to exit...")