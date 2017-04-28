import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

a = 0.01
b = 1*a
D = 1e-1*a
N = 50

def func(R):
	I = 1/np.sqrt((D)**2+(R+b)**2)
	return I


###################################################################################

u_vec = np.arange(-int(N/2),int(N/2))
# u_vec = np.arange(0,N)
R_vec = a*u_vec
ik_vec = np.arange(0,2*N)
k_vec = (2*np.pi/a/N)*(ik_vec-int(N/2))

K, Kp = np.meshgrid(k_vec,k_vec)

V = np.zeros(K.shape,dtype=type(1j))

for R in R_vec:
	for Rp in R_vec:
		I = func(R-Rp)
		V = V + np.exp(1j*R*K)*np.exp(-1j*Rp*Kp)*I

V = V/N

u2_vec = np.arange(0,N)
R2_vec = a*u2_vec
ik_vec = np.arange(0,2*N)
k_vec = (2*np.pi/a/N)*(ik_vec-int(N/2))

K, Kp = np.meshgrid(k_vec,k_vec)

V2 = np.zeros(K.shape,dtype=type(1j))

for R in R2_vec:
	for Rp in R2_vec:
		I = func(R-Rp)
		V2 = V2 + np.exp(1j*R*K)*np.exp(-1j*Rp*Kp)*I

V2 = V2/N

plt.ion() # enables interactive mode for ploting. No plt.show() is needed!

# fig = plt.figure()
# ax = fig.gca(projection='3d')
# surf = ax.plot_surface(K, Kp, np.absolute(V), cmap=cm.coolwarm, linewidth=0, antialiased=True)
# fig.colorbar(surf, shrink=0.5, aspect=5)

fig = plt.figure()
ax = fig.gca()
ax.plot(K.diagonal(0), np.absolute(V.diagonal(0)), color='blue', linestyle='solid', linewidth=6)
ax.plot(K.diagonal(0), np.absolute(V2.diagonal(0)), color='black', linestyle='none', marker='o', linewidth=6, markersize=10)
# ax.plot(K.diagonal(1), np.absolute(V.diagonal(1)), linestyle='dashed', linewidth=3)
plt.ylim((0,2000))

V1d = np.zeros(k_vec.shape, dtype=type(1j))
for R in R_vec:
	I = func(R)
	V1d = V1d + np.exp(+1j*R*k_vec)*I
ax.plot(k_vec, np.absolute(V1d), color='red', linestyle='solid', linewidth=3)

plt.axis('auto')

########################################

fig = plt.figure()
ax = fig.gca()
ax.plot(K.diagonal(0), V.diagonal(0).imag, color='blue', linestyle='solid', linewidth=6)
ax.plot(K.diagonal(0), V2.diagonal(0).imag, color='black', linestyle='none', marker='o', linewidth=6, markersize=10)
# ax.plot(K.diagonal(1), V.diagonal(1).imag, linestyle='dashed', linewidth=3)
plt.ylim((0,2000))

V1d = np.zeros(k_vec.shape, dtype=type(1j))
for R in R_vec:
	I = func(R)
	V1d = V1d + np.exp(+1j*R*k_vec)*I
ax.plot(k_vec, V1d.imag, color='red', linestyle='solid', linewidth=3)

plt.axis('auto')

########################################

fig = plt.figure()
ax = fig.gca()
ax.plot(K.diagonal(0), V.diagonal(0).real, color='blue', linestyle='solid', linewidth=6)
ax.plot(K.diagonal(0), V2.diagonal(0).real, color='black', linestyle='none', marker='o', linewidth=6, markersize=10)
# ax.plot(K.diagonal(1), V.diagonal(1).real, linestyle='dashed', linewidth=3)
plt.ylim((0,2000))

V1d = np.zeros(k_vec.shape, dtype=type(1j))
for R in R_vec:
	I = func(R)
	V1d = V1d + np.exp(+1j*R*k_vec)*I
ax.plot(k_vec, V1d.real, color='red', linestyle='solid', linewidth=3)


plt.axis('auto')


input("Press Enter to exit...")