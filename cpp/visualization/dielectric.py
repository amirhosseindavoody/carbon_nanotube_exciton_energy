import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

plt.ion() # enables interactive mode for ploting. No plt.show() is needed!

directory = "/home/amirhossein/research/test/"

data = np.loadtxt(directory+"cnt1.vq.real.dat", skiprows=0)
wave_vec = data[0,:]
vq_real = data[1:,:]
data = np.loadtxt(directory+"cnt1.vq.imag.dat", skiprows=1)
vq_imag = data[:,:]

Nb = int(np.sqrt(vq_real.shape[0]))
Nq = wave_vec.shape[0];

fig = plt.figure()
ax = fig.gca()
for ib in range(0,Nb):
	for ibp in range(ib, Nb):
		ax.plot(wave_vec[:],vq_imag[ib*Nb+ibp,:], linestyle="solid", linewidth=2, marker="")

fig = plt.figure()
ax = fig.gca()
for ib in range(0,Nb):
	for ibp in range(ib, Nb):
		ax.plot(wave_vec[:],vq_real[ib*Nb+ibp,:], linestyle="solid", linewidth=2, marker="")

# diff_real = np.zeros((Nb*Nb,int(Nq/2)-1))
# diff_imag = np.zeros((Nb*Nb,int(Nq/2)-1))
# for i in range(0,int(Nq/2)-1):
# 	diff_real[:,i] = vq_real[:,i+1] - vq_real[:,-1-i]
# 	diff_imag[:,i] = vq_imag[:,i+1] + vq_imag[:,-1-i]

# fig = plt.figure()
# ax = fig.gca()
# for ib in range(0,Nb):
# 	for ibp in range(ib, Nb):
# 		ax.plot(wave_vec[0:int(Nq/2)-1],diff_imag[ib*Nb+ibp,:], linestyle="solid", linewidth=2, marker="")

# fig = plt.figure()
# ax = fig.gca()
# for ib in range(0,Nb):
# 	for ibp in range(ib, Nb):
# 		ax.plot(wave_vec[0:int(Nq/2)-1],diff_real[ib*Nb+ibp,:], linestyle="solid", linewidth=2, marker="")

input("Press Enter to exit...")