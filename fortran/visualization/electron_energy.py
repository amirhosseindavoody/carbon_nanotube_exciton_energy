import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

plt.ion() # enables interactive mode for ploting. No plt.show() is needed!


directory = "/home/amirhossein/research/test/"

data = np.loadtxt(directory+"cnt1.c_band.dat", skiprows=0)
wave_vec = data[0,:]
c_energy = np.transpose(data[1:,:])

data = np.loadtxt(directory+"cnt1.v_band.dat", skiprows=0)
v_energy = np.transpose(data[1:,:])


fig = plt.figure()
ax = fig.gca()
ax.plot(wave_vec,c_energy, linestyle="solid", linewidth=2, marker="")
ax.plot(wave_vec,v_energy, linestyle="solid", linewidth=2, marker="")
ax.plot(wave_vec,c_energy[:,3], linestyle="solid", color="black", linewidth=4, marker="")
ax.autoscale(enable=True, axis='x', tight=True)

# fig = plt.figure()
# ax = fig.gca()
# for i in range(1,energy.shape[0]):
# 	ax.plot(wave_vec,energy[i,:], linestyle="solid", linewidth=2, marker="")


# fig.tight_layout()
# ax.set_ylim([-1, 1])

# fig = plt.figure()
# ax = fig.gca()
# energy = np.loadtxt(directory+"cnt1.pi.dat", skiprows=0)
# for i in range(1,energy.shape[1]):
# 	ax.plot(energy[:,0],energy[:,i], linestyle="solid", linewidth=2, marker="")

# fig = plt.figure()
# ax = fig.gca()
# energy = np.loadtxt(directory+"cnt1.vq.dat", skiprows=0)
# for i in range(1,energy.shape[1]):
# 	ax.plot(energy[:,0],energy[:,i], linestyle="solid", linewidth=2, marker="")

# plt.yscale('log')


# fig = plt.figure()
# ax = fig.gca()

# directory = "/home/amirhossein/research/test/"

# psi = np.loadtxt(directory+"cnt1.electron_psi.dat", skiprows=0)
# for i in range(1,psi.shape[1]):
# 	ax.cla();
# 	ax.plot(psi[:,0],psi[:,i], linestyle="solid", linewidth=2, marker="")
# 	# input()

input("Press Enter to exit...")