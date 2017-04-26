import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

plt.ion() # enables interactive mode for ploting. No plt.show() is needed!

fig = plt.figure()
ax = fig.gca()

directory = "/home/amirhossein/research/test/"

energy = np.loadtxt(directory+"cnt1.electron_energy.dat", skiprows=0)
for i in range(1,energy.shape[1]):
	ax.plot(energy[:,0],energy[:,i], linestyle="solid", linewidth=2, marker="")
# fig.tight_layout()
# ax.set_ylim([-1, 1])

fig = plt.figure()
ax = fig.gca()
energy = np.loadtxt(directory+"cnt1.pi.dat", skiprows=0)
for i in range(1,energy.shape[1]):
	ax.plot(energy[:,0],energy[:,i], linestyle="solid", linewidth=2, marker="")

fig = plt.figure()
ax = fig.gca()
energy = np.loadtxt(directory+"cnt1.vq.dat", skiprows=0)
for i in range(1,energy.shape[1]):
	ax.plot(energy[:,0],energy[:,i], linestyle="solid", linewidth=2, marker="")

# plt.yscale('log')

input("Press Enter to exit...")