import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

plt.ion() # enables interactive mode for ploting. No plt.show() is needed!

directory = "/home/amirhossein/research/test/"

data = np.loadtxt(directory+"cnt1.ex_energy.dat", skiprows=0)

ex_energy = np.zeros(data.shape)
ex_energy[:,0] = data[:,0]
ex_energy[:,1] = data[:,0]

eh_energy = np.zeros(data.shape)
eh_energy[:,0] = np.sort(data[:,1])
eh_energy[:,1] = np.sort(data[:,1])


ex_vec = [0, 0.4]
eh_vec = [0.6, 1]

fig = plt.figure()
ax = fig.gca()
ax.plot(eh_vec,np.transpose(eh_energy), linestyle="solid", color="blue", linewidth=2, marker="")
ax.plot(ex_vec,np.transpose(ex_energy), linestyle="solid", color="red", linewidth=2, marker="")

fig = plt.figure()
ax = fig.gca()
ax.plot(eh_energy[:,0], linestyle="solid", color="blue", linewidth=2, marker="o")
ax.plot(ex_energy[:,0], linestyle="solid", color="red", linewidth=2, marker="o")


data = np.loadtxt(directory+"cnt1.ex_psi.dat", skiprows=0)
ex_psi = data

fig = plt.figure()
ax = fig.gca()
for j in range(0,10):
	ax.plot(ex_psi[:,j], linestyle="solid", linewidth=2, marker="")
	print(j, ex_energy[j,0])
	input()

input("Press Enter to exit...")