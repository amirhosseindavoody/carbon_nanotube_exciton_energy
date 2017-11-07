# import matplotlib as mpl
# import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
#
# plt.ion() # enables interactive mode for ploting. No plt.show() is needed!
#
#
# directory = "/home/amirhossein/research/test/"
#
# data = np.loadtxt(directory+"cnt1.electron_energy.dat", skiprows=0)
# wave_vec = data[:,0]
# energy = np.transpose(data[:,1:])
#
#
# fig = plt.figure()
# ax = fig.gca()
# for i in range(1,energy.shape[0]):
# 	ax.plot(wave_vec,energy[i,:], linestyle="solid", linewidth=2, marker="")
#
# Nu = int(energy.shape[0]/2)
# bands = Nu+np.arange(-3,3)
# fig = plt.figure()
# ax = fig.gca()
# for i in bands:
# 	ax.plot(wave_vec,energy[i,:], linestyle="solid", linewidth=2, marker="")
#
# input("Press Enter to exit...")

import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt

plt.ion() # enables interactive mode for ploting. No plt.show() is needed!


directory = "/Users/amirhossein/research/exciton_energy/"

data = np.loadtxt(directory+"cnt1.el_energy_redu.dat", skiprows=2)
# data = np.loadtxt(directory+"cnt1.el_energy_full.dat", skiprows=2)
data = np.transpose(data)


fig = plt.figure()
ax = fig.gca()
ax.plot(data, linestyle="solid", linewidth=2, marker="")

input("Press Enter to exit...")
