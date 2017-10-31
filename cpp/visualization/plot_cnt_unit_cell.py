# import matplotlib as mpl
# import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
#
# plt.ion() # enables interactive mode for ploting. No plt.show() is needed!
#
# fig = plt.figure()
# ax = fig.gca()
#
# directory = "/Users/amirhossein/Desktop/"
#
# pos_2d = np.loadtxt(directory+"cnt1.pos_2d.dat", skiprows=1)
# ax.plot(pos_2d[:,0],pos_2d[:,1], linestyle="none", linewidth=2, color="blue", marker="o")
# ax.set_aspect('equal', 'datalim')
#
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# pos_3d = np.loadtxt(directory+"cnt1.pos_3d.dat", skiprows=1)
# ax.plot(pos_3d[:,0],pos_3d[:,1],pos_3d[:,2], linestyle="none", linewidth=2, color="red", marker="o")
# ax.set_aspect('equal', 'datalim')
#
# input("Press Enter to exit...")


import numpy as np
import matplotlib.pyplot as plt

directory = "/Users/amirhossein/Desktop/"

pos_2d = np.loadtxt(directory+"coordinates.dat", skiprows=0)
plt.plot(pos_2d[:,0],pos_2d[:,1], linestyle="none", linewidth=2, color="blue", marker="o")
plt.axes().set_aspect('equal', 'datalim')

plt.show()

# fig = plt.figure()
# ax = fig.gca(projection='3d')
# pos_3d = np.loadtxt(directory+"cnt1.pos_3d.dat", skiprows=1)
# ax.plot(pos_3d[:,0],pos_3d[:,1],pos_3d[:,2], linestyle="none", linewidth=2, color="red", marker="o")
# ax.set_aspect('equal', 'datalim')

# input("Press Enter to exit...")
