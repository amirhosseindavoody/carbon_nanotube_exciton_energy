import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

plt.ion() # enables interactive mode for ploting. No plt.show() is needed!

directory = "/home/amirhossein/research/test/"

data = np.loadtxt(directory+"cnt1.v_dir.real.dat", skiprows=0)
kr_vec = data[1:,0]
krp_vec = data[0,1:]
V_dir_real = data[1:,1:]
data = np.loadtxt(directory+"cnt1.v_dir.imag.dat", skiprows=0)
V_dir_imag = data[1:,1:]

V_dir_mag = np.sqrt(np.power(V_dir_real,2)+np.power(V_dir_imag,2))

fig=plt.figure()
ax = fig.gca()
im = ax.imshow(V_dir_real, interpolation='bilinear', cmap=cm.coolwarm, origin='lower', extent=[kr_vec[0], kr_vec[-1], krp_vec[0], krp_vec[-1]], vmax=abs(V_dir_real).max(), vmin=-abs(V_dir_real).max())
fig.colorbar(im)

fig=plt.figure()
ax = fig.gca()
im = ax.imshow(V_dir_imag, interpolation='bilinear', cmap=cm.coolwarm, origin='lower', extent=[kr_vec[0], kr_vec[-1], krp_vec[0], krp_vec[-1]], vmax=abs(V_dir_imag).max(), vmin=-abs(V_dir_imag).max())
fig.colorbar(im)

fig=plt.figure()
ax = fig.gca()
im = ax.imshow(V_dir_mag, interpolation='bilinear', cmap=cm.coolwarm, origin='lower', extent=[kr_vec[0], kr_vec[-1], krp_vec[0], krp_vec[-1]], vmax=V_dir_mag.max(), vmin=V_dir_mag.min())
fig.colorbar(im)

fig=plt.figure()
ax = fig.gca()
ax.plot(V_dir_mag.diagonal())


input("Press Enter to exit...")