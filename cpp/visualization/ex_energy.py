import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

plt.ion() # enables interactive mode for ploting. No plt.show() is needed!

directory = "/home/amirhossein/research/test/"

data = np.loadtxt(directory+"cnt1.ex_energy.dat", skiprows=0)
ex_energy = np.zeros(data.shape)
elec_hole = np.zeros(data.shape)
ex_energy[:,0] = data[:,0]
ex_energy[:,1] = data[:,0]
elec_hole[:,0] = data[:,1]
elec_hole[:,1] = data[:,1]
# elec_hole = data
# ex_energy[:,1] = ex_energy[:,0]
# elec_hole[:,0] = elec_hole[:,1]

xvec = [0, 1]

fig = plt.figure()
ax = fig.gca()
for i in range(0,data.shape[0]):
	# print(data[i,0], data[i,1])
	# input()
	ax.plot(xvec[:],ex_energy[i,:], linestyle="solid", color="red", linewidth=2, marker="")
	# print(ex_energy[i,0], elec_hole[i,0])
	# input("")
	ax.plot(xvec[:],elec_hole[i,:], linestyle="solid", color="blue", linewidth=2, marker="")

# # fig=plt.figure()
# # ax = fig.gca()
# # im = ax.imshow(V_dir_imag, interpolation='bilinear', cmap=cm.coolwarm, origin='lower', extent=[kr_vec[0], kr_vec[-1], krp_vec[0], krp_vec[-1]], vmax=abs(V_dir_imag).max(), vmin=-abs(V_dir_imag).max())
# # fig.colorbar(im)

# fig=plt.figure()
# ax = fig.gca()
# im = ax.imshow(V_dir_mag, interpolation='bilinear', cmap=cm.coolwarm, origin='lower', extent=[kr_vec[0], kr_vec[-1], krp_vec[0], krp_vec[-1]], vmax=V_dir_mag.max(), vmin=V_dir_mag.min())
# fig.colorbar(im)

# # fig=plt.figure()
# # ax = fig.gca()
# # ax.plot(V_dir_mag.diagonal())


# data = np.loadtxt(directory+"cnt1.v_xch.real.dat", skiprows=0)
# kr_vec = data[1:,0]
# krp_vec = data[0,1:]
# V_xch_real = data[1:,1:]
# data = np.loadtxt(directory+"cnt1.v_xch.imag.dat", skiprows=0)
# V_xch_imag = data[1:,1:]

# V_xch_mag = np.sqrt(np.power(V_xch_real,2)+np.power(V_xch_imag,2))

# # fig=plt.figure()
# # ax = fig.gca()
# # im = ax.imshow(V_xch_real, interpolation='bilinear', cmap=cm.coolwarm, origin='lower', extent=[kr_vec[0], kr_vec[-1], krp_vec[0], krp_vec[-1]], vmax=abs(V_xch_real).max(), vmin=-abs(V_xch_real).max())
# # fig.colorbar(im)

# # fig=plt.figure()
# # ax = fig.gca()
# # im = ax.imshow(V_xch_imag, interpolation='bilinear', cmap=cm.coolwarm, origin='lower', extent=[kr_vec[0], kr_vec[-1], krp_vec[0], krp_vec[-1]], vmax=abs(V_xch_imag).max(), vmin=-abs(V_xch_imag).max())
# # fig.colorbar(im)

# fig=plt.figure()
# ax = fig.gca()
# im = ax.imshow(V_xch_mag, interpolation='bilinear', cmap=cm.coolwarm, origin='lower', extent=[kr_vec[0], kr_vec[-1], krp_vec[0], krp_vec[-1]], vmax=V_xch_mag.max(), vmin=V_xch_mag.min())
# fig.colorbar(im)

input("Press Enter to exit...")