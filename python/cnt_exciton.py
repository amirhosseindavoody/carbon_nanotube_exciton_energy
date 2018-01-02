
# coding: utf-8

# In[8]:


import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os

# get_ipython().run_line_magic('matplotlib', 'inline')


def gcd(x, y):
    gcd = 1
    
    if x % y == 0:
        return y
    
    for k in range(int(y / 2), 0, -1):
        if x % k == 0 and y % k == 0:
            gcd = k
            break  
    return gcd


# # Physical constants

# In[9]:


eV = 1.6e-19


# # Carbon structure constants

# In[10]:


n, m = 13, 2
print("chirality: ({},{})".format(n,m))

a_cc = 1.42e-10
print("carbon-carbon distance: {}".format(a_cc))

a_l = np.sqrt(3)*a_cc
print("lattice constant: {}".format(a_l))

# tight-binding parameters
e2p = 0
t0 = 2.7*eV
s0 = 0

# Coulomb interaction parameters
Upp = 11.3*eV
kappa = 2.0

# length of the CNT
number_of_unit_cells = 100


# # Graphene unit cells, chirality vector, translational vector
# 

# In[11]:


a1 = np.array([np.sqrt(3)/2, +1/2])*a_l
a2 = np.array([np.sqrt(3)/2, -1/2])*a_l
print("a1: {}".format(a1))
print("a2: {}".format(a2))

b1 = np.array([1/np.sqrt(3)*2*np.pi, 2*np.pi])/a_l
b2 = np.array([1/np.sqrt(3)*2*np.pi, -2*np.pi])/a_l
print("b1: {}".format(b1))
print("b2: {}".format(b2))

aCC_vec = 1./3.*(a1+a2)
print("aCC_vec: {}".format(aCC_vec))

ch_vec = n*a1+m*a2
print("ch_vec: {}".format(ch_vec))

ch_len = np.linalg.norm(ch_vec)
print("ch_len: {}".format(ch_len))

radius = ch_len/2/np.pi
print("radius: {}".format(radius))

dR = gcd(2*n+m, n+2*m)
t1 = int((2*m+n)/dR)
t2 = int(-(2*n+m)/dR)
t_vec = t1*a1+ t2*a2
print("t_vec: {}".format(t_vec))

Nu = int(2*(n**2+m**2+n*m)/dR)
print("Total number of hexagons in CNT unit cell: {}".format(Nu))
print()

cos_theta = ch_vec[0]/ch_len
sin_theta = ch_vec[1]/ch_len
rot_mat = np.array([[cos_theta, sin_theta], [-sin_theta, cos_theta]])

ch_vec = np.matmul(rot_mat, ch_vec)
t_vec = np.matmul(rot_mat, t_vec)
a1 = np.matmul(rot_mat, a1)
a2 = np.matmul(rot_mat, a2)
b1 = np.matmul(rot_mat, b1)
b2 = np.matmul(rot_mat, b2)
aCC_vec = np.matmul(rot_mat, aCC_vec)
print("ch_vec: {}".format(ch_vec))
print("t_vec: {}".format(t_vec))
print("a1: {}".format(a1))
print("a2: {}".format(a2))
print("b1: {}".format(b1))
print("b2: {}".format(b2))
print("aCC_vec: {}".format(aCC_vec))

# # Reciprocal lattice of CNT

# In[12]:


K1 = (-t2*b1+t1*b2)/Nu
K2 = (+m*b1 -n*b2)/Nu
K2_normed = K2/np.linalg.norm(K2)
nk = number_of_unit_cells
dk_l = K2/nk

print("K1: {}".format(K1))
print("K2: {}".format(K2))
print("K2_normed: {}".format(K2_normed))
print("nk: {}".format(nk))
print("dk_l: {}".format(dk_l))


# # Atom coordinates

# In[14]:


def get_atom_coordinates():
    pos_a = np.zeros((Nu,2))
    pos_b = np.zeros((Nu,2))
    
    count = 0
    for i in range(0, t1+n+1):
        for j in range(t2, m+1):
            
            flag1 = (float(t2*i)/float(t1) <= j)
            flag2 = (float(m*i)/float(n) >= j)
            flag3 = (float(t2*(i-n))/float(t1) > (j-m))
            flag4 = (float(m*(i-t1))/float(n) < (j-t2))
            
            if (flag1 and flag2 and flag3 and flag4):
                pos_a[count,:] = i*a1 + j*a2
                pos_b[count,:] = pos_a[count,:] + aCC_vec
                
                if(pos_a[count,0] > ch_vec[0]):
                    pos_a[count,0] -= ch_vec[0]
                if(pos_a[count,0] < 0.0):
                    pos_a[count,0] += ch_vec[0]
                if (pos_a[count,1] > ch_vec[1]):
                    pos_a[count,1] -= ch_vec[1]
                if(pos_a[count,1] < 0.0):
                    pos_a[count,1] += ch_vec[1]

                if(pos_b[count,0] > ch_vec[0]):
                    pos_b[count,0] -= ch_vec[0]
                if(pos_b[count,0] < 0.0):
                    pos_b[count,0] += ch_vec[0]
                if(pos_b[count,0] > ch_vec[1]):
                    pos_b[count,1] -= ch_vec[1]
                if(pos_b[count,1] < 0.0):
                    pos_b[count,1] += ch_vec[1]
                
                count += 1
    
    assert (count==Nu), "count: {}, Nu: {}".format(count, Nu)
    return pos_a, pos_b

pos_a, pos_b = get_atom_coordinates()

# fig=plt.figure(figsize=(18, 9), dpi= 80, facecolor='w', edgecolor='k')
plt.plot(pos_a[:,0], pos_a[:,1], linestyle='none', marker='o')
plt.plot(pos_b[:,0], pos_b[:,1], linestyle='none', marker='o')
plt.axes().set_aspect('equal')
plt.show()
