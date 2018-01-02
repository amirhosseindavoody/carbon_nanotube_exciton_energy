import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os

print("\n*****\n")

# physical constants
eV = 1.6e-19

# carbon structure constants
n, m = 13, 0
print("chirality: ({},{})".format(n,m))

a_cc = 1.42e-10
print("carbon-carbon distance: {}".format(a_cc))

a_l = np.sqrt(3)*a_cc
print("lattice constant: {}".format(a_l))

e2p = 0
t0 = 2.7*eV
s0 = 0

Upp = 11.3*eV
kappa = 2.0

number_of_unit_cells = 100

