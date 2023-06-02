import numpy as np
from qutip import * 
import qutip as q
import matplotlib.pyplot as plt
from matplotlib import cm
from States import plot_wigner3d
from States import wigner_cmap 

N = 20 #why not : Hilbert Space dimensions
# xvec = np.linspace(-6, 6, 1500)
# fock = q.fock(N, 4)
# W_fock = q.wigner(fock, xvec, xvec)
# q.plot_wigner_fock_distribution(fock, figsize=(10, 4), colorbar=True)

# States.plot_wigner3d(fock,xvec)
# plt.show()


# psi = (basis(10, 8) + basis(10,3)+ basis(10,5)).unit()
xvec = np.linspace(-5, 5, 1500)
# States.plot_wigner3d(psi,xvec) 
# plt.show()


# # Create a new figure with vertical subplots
# fig_vertical = plt.figure(figsize=(8, 12))

# # Call plot_wigner3d for the first subplot
# fig_vertical = plot_wigner3d(psi, xvec, fig_vertical, 311)

# # Call plot_wigner3d for the second subplot
# fig_vertical = plot_wigner3d(psi, xvec, fig_vertical, 312)

# # Call plot_wigner3d for the third subplot
# fig_vertical = plot_wigner3d(psi, xvec, fig_vertical, 313)

# # Show the figure with subplots
# plt.show()
# #############

psi_list = [(basis(10, 8) + basis(10,3)+ basis(10,5)).unit(), q.fock(N, 3), q.coherent(N, np.sqrt(10)) + q.coherent(N, -np.sqrt(10))]

fig_vertical = plt.figure(figsize=(8, 12))

num_subplots = len(psi_list)
W_list = [wigner(psi, xvec, xvec) for psi in psi_list]
wmap_list = [wigner_cmap(W) for W in W_list]
subplot_positions = [311, 312, 313]

ax = fig_vertical.add_subplot(subplot_positions[0])
plot_wigner3d(psi_list[0],xvec)

ax = fig_vertical.add_subplot(subplot_positions[1])
plot_wigner3d(psi_list[1], xvec)

ax = fig_vertical.add_subplot(subplot_positions[2])
plot_wigner3d(psi_list[2], xvec)

# Call plot_wigner3d to add subplots to the figure
# fig_vertical = plot_wigner3d(psi_list, xvec, fig_vertical, subplot_positions)

# Show the figure with subplots
plt.show()
