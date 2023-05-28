import numpy as np
from qutip import * 
import qutip as q
import matplotlib.pyplot as plt
from matplotlib import cm


N = 30
xvec = np.linspace(-5, 5, 500)

psi = (basis(10, 0) + basis(10, 3) + basis(10, 9)).unit()

W = wigner(psi, xvec, xvec)

wmap = wigner_cmap(W)  # Generate Wigner colormap

# # nrm = plt.colors.Normalize(-W.max(), W.max())

# fig, axes = plt.subplots(1, 2, figsize=(10, 4))

# plt1 = axes[0].contourf(xvec, xvec, W, 100, cmap=cm.RdBu)

# axes[0].set_title("Standard Colormap");

# cb1 = fig.colorbar(plt1, ax=axes[0])

# plt2 = axes[1].contourf(xvec, xvec, W, 100, cmap=wmap)  # Apply Wigner colormap

# axes[1].set_title("Wigner Colormap");
# cb2 = fig.colorbar(plt2, ax=axes[1])
# fig.tight_layout()
# plt.show()


def plot_wigner_2d_3d(psi):
    #fig, axes = plt.subplots(1, 2, subplot_kw={'projection': '3d'}, figsize=(12, 6))
    fig = plt.figure(figsize=(20, 8))

    ax = fig.add_subplot(1, 2, 1)        
    plot_wigner(psi, fig=fig, ax=ax, alpha_max=6, cmap=wmap, colorbar = True);

    ax = fig.add_subplot(1, 2, 2, projection='3d')
    plot_wigner(psi, fig=fig, ax=ax, projection='3d', alpha_max=6, cmap=wmap);
    
    plt.close(fig)
    return fig

fock = q.fock(N, 3)
W_fock = q.wigner(fock, xvec, xvec)
q.plot_wigner_fock_distribution(fock, figsize=(10, 4), colorbar=True)

plot_wigner_2d_3d(fock)
plt.show()