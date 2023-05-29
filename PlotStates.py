import numpy as np
from qutip import * 
import qutip as q
import matplotlib.pyplot as plt
from matplotlib import cm
import States 

N = 35 #why not : Hilbert Space dimensions
xvec = np.linspace(-6, 6, 1500)

psiex = (basis(10, 0) + basis(10, 3) + basis(10, 9)).unit() #example state

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


fock = q.fock(N, 3)
W_fock = q.wigner(fock, xvec, xvec)
# q.plot_wigner_fock_distribution(fock, figsize=(10, 4), colorbar=True)

States.plot_wigner3d(fock)
plt.show()