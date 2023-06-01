import numpy as np
from qutip import * 
import qutip as q
import matplotlib.pyplot as plt
from matplotlib import cm
import States 

N = 20 #why not : Hilbert Space dimensions
# xvec = np.linspace(-6, 6, 1500)

# fock = q.fock(N, 4)
# W_fock = q.wigner(fock, xvec, xvec)
# q.plot_wigner_fock_distribution(fock, figsize=(10, 4), colorbar=True)

# States.plot_wigner3d(fock,xvec)
# plt.show()


psi = (basis(10, 0) + basis(10, 3) + basis(10, 9)).unit()
xvec = np.linspace(-5, 5, 500)
States.plot_wigner3d(psi,xvec) 
plt.show()