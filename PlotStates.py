import numpy as np
from qutip import * 
import qutip as q
import matplotlib.pyplot as plt
from matplotlib import cm
import States

N = 40 #why not : Hilbert Space dimensions
# xvec = np.linspace(-6, 6, 1500)
# fock = q.fock(N, 4)
# W_fock = q.wigner(fock, xvec, xvec)
# q.plot_wigner_fock_distribution(fock, figsize=(10, 4), colorbar=True)

# States.plot_wigner3d(fock,xvec)
# plt.show()


psi = [(basis(10, 8) + basis(10,3)+ basis(10,5)).unit(), q.coherent(N, np.sqrt(10)) + q.coherent(N, -np.sqrt(10)), q.fock(N, 3) ]
xvec = np.linspace(-5, 5, 1500)
# States.plot_wigner3d(psi,xvec) 
# plt.show()
# States.figuremaker(psi,xvec,3,2)
# plt.show()



psi2 = [((coherent_dm(N, -2.0) + coherent_dm(N, 2.0)).unit()), q.basis(N,10), (ket2dm(squeeze(N, 0.75j) * basis(N, 0)) +  ket2dm(squeeze(N, -0.75j) * basis(N, 0))).unit()]

States.figuremaker(psi2,xvec,3,2)
plt.show()
