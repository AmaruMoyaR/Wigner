# %%
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


psi = [(basis(10, 2) + basis(10,3)+ basis(10,5)).unit(), q.coherent(N, np.sqrt(10)) + q.coherent(N, -np.sqrt(10)), q.fock(N, 3) ]
xvec = np.linspace(-6, 6, 1500)
# States.plot_wigner3d(psi,xvec) 
# plt.show()
States.figuremaker(psi,xvec,3,2)
plt.show()

# %%

psi2 = [((coherent_dm(N, -2.0) + coherent_dm(N, 2.0)).unit()), q.basis(N,10), (ket2dm(squeeze(N, 0.75j) * basis(N, 0)) +  ket2dm(squeeze(N, -0.75j) * basis(N, 0))).unit()]

States.figuremaker(psi2,xvec,3,2)
plt.show()

# %%








# listentropy = []
# listentropy2 = []
# for i in range(len(psi)):
#     listentropy.append(entropy_vn(psi[i]))
#     listentropy2.append(entropy_vn(psi2[i]))

# plt.plot(listentropy, 'o')
# plt.plot(listentropy2, 'o')

# # %%

# fig, axes = plt.subplots(1, 3, figsize=(12,3))
# plot_fock_distribution(psi2[0], fig=fig, ax=axes[0], title="2 Coherent states")
# plot_fock_distribution(psi2[1], fig=fig, ax=axes[1], title="Fock state n = 10")
# plot_fock_distribution(psi2[2], fig=fig, ax=axes[2], title="Squeezed State")
# fig.tight_layout()
# plt.show()

# # %%
# fig, axes = plt.subplots(1, 3, figsize=(12,3))
# plot_fock_distribution(psi[0], fig=fig, ax=axes[0], title=" Superposition of basis")
# plot_fock_distribution(psi[1], fig=fig, ax=axes[1], title="Schrödinger cat")
# plot_fock_distribution(psi[2], fig=fig, ax=axes[2], title="Fock state n = 3")
# fig.tight_layout()
# plt.show()
# %%
