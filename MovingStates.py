import numpy as np
from qutip import * 
# import qutip as q
import matplotlib.pyplot as plt
from States import plot_wigner3d
from scipy.linalg import expm
import multiprocessing

import matplotlib.animation as ani
import matplotlib.colors as clr
import types









# if __name__ == '__main__':
#     multiprocessing.freeze_support()
#     # Rest of your code

#     w_cav = 1.0 * 2 * np.pi  # cavity frequency
#     w_at = 1.0 * 2 * np.pi  # atom frequency
#     g = 0.05 * 2 * np.pi  # coupling strength

#     kappa = 0.005  # cavity dissipation rate
#     gamma = 0.05  # atom dissipation rate
#     N = 15  # number of cavity fock states
#     n_th_a = 0.0  # avg number of thermal bath excitation
#     use_rwa = True

#     tlist = np.linspace(0, 25, 150)


#     # intial state
#     psi0 = tensor(basis(N, 0), basis(2, 1))  # start with an excited atom

#     # operators
#     a = tensor(destroy(N), qeye(2))
#     sm = tensor(qeye(N), destroy(2))

#     # Hamiltonian
#     if use_rwa:
#         H = w_cav * a.dag() * a + w_at * sm.dag() * sm + \
#             g * (a.dag() * sm + a * sm.dag())
#     else:
#         H = w_cav * a.dag() * a + w_at * sm.dag() * sm + \
#             g * (a.dag() + a) * (sm + sm.dag())


#     #create a list of operators that describe the 
#     # dissipation
#     c_ops = []
#     # cavity relaxation
#     rate = kappa * (1 + n_th_a)
#     if rate > 0.0:
#         c_ops.append(np.sqrt(rate) * a)
#     # cavity excitation, if temperature > 0
#     rate = kappa * n_th_a
#     if rate > 0.0:
#         c_ops.append(np.sqrt(rate) * a.dag())
#     # qubit relaxation
#     rate = gamma
#     if rate > 0.0:
#         c_ops.append(np.sqrt(rate) * sm)

#     output = mesolve(H, psi0, tlist, c_ops, [])
#     #the fith parameter allows the calculation of expectation values for the
#     # operators in the list.
#     rho = output.states
#     rho_zz = steadystate(H, c_ops)

#     xvec = np.linspace(-7, 7, 150)

#     # Parallelized Wigner computation
#     def compute_wigner(rho):
#         return wigner(rho, xvec, xvec)

#     W_list = parfor(compute_wigner, rho)

#     # displaying animation with auto play and looping
#     fname = "rabi_oscill.gif"
#     fig, ax = plt.subplots(1,2, figsize=(12,5))

#     # Plot steady state
#     ax[1].contourf(xvec, xvec, wigner(rho_zz, xvec, xvec), 100, norm=clr.Normalize(-0.25,0.25), cmap=plt.get_cmap("RdBu"))
#     ax[1].set_aspect("equal"); ax[1].set_title("Steady State", fontsize=14); ax[1].tick_params(labelsize=14)
#     ax[1].set_xlabel(r"$x$", fontsize = 14); ax[1].set_ylabel(r"$p$", fontsize=14)

#     # Animate evolution
#     def animate(n):
#         ax[0].cla(); ax[0].set_aspect("equal"); ax[0].tick_params(labelsize=14)
#         ax[0].set_title("Time: %.2f"%(output.times[n]), fontsize=14);
#         ax[0].set_xlabel(r"$x$", fontsize=14); ax[0].set_ylabel(r"$p$", fontsize=14)
#         im = ax[0].contourf(xvec, xvec, W_list[n], 100, norm=clr.Normalize(-0.25,0.25), cmap=plt.get_cmap("RdBu"))
#         def setvisible(self, vis): # Work-around for visibility bug in contourf
#             for c in self.collections: c.set_visible(vis)
#         im.set_visible = types.MethodType(setvisible, im)
#     anim = ani.FuncAnimation(fig, animate, frames=len(output.times))
#     anim.save(fname, writer="imagemagick", fps=20)
#     plt.show()

