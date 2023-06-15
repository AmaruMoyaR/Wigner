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

if __name__ == '__main__':
    multiprocessing.freeze_support()
    # Rest of your code

    w_cav = 1.0 * 2 * np.pi  # cavity frequency
    w_at = 1.0 * 2 * np.pi  # atom frequency
    g = 0.05 * 2 * np.pi  # coupling strength

    kappa = 0.005  # cavity dissipation rate
    gamma = 0.05  # atom dissipation rate
    N = 15  # number of cavity fock states
    n_th_a = 0.0  # avg number of thermal bath excitation
    use_rwa = True

    tlist = np.linspace(0, 25, 150)


    # intial state
    psi0 = tensor(basis(N, 0), basis(2, 1))  # start with an excited atom

    # operators
    a = tensor(destroy(N), qeye(2))
    sm = tensor(qeye(N), destroy(2))

    # Hamiltonian
    if use_rwa:
        H = w_cav * a.dag() * a + w_at * sm.dag() * sm + \
            g * (a.dag() * sm + a * sm.dag())
    else:
        H = w_cav * a.dag() * a + w_at * sm.dag() * sm + \
            g * (a.dag() + a) * (sm + sm.dag())


    #create a list of operators that describe the 
    # dissipation
    c_ops = []
    # cavity relaxation
    rate = kappa * (1 + n_th_a)
    if rate > 0.0:
        c_ops.append(np.sqrt(rate) * a)
    # cavity excitation, if temperature > 0
    rate = kappa * n_th_a
    if rate > 0.0:
        c_ops.append(np.sqrt(rate) * a.dag())
    # qubit relaxation
    rate = gamma
    if rate > 0.0:
        c_ops.append(np.sqrt(rate) * sm)

    output = mesolve(H, psi0, tlist, c_ops, [])
    #the fith parameter allows the calculation of expectation values for the
    # operators in the list.
    rho = output.states
    rho_zz = steadystate(H, c_ops)

    xvec = np.linspace(-7, 7, 150)

    # Parallelized Wigner computation
    def compute_wigner(rho):
        return wigner(rho, xvec, xvec)

    W_list = parfor(compute_wigner, rho)

    # displaying animation with auto play and looping
    fname = "rabi_oscill.gif"
    fig, ax = plt.subplots(1,2, figsize=(12,5))

    # Plot steady state
    ax[1].contourf(xvec, xvec, wigner(rho_zz, xvec, xvec), 100, norm=clr.Normalize(-0.25,0.25), cmap=plt.get_cmap("RdBu"))
    ax[1].set_aspect("equal"); ax[1].set_title("Steady State", fontsize=14); ax[1].tick_params(labelsize=14)
    ax[1].set_xlabel(r"$x$", fontsize=14); ax[1].set_ylabel(r"$p$", fontsize=14)

    # Animate evolution
    def animate(n):
        ax[0].cla(); ax[0].set_aspect("equal"); ax[0].tick_params(labelsize=14)
        ax[0].set_title("Time: %.2f"%(output.times[n]), fontsize=14);
        ax[0].set_xlabel(r"$x$", fontsize=14); ax[0].set_ylabel(r"$p$", fontsize=14)
        im = ax[0].contourf(xvec, xvec, W_list[n], 100, norm=clr.Normalize(-0.25,0.25), cmap=plt.get_cmap("RdBu"))
        def setvisible(self, vis): # Work-around for visibility bug in contourf
            for c in self.collections: c.set_visible(vis)
        im.set_visible = types.MethodType(setvisible, im)
    anim = ani.FuncAnimation(fig, animate, frames=len(output.times))
    anim.save(fname, writer="imagemagick", fps=20)
    plt.show()


# plt.close()
# for idx, rho_ss in enumerate(rho):
#     # trace out the cavity density matrix
#     rho_cavity = ptrace(rho_ss, 0)

#     # calculate its wigner function
#     W = wigner(rho_cavity, xvec, xvec)
#     W_list = parfor(W, output.states)
#     # plot its wigner function
#     wlim = abs(W).max()
#     axes[0, idx].contourf(
#         xvec,
#         xvec,
#         W,
#         100,
#         norm = mpl.colors.Normalize(-wlim, wlim),
#         cmap = plt.get_cmap("RdBu"),
#     )
#     axes[0, idx].set_title(r"$t = %.1f$" % tlist[idx])

#     # plot its fock-state distribution
#     axes[1, idx].bar(np.arange(0, N), np.real(rho_cavity.diag()),
#                      color="blue", alpha=0.8)
#     axes[1, idx].set_ylim(0, 1)
#     axes[1, idx].set_xlim(0, 15)

# plt.show()









# N = 40 #we define the Hilbert Space Dimensions # number of cavity fock states
# n = np.arange(0, N-5)
# w_atom = 1.0 * 2 * np.pi #atom frequency
# w_cavity = 1.0 * 2 * np.pi
# delta_w = w_atom - w_cavity

# g = 0.05 * 2 * np.pi  # coupling strength

# kappa = 0.04  # cavity dissipation rate
# gamma = 0.00  # atom dissipation rate
# Gamma = 0.35  # atom pump rate
# use_rwa = False # idk what this does, I trust!
# # N = 50  
# n_th_a = 0.0  # avg number of thermal bath excitation
# tlist = np.linspace(0, 150, 101) #time steps!

# psi0 = tensor(basis(N, 0), basis(2, 0))  # start without excitations
# xvec = np.linspace(-5, 5, 200)

# # operators
# a = tensor(destroy(N), identity(2)) #anihilation x identity ##a_a
# sm = tensor(identity(N), destroy(2)) #identity x anihilation ##a_b
# sz = tensor(identity(N), sigmaz()) #sigma 
# # Hamiltonian
# # H = w_cavity * a.dag() * a + w_atom * sm.dag() * sm + g * (a.dag() + a) * sx
# # Hamiltonian
# if use_rwa:
#     H = w_cavity * a.dag() * a + w_atom * sm.dag() * sm + g * (a.dag() * sm + a * sm.dag())
# else:
#     # H = w_cavity * a.dag() * a + w_atom * sm.dag() * sm + g * (a.dag() + a) * (sm + sm.dag())
#     H = 1/2 * w_cavity * sz + w_atom*a.dag()*a + Gamma*(sm.dag()*a + sm*a.dag())

# # collapse operators
# c_ops = []
# # cavity relaxation
# rate = kappa * (1 + n_th_a)
# if rate > 0.0:
#     c_ops.append(np.sqrt(rate) * a)
# # cavity excitation, if temperature > 0
# rate = kappa * n_th_a
# if rate > 0.0:
#     c_ops.append(np.sqrt(rate) * a.dag())
# # qubit relaxation
# rate = Gamma
# if rate > 0.0:
#     c_ops.append(np.sqrt(rate) * sm.dag())

# opt = Options(nsteps=2000)  # allow extra time-steps
# #idk what it does imma keep it for now! ^

# #evolve the system
# output = mesolve(H, psi0, tlist, c_ops, [a.dag() * a, sm.dag() * sm],
#                  options=opt)


# n_c = output.expect[0] #we can get the occupation pbb!
# n_a = output.expect[1]

# e_i = coherent(N,np.sqrt(5)) #estado coherente con 5 fotones promedio
# #estado basal de 2 niveles
# #ground state
# g_atomo = basis(2,0)
# #excited state
# e_atomo = basis(2,1)

# #Estado inicial es el producto tensorial
# estado_e = tensor(e_atomo,e_i)
# estado_g = tensor(g_atomo,e_i)

# #Calculamos la probabilidad atómica para ambos niveles de energía según el tiempo
# #definimos el tiempo
# tiempo = np.linspace(0,25*Gamma/w_atom,500)

# #Dado el estado inicial estado_g, la evolución se calcula usando mesolve
# estado_final = mesolve(H,estado_g,tiempo)

# evolucion_temporal_estado = estado_final.states
# #QUeremos ver la tasa de inversion definida como el bracket del estado con el operador sig_z
# tasa_inversion = expect(sz,evolucion_temporal_estado)

# # tlist = np.linspace(0, 25, 5) #list of time
# # output = mesolve(H, psi0, tlist, c_ops, [],
# #                  options=Options(nsteps=5000))
# # rho = output.states #density matix yeaaah
# # print(rho)
# # fig, axes = plt.subplots(2, len(rho), figsize=(3 * len(rho), 6))
# wigners=[]
# x_vec = np.linspace(-10,10,300)

# for i in range(len(evolucion_temporal_estado)):
#     wigners.append(wigner(evolucion_temporal_estado[i],x_vec,x_vec))


# from moviepy.editor import VideoClip
# from moviepy.video.io.bindings import mplfig_to_npimage
 
# #duracion del video
# duration = 25*Gamma/w_atom
# x_vec = np.linspace(-10,10,300)
# fig,ax = plt.subplots()
# # method to get frames
# # t recorre desde 0 hasta duration.
# def make_frame(t):
     
#     # clear
#     ax.clear()
#     #cada tiempo se define como
#     tiempo = int(t*500/duration)
#     # Now we plot the energy spectrum for inter-cell coupling 
#     cm1 = ax.pcolormesh(x_vec,x_vec,wigners[tiempo],cmap='RdBu')
#     ax.set_title('Wigner evolution')
#     # returning numpy image
#     return mplfig_to_npimage(fig)

# # fig.colorbar(ax.pcolormesh(x_vec,x_vec,wigners[0],cmap='RdBu'))
# # creating animation
# animation = VideoClip(make_frame, duration = duration)
 