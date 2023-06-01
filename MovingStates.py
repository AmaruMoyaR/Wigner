import numpy as np
from qutip import * 
# import qutip as q
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from States import plot_wigner3d
import matplotlib.animation as ani
import matplotlib.colors as clr
import types


N = 30 #we define the Hilbert Space Dimensions
n = np.arange(0, N-5)
w_atom = 1.0 * 2 * np.pi #atom frequency
w_cavity = 1.0 * 2 * np.pi
delta_w = w_atom - w_cavity

g = 0.05 * 2 * np.pi  # coupling strength

kappa = 0.04  # cavity dissipation rate
gamma = 0.00  # atom dissipation rate
Gamma = 0.35  # atom pump rate

N = 30  # number of cavity fock states
n_th_a = 0.0  # avg number of thermal bath excitation

tlist = np.linspace(0, 150, 101)

# intial state
psi0 = tensor(basis(N, 0), basis(2, 0))  # start without excitations

# operators
a = tensor(destroy(N), qeye(2)) #anihilation x identity ##a_a
sm = tensor(qeye(N), destroy(2)) #identity x anihilation ##a_b
sx = tensor(qeye(N), sigmax()) #sigma 

# Hamiltonian
H = w_cavity * a.dag() * a + w_atom * sm.dag() * sm + g * (a.dag() + a) * sx

# collapse operators
c_ops = []

rate = kappa * (1 + n_th_a)
if rate > 0.0:
    c_ops.append(np.sqrt(rate) * a)

rate = kappa * n_th_a
if rate > 0.0:
    c_ops.append(np.sqrt(rate) * a.dag())

rate = gamma
if rate > 0.0:
    c_ops.append(np.sqrt(rate) * sm)

rate = Gamma
if rate > 0.0:
    c_ops.append(np.sqrt(rate) * sm.dag())
    
opt = Options(nsteps=2000)  # allow extra time-steps
#evolve the system
output = mesolve(H, psi0, tlist, c_ops, [a.dag() * a, sm.dag() * sm],
                 options=opt)

tlist = np.linspace(0, 25, 5) #list of time

output = mesolve(H, psi0, tlist, c_ops, [],
                 options=Options(nsteps=5000))

rho_ss_sublist = output.states

xvec = np.linspace(-5, 5, 200)

fig, axes = plt.subplots(2, len(rho_ss_sublist),
                         figsize=(3 * len(rho_ss_sublist), 6))

for idx, rho_ss in enumerate(rho_ss_sublist):

    # trace out the cavity density matrix
    rho_ss_cavity = ptrace(rho_ss, 0)

    # calculate its wigner function
    W = wigner(rho_ss_cavity, xvec, xvec)
    # W_list = parfor(W, output.states)
    # plot its wigner function
    wlim = abs(W).max()
    axes[0, idx].contourf(
        xvec,
        xvec,
        W,
        100,
        norm = mpl.colors.Normalize(-wlim, wlim),
        cmap = plt.get_cmap("RdBu"),
    )
    axes[0, idx].set_title(r"$t = %.1f$" % tlist[idx])

    # plot its fock-state distribution
    axes[1, idx].bar(np.arange(0, N), np.real(rho_ss_cavity.diag()),
                     color="blue", alpha=0.8)
    axes[1, idx].set_ylim(0, 1)
    axes[1, idx].set_xlim(0, 15)

plt.show()

# Prepare figure
fname = "wignertime.gif"
fig, ax = plt.subplots(1,2, figsize=(12,5))

# Plot steady state
ax[1].contourf(xvec, xvec, wigner(rho_ss, xvec, xvec), 100, norm=clr.Normalize(-0.25,0.25), cmap=plt.get_cmap("RdBu"))
ax[1].set_aspect("equal");
ax[1].set_title("Steady State", fontsize=14);
ax[1].tick_params(labelsize=14)
ax[1].set_xlabel(r"$x$", fontsize=14); ax[1].set_ylabel(r"$p$", fontsize=14)

# Animate evolution
def animate(n):
    ax[0].cla(); 
    ax[0].set_aspect("equal"); 
    ax[0].tick_params(labelsize=14)
    ax[0].set_title("Time: %.2f"%(output.times[n]), fontsize=14);
    ax[0].set_xlabel(r"$x$", fontsize=14); ax[0].set_ylabel(r"$p$", fontsize=14)
    im = ax[0].contourf(xvec, xvec, W_list[n], 100, norm=clr.Normalize(-0.25,0.25), cmap=plt.get_cmap("RdBu"))
    def setvisible(self, vis): # Work-around for visibility bug in contourf
        for c in self.collections: c.set_visible(vis)
    im.set_visible = types.MethodType(setvisible, im)
anim = ani.FuncAnimation(fig, animate, frames=len(output.times))
anim.save(fname, writer="imagemagick", fps=20)
plt.show()

# import IPython.display as ipyd
# ipyd.Image(url=fname)