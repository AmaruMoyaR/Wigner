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


N = 50 #we define the Hilbert Space Dimensions
n = np.arange(0, N-5)
w_atom = 1.0 * 2 * np.pi #atom frequency
w_cavity = 1.0 * 2 * np.pi
#resonant!!
delta_w = w_atom - w_cavity

g = 0.05  * np.pi  # coupling strength

kappa = 0.04  # cavity dissipation rate
Gamma = 0.05  # atom dissipation rate
# Gamma = 0.35  # atom pump rate

# N = 50  # number of cavity fock states
# n_th_a = 0.0  # avg number of thermal bath excitation

alpha = 50 #avg number of photons^2 coherent state eigenvalue

# ns |α|^ 2 = n = {5, 10, 20, 50}

tlist = np.linspace(0, 150, 101)
tau_C = 2 * np.pi * np.sqrt(np.abs(alpha)**2)

# intial state
psi0 = tensor(basis(N, 0), basis(2, 0))  # start without excitations

# operators
a = tensor(destroy(N), qeye(2)) #anihilation x identity ##a_a
sm = tensor(qeye(N), destroy(2)) #identity x anihilation ##a_b
sx = tensor(qeye(N), sigmax()) #sigma 

# Hamiltonian
H = w_atom/2 * sx + w_cavity*a.dag()*a + g*(a*sm.dag() + a.dag()*sm)
 
#since they are resonant, we can get the master equation:

H_int = g*(a*sm.dag() + a.dag*sm)#Hamiltonian interaction

'''
We first consider the case of a single atom
initially in the excited state interacting with a resonant
coherent field in the absence of any decoherence mechanism. 
'''
