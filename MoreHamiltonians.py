# %%
import numpy as np
from qutip import * 
# import qutip as q
import matplotlib.pyplot as plt
from States import *
 
N = 40
#quantum harmonic oscillator with a squeezed vacuum state as initial state
w = 1 * 2 * np.pi              # oscillator frequency
tlist = np.linspace(0, 4, 101) # periods

a = destroy(N)
n = num(N)
x = (a + a.dag())/np.sqrt(2)
p = -1j * (a - a.dag())/np.sqrt(2)

# the quantum harmonic oscillator Hamiltonian
H = w * a.dag() * a
c_ops = []
# uncomment to see how things change when disspation is included
# c_ops = [np.sqrt(0.25) * a]

def plot_expect_with_variance(N, op_list, op_title, states):
    """
    Plot the expectation value of an operator (list of operators)
    with an envelope that describes the operators variance.
    """
    
    fig, axes = plt.subplots(1, len(op_list), figsize=(14,3))

    for idx, op in enumerate(op_list):
        
        e_op = expect(op, states)
        v_op = variance(op, states)

        axes[idx].fill_between(tlist, e_op - np.sqrt(v_op), e_op + np.sqrt(v_op), color="green", alpha=0.5);
        axes[idx].plot(tlist, e_op, label="expectation")
        axes[idx].set_xlabel('Time')
        axes[idx].set_title(op_title[idx])

    return fig, axes
# %%
psi0 = squeeze(N, 1.0) * basis(N, 0)
result = mesolve(H, psi0, tlist, c_ops, [])
plot_expect_with_variance(N, [n, x, p], [r'$n$', r'$x$', r'p'], result.states)
# %%


ev_t = result.states
tasa_inversion = q.expect([n,x,p], ev_t)
# Plot_Population(tlist, tasa_inversion, 'teal')

Wigner = []
xvec = np.linspace(-6, 6, 500)
WignerEvolution(ev_t,xvec,Wigner)
# %%
AnimatedWigner(Wigner, xvec, 'wignerharmosc')
# %%
