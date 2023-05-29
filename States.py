import numpy as np
from qutip import * 
import qutip as q
import matplotlib.pyplot as plt
from matplotlib import cm


def plot_wigner_2d_3d(psi,xvec):
    #fig, axes = plt.subplots(1, 2, subplot_kw={'projection': '3d'}, figsize=(12, 6))
    W = wigner(psi, xvec, xvec)

    wmap = wigner_cmap(W)  # Generate Wigner colormap
    
    fig = plt.figure(figsize=(20, 8))

    ax = fig.add_subplot(1, 2, 1)        
    plot_wigner(psi, fig=fig, ax=ax, alpha_max=6, cmap=wmap, colorbar = True);

    ax = fig.add_subplot(1, 2, 2, projection='3d')
    plot_wigner(psi, fig=fig, ax=ax, projection='3d', alpha_max=6, cmap=wmap);
    
    # plt.close(fig)
    return fig

# plt.show()
