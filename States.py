import numpy as np
from qutip import * 
import qutip as q
import matplotlib.pyplot as plt
import matplotlib

# xvec = np.linspace(-6, 6, 1500)

# psi = (basis(10, 0) + basis(10, 2) + basis(10, 3)).unit() #example state

def plot_wigner3d(psi,xvec):
    #fig, axes = plt.subplots(1, 2, subplot_kw={'projection': '3d'}, figsize=(12, 6))
    
    W = wigner(psi, xvec, xvec)

    wmap = wigner_cmap(W)  # Generate Wigner colormap

    # nrm = matplotlib.colors.Normalize(-W.max(), W.max())
    
    fig = plt.figure(figsize=(20, 8))

    ax = fig.add_subplot(1, 2, 1)        
    plot_wigner(psi, fig=fig, ax=ax, alpha_max=6, cmap=wmap, colorbar = True);

    ax = fig.add_subplot(1, 2, 2, projection='3d')
    plot_wigner(psi, fig=fig, ax=ax, projection='3d', alpha_max=6, cmap=wmap);
    
    # plt.close(fig)
    return fig

# plt.show()
