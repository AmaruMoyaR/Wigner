import numpy as np
from qutip import * 
import qutip as q
import matplotlib.pyplot as plt
import matplotlib

# xvec = np.linspace(-6, 6, 1500)
def wigner_cmap(W, levels=1024, shift=0, max_color='#09224F',
                mid_color='#FFFFFF', min_color='#530017',
                neg_color='#FF97D4', invert=False):
    """A custom colormap that emphasizes negative values by creating a
    nonlinear colormap.

    Parameters
    ----------
    W : array
        Wigner function array, or any array.
    levels : int
        Number of color levels to create.
    shift : float
        Shifts the value at which Wigner elements are emphasized.
        This parameter should typically be negative and small (i.e -1e-5).
    max_color : str
        String for color corresponding to maximum value of data.  Accepts
        any string format compatible with the Matplotlib.colors.ColorConverter.
    mid_color : str
        Color corresponding to zero values.  Accepts any string format
        compatible with the Matplotlib.colors.ColorConverter.
    min_color : str
        Color corresponding to minimum data values.  Accepts any string format
        compatible with the Matplotlib.colors.ColorConverter.
    neg_color : str
        Color that starts highlighting negative values.  Accepts any string
        format compatible with the Matplotlib.colors.ColorConverter.
    invert : bool
        Invert the color scheme for negative values so that smaller negative
        values have darker color.

    Returns
    -------
    Returns a Matplotlib colormap instance for use in plotting.

    Notes
    -----
    The 'shift' parameter allows you to vary where the colormap begins
    to highlight negative colors. This is beneficial in cases where there
    are small negative Wigner elements due to numerical round-off and/or
    truncation.

    """
    cc = matplotlib.colors.ColorConverter()
    max_color = np.array(cc.to_rgba(max_color), dtype=float)
    mid_color = np.array(cc.to_rgba(mid_color), dtype=float)
    if invert:
        min_color = np.array(cc.to_rgba(neg_color), dtype=float)
        neg_color = np.array(cc.to_rgba(min_color), dtype=float)
    else:
        min_color = np.array(cc.to_rgba(min_color), dtype=float)
        neg_color = np.array(cc.to_rgba(neg_color), dtype=float)
    # get min and max values from Wigner function
    
    # nrm = matplotlib.colors.Normalize(-W.max(), W.max())
    bounds = [-W.max(), W.max()]
    # create empty array for RGBA colors
    adjust_RGBA = np.hstack((np.zeros((levels, 3)), np.ones((levels, 1))))
    zero_pos = int(np.round(levels * np.abs(shift - bounds[0])
                        / (bounds[1] - bounds[0])))
    num_pos = levels - zero_pos
    num_neg = zero_pos - 1
    # set zero values to mid_color
    adjust_RGBA[zero_pos] = mid_color
    # interpolate colors
    for k in range(0, levels):
        if k < zero_pos:
            interp = k / (num_neg + 1.0)
            adjust_RGBA[k][0:3] = (1.0 - interp) * \
                min_color[0:3] + interp * neg_color[0:3]
        elif k > zero_pos:
            interp = (k - zero_pos) / (num_pos + 1.0)
            adjust_RGBA[k][0:3] = (1.0 - interp) * \
                mid_color[0:3] + interp * max_color[0:3]
    # create colormap
    wig_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('wigner_cmap',
                                                            adjust_RGBA,
                                                            N=levels)
    return wig_cmap


def plot_wigner3d(psi,xvec):
    W = wigner(psi, xvec, xvec)
    wmap = wigner_cmap(W)  # Generate Wigner colormap
    
    fig = plt.figure(figsize=(16, 8))

    ax = fig.add_subplot(1, 2, 1)        
    plot_wigner(psi, fig=fig, ax=ax, alpha_max=5.5, cmap=wmap, colorbar = True);

    ax = fig.add_subplot(1, 2, 2, projection='3d')
    plot_wigner(psi, fig=fig, ax=ax, projection='3d', alpha_max=5.5, cmap=wmap);
    # plt.close(fig)
    return fig


def figuremaker(psi_list, xvec, n=None, m =None):
    
    W_list = [wigner(psi, xvec, xvec) for psi in psi_list] #list of states
    wmap_list = [wigner_cmap(W) for W in W_list] #map for colors of 
    # n   # Number of columns ; m   # Number of rows
    fig, axes = plt.subplots(n, m, figsize=(6, 8))  # Create the figure and subplots
    fig.patch.set_visible(False)
    # 
    # Iterate over the rows and columns
    for i in range(n):
        for j in range(m):
            # Customize the plotting based on the position of the subplot
            if j == 0:
                # for w in range(len(psi_list)):
                ax = axes[i, j]
                ax.axis('off')
                ax = fig.add_subplot(n,m,i*m+j+1)   
                plot_wigner(psi_list[i], fig=fig, ax=ax, alpha_max=5.5, cmap=wmap_list[i], colorbar = True);
                # 
                # fig.tight_layout()
            
            elif j == 1:
                # for w in range(len(psi_list)):
                ax = axes[i, j]
                ax.axis('off')
                ax = fig.add_subplot(n, m, i*m+j+1, projection='3d')
                plot_wigner(psi_list[i], fig=fig, ax=ax, projection='3d', alpha_max=5.5, cmap=wmap_list[i]);
                # ax.axis('off')
                # fig.tight_layout()
    # plt.close(fig)
    # Adjust the spacing between subplots
    fig.tight_layout()
    return fig
    
    
def Hamiltonian_2levels(psi_inic, t):
    #Constantes del sistema
    w0 = 1
    w = 1
    lamda = 1
    #Creamos los operadores para el hamiltoniano
    sig_atomo = q.tensor(q.destroy(2),q.identity(N)) #producto tensorial entre sigma 2 y la identidad N
    a_campo = q.tensor(q.identity(2),q.destroy(N)) #producto tensorial entre la identidad 2 y operador destruccion N
    sig_z = q.tensor(q.sigmaz(),q.identity(N)) #producto tensorial entre sigma z y la identidad N
    #Creamos el hamiltoniano para la evolución temporal
    H_int = 0.5 * w0 * sig_z + w*a_campo.dag()*a_campo + lamda*(sig_atomo.dag()*a_campo + sig_atomo*a_campo.dag())
    
    #estado basal de 2 niveles
    g_atomo = q.basis(2,0)#ground state
    e_atomo = q.basis(2,1)#excited state
    estado_e = q.tensor(e_atomo,psi_inic)
    estado_g = q.tensor(g_atomo,psi_inic)
    
    tiempo = np.linspace(0,25*lamda/w,300)
    #Dado el estado inicial estado_g, la evolución se calcula usando mesolve
    estado_final = q.mesolve(H_int,estado_g,tiempo)
    evolucion_temporal_estado = estado_final.states
    return evolucion_temporal_estado, tiempo
    
def PopulationInv(sig_z,estado): #calcula tasa de inversion
    tasa_inversion = q.expect(sig_z,estado)
    return tasa_inversion
    
def WignerEvolution(State_in_time,xvec,Wigners):
    for i in range(len(State_in_time)):
        Wigners.append(q.wigner(State_in_time[i],xvec,xvec))    
    return Wigners

def VonEntropy(State,t):
    rho = q.ket2dm(State)
    Entropy = [q.entropy_vn(rho) for i in range(len(t))]
    return Entropy 

#     def plot_wigner3d(psi_list, xvec, fig_vertical=None, subplot_positions=None):
   

# # Create a new figure with vertical subplots
# fig_vertical = plt.figure(figsize=(8, 12))

# # Example values for psi and xvec
# psi_list = [psi1, psi2, psi3]  # Replace with your psi values
# xvec = np.linspace(-5, 5, 100)  # Replace with appropriate xvec values

# # Define the subplot positions for each subplot
# subplot_positions = ['311', '312', '313']

# # Call plot_wigner3d to add subplots to the figure
# fig_vertical = plot_wigner3d(psi_list, xvec, fig_vertical, subplot_positions)

# # Show the figure with subplots
# plt.show()

# # Create a new figure with vertical subplots
# fig_vertical = plt.figure(figsize=(8, 12))

# # Call plot_wigner3d for the first subplot
# fig_vertical = plot_wigner3d(psi, xvec, fig_vertical, 311)

# # Call plot_wigner3d for the second subplot
# fig_vertical = plot_wigner3d(psi, xvec, fig_vertical, 312)

# # Call plot_wigner3d for the third subplot
# fig_vertical = plot_wigner3d(psi, xvec, fig_vertical, 313)

# # Show the figure with subplots
# plt.show()
# #############

# psi_list = [(basis(10, 8) + basis(10,3)+ basis(10,5)).unit(), q.fock(N, 3), q.coherent(N, np.sqrt(10)) + q.coherent(N, -np.sqrt(10))]

# fig_vertical = plt.figure(figsize=(8, 12))

# num_subplots = len(psi_list)
# W_list = [wigner(psi, xvec, xvec) for psi in psi_list]
# wmap_list = [wigner_cmap(W) for W in W_list]
# subplot_positions = [311, 312, 313]

# ax = fig_vertical.add_subplot(3,1,1)
# plot_wigner3d(psi_list[0],xvec)

# ax = fig_vertical.add_subplot(3,1,2)
# plot_wigner3d(psi_list[1], xvec)

# ax = fig_vertical.add_subplot(3,1,3)
# plot_wigner3d(psi_list[2], xvec)

# # Call plot_wigner3d to add subplots to the figure
# # fig_vertical = plot_wigner3d(psi_list, xvec, fig_vertical, subplot_positions)

# # Show the figure with subplots
# plt.show()