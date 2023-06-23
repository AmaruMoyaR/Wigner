#%%

import numpy as np
from qutip import * 
# import qutip as q
import matplotlib.pyplot as plt
from States import *

N = 40
psi_inicial_fock = q.basis(N,10) 
evolucion_temporal_estado_Fock, tiempo, tasa_inversion = Hamiltonian_2levels(psi_inicial_fock, N , 300)

wignerslist = []
xvec = np.linspace(-10,10,300)

Wigners = WignerEvolution(evolucion_temporal_estado_Fock, xvec, wignerslist)


# %%

AnimatedWigner(Wigners,xvec,'wignerfock')

AnimatePopulation(tasa_inversion,tiempo,'populationwignerfock')

# %%

# Entropy = VonEntropy(Wigners,tiempo)


# %%
Entropy = []
for i in range(len(evolucion_temporal_estado_Fock)):
    Entropy.append(entropy_vn(evolucion_temporal_estado_Fock[i]))

fig, ax = plt.subplots(figsize=(6, 6))    
ax.plot(tiempo,Entropy , color = 'firebrick')
ax.set_xlabel('Tiempo')
ax.set_ylabel('Entropia')
ax.set_title('Entropia de Von Neumann')


# %%
