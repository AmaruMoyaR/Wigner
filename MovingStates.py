#%%

import numpy as np
from qutip import * 
# import qutip as q
import matplotlib.pyplot as plt
from States import *
N = 40

# %%
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
psi_inicial_cat = (coherent(N, -3.0) + coherent(N, 3.0)) / np.sqrt(2)
evolucion_temporal_estado_cat, tiempo, tasa_inversion = Hamiltonian_2levels(psi_inicial_cat, N , 300)

wignerslist = []
xvec = np.linspace(-10,10,300)

Wigners = WignerEvolution(evolucion_temporal_estado_cat, xvec, wignerslist)

#%%
AnimatedWigner(Wigners,xvec,'wignercat')

AnimatePopulation(tasa_inversion,tiempo,'populationwignercat')

# %%
Entropy = []
for i in range(len(evolucion_temporal_estado_cat)):
    Entropy.append(entropy_vn(evolucion_temporal_estado_cat[i]))

fig, ax = plt.subplots(figsize=(6, 6))    
ax.plot(tiempo,Entropy , color = 'firebrick')
ax.set_xlabel('Tiempo')
ax.set_ylabel('Entropia')
ax.set_title('Entropia de Von Neumann')

# %%
psi_inicial_Sq = q.displace(N, 2) * q.squeeze(N, 1.0)* q.basis(N,0) 
evolucion_temporal_estado_Sq, tiempo, tasa_inversion = Hamiltonian_2levels(psi_inicial_Sq, N , 300)

wignerslist = []
xvec = np.linspace(-10,10,300)

Wigners = WignerEvolution(evolucion_temporal_estado_Sq, xvec, wignerslist)


# %%
AnimatedWigner(Wigners,xvec,'wignersq')

AnimatePopulation(tasa_inversion,tiempo,'populationwignersq')

# %%
Entropy = []
for i in range(len(evolucion_temporal_estado_Sq)):
    Entropy.append(entropy_vn(evolucion_temporal_estado_Sq[i]))

fig, ax = plt.subplots(figsize=(6, 6))    
ax.plot(tiempo,Entropy , color = 'firebrick')
ax.set_xlabel('Tiempo')
ax.set_ylabel('Entropia')
ax.set_title('Entropia de Von Neumann')
# %%
