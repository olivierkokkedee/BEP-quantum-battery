# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 16:25:29 2023

@author: olivi
"""

import os
os.chdir('private')
import BEP_code_correct_inclHb
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from scipy.optimize import curve_fit

#Important to note: this file changes a lot over time if I want to 
#plot different functions as a function of J. In this case I used this file
#to plot the first peak as a funciton of J

def find_first_peak_index(array):
    # Find indices of relative extrema (peaks and valleys)
    extrema_indices = argrelextrema(array, np.greater)[0]
    
    # Find the index of the first peak
    return extrema_indices[0] if len(extrema_indices) > 0 else None

def function_fit(x,a):
    return a/(x)

def tanh(x,a,b):
    return a*np.tan(x-b)

#allthough I've imported the main file above, I still define the parameters
#again, so I can use them to e.g. make the titles of the plots more accurate.
#Besides it is more convenient, because I don't need to take the paramters
#of the other evaluation files into consideration.
N = 20
Nb = 5
Nc = 60
m = 30
omega_b = 2
omega_c = 2
g = 1

Nb2 = 5


N_datapoints=30
t_end = 60
t_steps = 12000
t = np.linspace(0, t_end, t_steps + 1)
E_max = np.zeros(N_datapoints)
fluc_at_Emax = np.zeros(N_datapoints)
J_values = np.linspace(1, 50, N_datapoints)
E_t= np.empty(len(t))
fluc=np.empty(len(t))

for iii in range(N_datapoints):
    jjj = 0
    for t_step in t:
        E_t[jjj] = BEP_code_correct_inclHb.E_t(t_step, Nb, Nc, m, J_values[iii], omega_b, omega_c, g)
        fluc[jjj] = BEP_code_correct_inclHb.fluctuations(t_step, Nb, Nc, m, J_values[iii], omega_b, omega_c, g)
        jjj += 1
    
    energy_peak_index = np.argmax(E_t)
    E_max[iii] = E_t[energy_peak_index]
    fluc_at_Emax[iii]=fluc[energy_peak_index]
    
#E_fit_values,cov=curve_fit(tanh,J_values,E_max)
h_continuous=np.linspace(1,50,1000)
ftsize=16

# Plotting
plt.figure(1)
#plt.plot(h_continuous, tanh(h_continuous,*E_fit_values))
plt.plot(J_values,E_max ,'.',markersize=8,label='Values $N_b=4$')
plt.plot(h_continuous,np.full(len(h_continuous),omega_b*Nb/4),linestyle='--',color='orange')
plt.plot(h_continuous,np.full(len(h_continuous),-omega_b*Nb/4),linestyle='--',color='orange')
plt.axvline(23.7,color='red',linestyle='--')
plt.xlabel('$J\;\;[\hbar\omega]$', fontsize=ftsize)
plt.ylabel('$E_{max}\;\;[(\hbar\omega)]$', fontsize=ftsize)
#plt.legend()
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.grid()
plt.tight_layout()
plt.savefig('E_max_Nb={}_Nc={}_m={}_g={}_varyingJ.pdf'.format(Nb,Nc,m,g),bbox_inches='tight')
plt.show()

# Plotting
plt.figure(2)
plt.plot(J_values,fluc_at_Emax ,'.',markersize=8,label='Values $N_b=4$')
plt.plot(h_continuous,np.full(len(h_continuous),0),linestyle='--',color='orange')
plt.axvline(23.7,color='red',linestyle='--')
plt.xlabel('$J\;\;[\hbar\omega]$', fontsize=ftsize)
plt.ylabel('$\Sigma^2_{max}\;\;[(\hbar\omega)^2]$', fontsize=ftsize)
#plt.legend()
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.grid()
plt.tight_layout()
plt.savefig('fluc_max_Nb={}_Nc={}_m={}_g={}_varyingJ.pdf'.format(Nb,Nc,m,g),bbox_inches='tight')
plt.show()

'''
plt.figure(3)
plt.title('$\Sigma(t)^2$ as a function of $J$, $0 \leq J \leq 60$, $g=2$, $m=N=2$, Log scaled')
plt.plot(J_values,1-max_values_fluc/16.95)
plt.xlabel('$J$')
plt.ylabel('$\Sigma^2(t)$')
plt.yscale('log')
#plt.savefig('fluc_func_of_J_g={}_log.pdf'.format(g))
#plt.legend()
plt.show()'''
