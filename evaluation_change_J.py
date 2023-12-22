# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 16:25:29 2023

@author: olivi
"""

import os
os.chdir('C:/Users/olivi/OneDrive/Bureaublad/BEP')
import BEP_code_correct_inclHb
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

#Important to note: this file changes a lot over time if I want to 
#plot different functions as a function of J. In this case I used this file
#to plot the first peak as a funciton of J

def find_first_peak_index(array):
    # Find indices of relative extrema (peaks and valleys)
    extrema_indices = argrelextrema(array, np.greater)[0]
    
    # Find the index of the first peak
    return extrema_indices[0] if len(extrema_indices) > 0 else None


#allthough I've imported the main file above, I still define the parameters
#again, so I can use them to e.g. make the titles of the plots more accurate.
#Besides it is more convenient, because I don't need to take the paramters
#of the other evaluation files into consideration.
N = 20
Nb = 20
Nc = N
m = N
omega_b = 2
omega_c = 2
g = 1

N_datapoints = 46
J_largest = (N_datapoints-1)/3
t_end = 10
t_steps = 4000
t = np.linspace(0, t_end, t_steps + 1)
charging_times_energy = np.zeros(N_datapoints)
charging_times_fluc = np.zeros(N_datapoints)
J_values = np.linspace(0, J_largest, N_datapoints)
E_t_num = np.empty(len(t))
fluc_t_num = np.empty(len(t))

for iii in range(N_datapoints):
    jjj = 0
    for t_step in t:
        E_t_num[jjj] = BEP_code_correct_inclHb.E_t(t_step, Nb, Nc, m, J_values[iii], omega_b, omega_c, g)
        fluc_t_num[jjj] = BEP_code_correct_inclHb.fluctuations(t_step, Nb, Nc, m, J_values[iii], omega_b, omega_c, g)
        jjj += 1

    energy_peak_index = find_first_peak_index(E_t_num)
    fluc_peak_index = find_first_peak_index(fluc_t_num)
    time_energy_peak = t[energy_peak_index] if energy_peak_index is not None else None
    time_fluc_peak = t[fluc_peak_index] if fluc_peak_index is not None else None

    charging_times_energy[iii] = time_energy_peak
    charging_times_fluc[iii] = time_fluc_peak

# Plotting
plt.figure(1)
plt.title('Charging Time to First Peak of $E$ as a function of $J$, $ 0\leq J \leq {}$, $g={}$ $\omega=2$, $N_b=10$ $m=N_c={}$'.format(J_largest, g, N))
plt.plot(J_values, charging_times_energy, label='Energy')
plt.xlabel('$J$')
plt.ylabel('Charging Time to First Peak')
plt.legend()
plt.savefig('E_charge_time_smallbeat_Nb=10_Ncm=20.pdf',bbox_inches='tight')
plt.show()

plt.figure(2)
plt.title('Charging Time to First Peak of $\Sigma^2$ as a function of $J$, $0 \leq J \leq {}$, $g={}$, $\omega=2$, $N_b=10$ $m=N_c={}$'.format(J_largest, g, N))
plt.plot(J_values, charging_times_fluc, label='Fluctuations')
plt.xlabel('$J$')
plt.ylabel('Charging Time to First Peak')
plt.legend()
plt.savefig('f_charge_time_smallbeat_Nb=10_Ncm=20.pdf',bbox_inches='tight')
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
