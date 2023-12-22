# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 18:18:36 2023

@author: olivi
"""

import os
os.chdir('C:/Users/olivi/OneDrive/Bureaublad/BEP')
import BEP_code_correct
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

#finds the first global maximum within a certain margin
def find_first_max_index(array, margin=0.05):
    max_value = np.amax(array)
    
    candidate_indices = np.where(np.isclose(array, max_value, rtol=margin))[0]

    return candidate_indices[0] if len(candidate_indices) > 0 else None

N = 20
Nb = 2
Nc = N
m = N
omega_b = 2
omega_c = 2
g = 1

N_datapoints = 30
#J_largest = (N_datapoints-1)/3
t_end = 30
t_steps = 4000
t = np.linspace(0, t_end, t_steps + 1)
charging_times_energy = np.zeros(N_datapoints)
charging_times_fluc = np.zeros(N_datapoints)
J_values = np.linspace(0, N_datapoints, N_datapoints)
E_t_num = np.empty(len(t))
fluc_t_num = np.empty(len(t))

for iii in range(N_datapoints):
    jjj = 0
    for t_step in t:
        E_t_num[jjj] = BEP_code_correct.E_t(t_step, Nb, Nc, m, J_values[iii], omega_b, omega_c, g)
        fluc_t_num[jjj] = BEP_code_correct.fluctuations(t_step, Nb, Nc, m, J_values[iii], omega_b, omega_c, g)
        jjj += 1

    energy_max_index = find_first_max_index(E_t_num)
    fluc_max_index = find_first_max_index(fluc_t_num)
    time_energy_max = t[energy_max_index] if energy_max_index is not None else None
    time_fluc_max = t[fluc_max_index] if fluc_max_index is not None else None

    charging_times_energy[iii] = time_energy_max
    charging_times_fluc[iii] = time_fluc_max

# Plotting
plt.figure(1)
plt.title('Charging Time to True First Maximum of $E(t)$ as a function of $J$, $ 0\leq J \leq {}$, $g={}$ $\omega=2$, $N_b=2$, $m=N_c=10$'.format(N_datapoints, g))
plt.plot(J_values, charging_times_energy, label='Energy')
plt.xlabel('$J$')
plt.ylabel('Charging Time to True First Maximum')
plt.legend()
plt.savefig('E_charge_time_truemax_Nb=2_Ncm=10.pdf', bbox_inches='tight')
plt.show()

plt.figure(2)
plt.title('Charging Time to True First Maximum of $\Sigma(t)^2$ as a function of $J$, $0 \leq J \leq {}$, $g={}$, $\omega=2$, $N_b=2$, $m=N_c=10$'.format(N_datapoints, g))
plt.plot(J_values, charging_times_fluc, label='Fluctuations')
plt.xlabel('$J$')
plt.ylabel('Charging Time to True First Maximum')
plt.legend()
plt.savefig('f_charge_time_truemax_Nb=2_Ncm=10.pdf', bbox_inches='tight')
plt.show()
