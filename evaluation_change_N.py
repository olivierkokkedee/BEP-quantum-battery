# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 16:46:37 2023

@author: olivi
"""

import os

os.chdir('C:/Users/private/Bureaublad/BEP')
import BEP_code_correct_inclHb
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

#this file is quite similar to "evaulation_change_J.py", so please
#take a look at that file if something is unclear.

def fit_fluc(x,a,b):
    return a*(1-np.cos(b/x))

def fit_E(x,b):
    return 2.5*np.cos(b/x)
Nb=5

h=1
J=h
omega_b=2
omega_c=2
g=h
ftsize=16
N_datapoints=21
N_start=5
N_end=41
N_values=np.linspace(N_start,N_end,N_datapoints)
'''
max_values_fluc=np.zeros(N_datapoints)
max_values_energy=np.zeros(N_datapoints)
E_t_num=np.empty(len(t))
fluc_t_num=np.empty(len(t))
for iii in range(N_datapoints):
    jjj=0
    for t_step in t:
        E_t_num[jjj]=E_t(t_step, iii+2, iii+2, iii+2, J, omega_b, omega_c, g)
        fluc_t_num[jjj]=fluctuations(t_step, iii+2, iii+2, iii+2, J, omega_b, omega_c, g)
        jjj+=1
    max_values_energy[iii]=np.amax(E_t_num)
    max_values_fluc[iii]=np.amax(fluc_t_num)
    
    '''
t_end=10
t_steps=8000
t=np.linspace(0,t_end,t_steps)

'''
for jjj in range(N_graphs):
    E_temp=np.empty(len(t))
    fluc_temp=np.empty(len(t))
    max_fluc_temp=np.zeros(N_datapoints)
    max_E_temp=np.zeros(N_datapoints)
    for iii in range(N_datapoints):
        kkk=0
        for t_step in t:
            E_temp[kkk]=BEP_code_correct_inclHb.E_t(t_step, Nb, iii+Nb, iii+Nb, 0.5*jjj, omega_b, omega_c, g)
            fluc_temp[kkk]=BEP_code_correct_inclHb.fluctuations(t_step, Nb, iii+Nb, iii+Nb, 0.5*jjj, omega_b, omega_c, g)
            kkk+=1
        max_E_temp[iii]=np.amax(E_temp)
        max_fluc_temp[iii]=np.amax(fluc_temp)
    E_graphs.append(max_E_temp)
    fluc_graphs.append(max_fluc_temp)
'''

E_num=np.empty(len(t))
fluc_num=np.empty(len(t))
E_max=np.zeros(N_datapoints)
fluc_at_max_E=np.zeros(N_datapoints)
for iii in range(N_datapoints):
    kkk=0
    for t_step in t:
        E_num[kkk]=BEP_code_correct_inclHb.E_t(t_step, Nb, int(N_values[iii]*2), int(N_values[iii]), h, omega_b, omega_c, h)
        fluc_num[kkk]=BEP_code_correct_inclHb.fluctuations(t_step, Nb, int(N_values[iii]*2), int(N_values[iii]), h, omega_b, omega_c, h)
        kkk+=1
    index_max_E=np.argmax(E_num)
    print(index_max_E)
    E_max[iii]=E_num[index_max_E]
    fluc_at_max_E[iii]=fluc_num[index_max_E]
    


E_fit_values,cov=curve_fit(fit_E,N_values,E_max)
fluc_fit_values,cov=curve_fit(fit_fluc,N_values,fluc_at_max_E)

N_continuous=np.linspace(N_start,N_end,1000)
plt.figure(1)
plt.plot(N_continuous,fit_E(N_continuous,*E_fit_values),color='red',label='Fit')
plt.plot(N_values,E_max,'b.',label='Found values')
plt.plot(N_values,np.full(len(N_values),omega_b*Nb/4),linestyle='--',label='$E_{max}$')
plt.xlabel('$m$', fontsize=ftsize)
plt.ylabel('$E_b(t)\;\; [\hbar \omega]$', fontsize=ftsize)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.legend(fontsize=12)
plt.grid()
plt.tight_layout()
plt.savefig('E_max_Nb=5_h=1_fitat2_5.pdf')
plt.show()


plt.figure(2)
plt.plot(N_continuous,fit_fluc(N_continuous,*fluc_fit_values),color='red',label='Fit')
plt.plot(N_values,fluc_at_max_E,'b.',label='Found values')
plt.plot(N_values,np.full(len(N_values),0),linestyle='--')
plt.xlabel('$m$', fontsize=ftsize)
plt.ylabel('$\Sigma^2_b(t)\;\; [(\hbar \omega)^2]$', fontsize=ftsize)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.legend(fontsize=12)
plt.grid()
plt.tight_layout()
plt.savefig('fluc_at_E_max_Nb=5_h=1.pdf')
plt.show()


'''
max_values_fluc_constNc=np.zeros(N_datapoints)
max_values_energy_constNc=np.zeros(N_datapoints)
N_values=np.linspace(0,N_datapoints-1,N_datapoints)
E_t_num_constNc=np.empty(len(t))
fluc_t_num_constNc=np.empty(len(t))
for iii in range(N_datapoints):
    jjj=0
    for t_step in t:
        E_t_num_constNc[jjj]=E_t(t_step, iii+2, iii+4, iii+3, J, omega_b, omega_c, g)
        fluc_t_num[jjj]=fluctuations(t_step, iii+2, iii+4, iii+3, J, omega_b, omega_c, g)
        jjj+=1
    max_values_energy[iii]=np.amax(E_t_num)
    max_values_fluc[iii]=np.amax(fluc_t_num)
    '''
