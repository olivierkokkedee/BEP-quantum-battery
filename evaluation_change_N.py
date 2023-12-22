# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 16:46:37 2023

@author: olivi
"""

import os

os.chdir('C:/Users/olivi/OneDrive/Bureaublad/BEP')
import BEP_code_correct_inclHb
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

def linear(x,a,b):
    return a*x+b

def quadratic(x,a,b,c):
    return a*x**2+b*x+c

Nb=4

J=1
omega_b=2
omega_c=2
g=15

N_datapoints=40
N_values=np.linspace(0,N_datapoints-1,N_datapoints)+Nb
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
t_end=5e14
t_steps=4000
t=np.linspace(0,t_end,t_steps+1)
N_graphs=5
E_graphs=[]
fluc_graphs=[]
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
    
#E_fit,lin_cov=curve_fit(linear,N_values[8:22],max_values_energy[8:22])
#sigma_fit,quad_cov=curve_fit(quadratic,N_values[8:22],max_values_fluc[8:22])
'''
plt.figure(1)
plt.title('$E(t)$ as a function of $N$, $0\leq N\leq 24$, $g={}$, J={}'.format(g,J))
plt.plot(N_values,max_values_energy)
#plt.plot(N_values,linear(N_values,*E_fit),'g--',label='a={}, b={}'.format(*list(np.around(np.array(E_fit),2))))
plt.xlabel('$N$')
plt.ylabel('$E(t)$')
plt.legend()
#plt.savefig('E_func_of_N_J={}_g={}.pdf'.format(J,g))
plt.show()


plt.figure(2)
plt.title('$\Sigma^2$ as a function of $N$, $2\leq N_b\leq 24$,$N_c=m=N_b$, $g=\omega =J=2$'.format())
plt.plot(N_values,max_values_fluc)
#plt.plot(N_values,quadratic(N_values,*sigma_fit),'g--',label='a={}, b={}, c={}'.format(*list(np.around(np.array(sigma_fit),2))))
plt.xlabel('$N$')
plt.ylabel('$\Sigma^2(t)$')
#plt.savefig('fluc_func_of_N_J={}_g={}.pdf'.format(J,g))
plt.legend()
plt.show()
'''
plt.figure(1)
plt.title('$E$ as a function of $N$, ${}\leq N_c,m\leq {}$, $N_b={}$, $g={}$'.format(Nb,Nb+N_datapoints-1,Nb,g))
plt.plot(N_values,E_graphs[0],label='J=0')
plt.plot(N_values,E_graphs[1],label='J=0.5')
plt.plot(N_values,E_graphs[2],label='J=1')
plt.plot(N_values,E_graphs[3],label='J=1.5')
plt.plot(N_values,E_graphs[4],label='J=2')
plt.xlabel('$N$')
plt.ylabel('$E$')
plt.legend()
plt.savefig('E_func_N_Nb=4_inclJ=0_g=15.pdf')
plt.show()


plt.figure(2)
plt.title('$\Sigma^2$ as a function of $N$, ${}\leq N_c,m\leq {}$, $N_b={}$, $g={}$'.format(Nb,Nb+N_datapoints-1,Nb,g))
plt.plot(N_values,fluc_graphs[0],label='J=0')
plt.plot(N_values,fluc_graphs[1],label='J=0.5')
plt.plot(N_values,fluc_graphs[2],label='J=1')
plt.plot(N_values,fluc_graphs[3],label='J=1.5')
plt.plot(N_values,fluc_graphs[4],label='J=2')
plt.xlabel('$N$')
plt.ylabel('$\Sigma^2$')
plt.savefig('fluc_func_of_N_Nb=4_inclJ=0_g=15.pdf')
plt.legend()
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