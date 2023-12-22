# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 14:42:35 2023

@author: olivi
"""

import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
import pandas as pd
import math
import os
N=2
Nb=2
Nc=N
m=N
J=0.333333
omega_b=2
omega_c=2
g=1

df = pd.read_excel('C:/Users/olivi/OneDrive/Bureaublad/BEP/nJ_contribution_per_state.xlsx')
nJ_contribution = df.iloc[:, 1:].to_numpy()


def u(j, Nb, Nc, m, J, omega_b, omega_c, g):
    return g * np.sqrt(j * (Nb - j + 1) * (Nc - m + j) * (m - j + 1))

def b(j, Nb, Nc, m, J, omega_b, omega_c, g):
    return omega_b / 2 * (j - Nb / 2) + omega_c / 2 * (m - j - Nc / 2)

def generate_spin_states(N, m):
    states = []
    for indices in combinations(range(N), m):
        state = [1 if i in indices else 0 for i in range(N)]
        states.append(state)
    return states

def calculate_nJ(N,m):
    nJ=0
    states=generate_spin_states(N, m)
    for lst in states:
        nJ_temp=0
        for i in range(N-1):
            if (lst[i] !=lst[i+1]):
                nJ_temp+=1
        nJ+=nJ_temp
    return nJ

def get_nJ_diagonal(Nb):
    diagonal = np.array([0])  # Start with a single element, 0

    if (Nb - 2) % 2 == 0:
        for i in range(1, int(Nb / 2) + 1):
            diagonal = np.append(diagonal, nJ_contribution[Nb - 2, i - 1])
        diagonal = np.concatenate((diagonal, diagonal[-2::-1]))
    else:
        for i in range(1, int(math.floor(Nb / 2)) + 1):
            diagonal = np.append(diagonal, nJ_contribution[Nb - 2, i - 1])
        diagonal = np.concatenate((diagonal, diagonal[::-1]))

    return diagonal


def total_hamiltonian(Nb, Nc, m, J, omega_b, omega_c, g):
    nJ_diag=get_nJ_diagonal(Nb)
    u_arr = [u(j + 1, Nb, Nc, m, J, omega_b, omega_c, g) for j in range(Nb)]
    b_arr = [b(j, Nb, Nc, m, J, omega_b, omega_c, g)  for j in range(Nb + 1)]
    b_arr=b_arr+nJ_diag*J
    H = np.diag(b_arr) + np.diag(u_arr, k=1) + np.diag(u_arr, k=-1)
    return H

def H_b(Nb,omega_b):
    diagonal=[omega_b/2*(j-Nb/2) for j in range(Nb+1)]
    Hb=np.diag(diagonal)    
    return Hb

def state_psi(t, Nb, Nc, m, J, omega_b, omega_c, g):
    H = total_hamiltonian(Nb, Nc, m, J, omega_b, omega_c, g)
    eigenvalues, eigenvectors = np.linalg.eig(H)
    normalized_eigenvectors = np.array([v / np.linalg.norm(v) for v in eigenvectors.T]).T
    U = normalized_eigenvectors
    D = np.diag(np.exp(-1j * eigenvalues * t))
    A = np.dot(U, np.dot(D, U.conj().T))
    
    x = np.zeros((Nb + 1, 1))
    x[0] = 1
    return np.dot(A, x)


def reduced_density(t, Nb, Nc, m, J, omega_b, omega_c, g):
    rho_b=np.zeros((Nb+1,Nb+1))
    psi=state_psi(t, Nb, Nc, m, J, omega_b, omega_c, g)
    for i in range(Nb+1):
        rho_b[i,i]=abs(psi[i,0])**2
    return rho_b

def E_t(t, Nb, Nc, m, J, omega_b, omega_c, g):
    rho_b=reduced_density(t, Nb, Nc, m, J, omega_b, omega_c, g)
    Hb=H_b(Nb,omega_b)
    return np.dot(Hb,rho_b).trace()

def fluctuations(t, Nb, Nc, m, J, omega_b, omega_c, g):
    Et=E_t(t, Nb, Nc, m, J, omega_b, omega_c, g)
    Hb=H_b(Nb, omega_b)
    rho_b=reduced_density(t, Nb, Nc, m, J, omega_b, omega_c, g)
    Hb_sq=np.dot(Hb,Hb)
    A=np.dot(Hb_sq,rho_b)
    return A.trace()-Et**2

def E_t_Nb2_ex(t, Nb, Nc, m, J, omega_b, omega_c, g):
    u1=u(1,Nb,Nc,m,J,omega_b,omega_c,g)
    u2=u(2,Nb,Nc,m,J,omega_b,omega_c,g)
    A=u1**2+u2**2+J**2
    rho11=1/(A*(A-J**2)**2) * (A * u2**4 +u1**4*(J**2+(A-J**2)*np.cos(np.sqrt(A)*t)**2)+np.sqrt(A)*u1**2*u2**2*((np.sqrt(A)+J)*np.cos((np.sqrt(A)-J)*t)+(np.sqrt(A)-J)*np.cos((np.sqrt(A)+J)*t)))
    rho22=u1**2/A*(1-np.cos(np.sqrt(A)*t)**2)
    rho33=u1**2*u2**2/(A*(A-J**2)**2)*(A+J**2-np.sqrt(A)*((np.sqrt(A)+J)*np.cos((np.sqrt(A)-J)*t)+(np.sqrt(A)-J)*np.cos((np.sqrt(A)+J)*t))+(A-J**2)*np.cos(np.sqrt(A)*t)**2)
    return omega_b/2*(rho33-rho11)

def fluc_t_Nb2_ex(t, Nb, Nc, m, J, omega_b, omega_c, g):
    u1=u(1,Nb,Nc,m,J,omega_b,omega_c,g)
    u2=u(2,Nb,Nc,m,J,omega_b,omega_c,g)
    A=u1**2+u2**2+J**2
    rho11=1/(A*(A-J**2)**2) * (A * u2**4 +u1**4*(J**2+(A-J**2)*np.cos(np.sqrt(A)*t)**2)+np.sqrt(A)*u1**2*u2**2*((np.sqrt(A)+J)*np.cos((np.sqrt(A)-J)*t)+(np.sqrt(A)-J)*np.cos((np.sqrt(A)+J)*t)))
    rho22=u1**2/A*(1-np.cos(np.sqrt(A)*t)**2)
    rho33=u1**2*u2**2/(A*(A-J**2)**2)*(A+J**2-np.sqrt(A)*((np.sqrt(A)+J)*np.cos((np.sqrt(A)-J)*t)+(np.sqrt(A)-J)*np.cos((np.sqrt(A)+J)*t))+(A-J**2)*np.cos(np.sqrt(A)*t)**2)
    return omega_b**2/4*(rho11+rho33)-(omega_b/2*(rho33-rho11))**2

t_end=15
t_steps=2000
t=np.linspace(0,t_end,t_steps+1)

E_t_num=np.empty(len(t))
#E_t_ex=np.empty(len(t))
fluc_t_num=np.empty(len(t))
#fluc_t_ex=np.empty(len(t))
jjj=0
for t_step in t:
    E_t_num[jjj]=E_t(t_step, Nb, Nc, m, J, omega_b, omega_c, g)
    fluc_t_num[jjj]=fluctuations(t_step, Nb, Nc, m, J, omega_b, omega_c, g)
    #E_t_ex[jjj]=E_t_Nb2_ex(t_step, Nb, Nc, m, J, omega_b, omega_c, g)
    #fluc_t_ex[jjj]=fluc_t_Nb2_ex(t_step, Nb, Nc, m, J, omega_b, omega_c, g)
    jjj+=1

plt.figure(1)
plt.title('$E(t)$, $N_b={}$, $N_c={}$, $m={}$, $\omega_b={}$, $\omega_c={}$, $J={}$, $g={}$'.format(Nb,Nc,m,omega_b,omega_c,J,g))
plt.plot(t,E_t_num,label='')
#plt.plot(t,E_t_ex,label='ex')
plt.xlabel('$t$')
plt.ylabel('$E(t)$')
#plt.legend()
plt.show()

plt.figure(2)
plt.title('$\Sigma^2(t)$, $N_b={}$, $N_c={}$, $m={}$, $\omega_b={}$, $\omega_c={}$, $J={}$, $g={}$'.format(Nb,Nc,m,omega_b,omega_c,J,g))
plt.plot(t,fluc_t_num)
#plt.plot(t,fluc_t_ex,label='exact')
plt.xlabel('$t$')
plt.ylabel('$\Sigma^2(t)$')
#plt.legend()
plt.show()
