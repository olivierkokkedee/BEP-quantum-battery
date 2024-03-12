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
os.chdir('C:/Users/private/Bureaublad/BEP')
from scipy.optimize import curve_fit

#Defining all parameters
h=5
n=60
Nb=10
Nc=200
m=100
k=2*n+1
J=1
omega_b=2
omega_c=2
g=1
M=2*m*(Nc-m+1)
A=2*M*g**2+J**2

t_end=2
t_steps=10000
t=np.linspace(0,t_end,t_steps+1)
ftsize=16

#%matplotlib inline
#Opening the excel file which contains all contribution values n for certain states
df = pd.read_excel('C:/Users/olivi/OneDrive/Bureaublad/BEP/nJ_contribution_per_state.xlsx')
nJ_contribution = df.iloc[:, 1:].to_numpy()

#u_j and b_j, as they are defined in Liu et al.
def u(j, Nb, Nc, m, J, omega_b, omega_c, g):
    return g * np.sqrt(j * (Nb - j + 1) * (Nc - m + j) * (m - j + 1))

def b(j, Nb, Nc, m, J, omega_b, omega_c, g):
    return omega_b / 2 * (j - Nb / 2) + omega_c / 2 * (m - j - Nc / 2)

#creating the array which contains all contributions nJ, so I can add it to the Hamiltonian later
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

#the following functions calculate the intermediate steps to calculate the energy and fluctuations
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


#functions that calculate the exact solutions for Nb=2
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

def E_J_eq_g(t,k,h,omega_b):
    return -omega_b/(4*k)*((k+1)*np.cos((k-1)*h*t)+(k-1)*np.cos((k+1)*h*t))

#omega_b/2*1/(A*(A-J**2)**2)*((u2**2-u1**2)*(-Au2+J**2*u1**2+u1**2*(A-J**2)*np.cos(np.sqrt(A)*t))-2*np.sqrt(A)*u1**2*u2**2*((np.sqrt(A)+J)*np.cos((np.sqrt(A)-J)*t)+(np.sqrt(A)-J)*np.cos((np.sqrt(A)+J)*t)))
def fluc_J_eq_g(t,k,h,omega_b):
    return omega_b**2/(32*k**2)*(4*k**2-2*(k**2-1)*np.cos(2*h*t)-(k+1)**2*np.cos(2*(k-1)*h*t)-(k-1)**2*np.cos(2*(k+1)*h*t))

def E_g_gg_J(t,M,g,omega_b):
    return -omega_b/2*np.cos(np.sqrt(2*M)*g*t)

def fluc_g_gg_J(t,M,g,omega_b):
    return omega_b**2/8*np.sin(np.sqrt(2*M)*g*t)**2

def E_global(k,omega):
    return -omega/(4*k)

def fluc_general_Jg_Nb2(t,J,g,omega,n):
    u2=2*(n**2+n)
    A=2*u2*g**2+J**2
    return omega**2/(32*A)*(4*A-2*(A-J**2)*np.cos(2*J*t)-(A**(1/2)+J)**2*np.cos(2*(A**(1/2)-J)*t)-(A**(1/2)-J)**2*np.cos(2*(A**(1/2)+J)*t))

'''
#
#Used to make the 6 plots in the report
#
E_eq=[]
fluc_eq=[]

k_values=[3,5,11]
t_E_eq=[np.empty(4),np.empty(4),np.empty(4)]
t_fluc_eq=[np.empty(4),np.empty(4),np.empty(4)]
for i in range(len(k_values)):
    E_eq.append(E_J_eq_g(t, k_values[i], J, omega_b))
    fluc_eq.append(fluc_J_eq_g(t, k_values[i], J, omega_b))
    t_E_temp=[]
    t_fluc_temp=[]
    t_E_temp.append(np.pi/J*(0+1/k_values[i]))
    t_E_temp.append(np.pi/J*(1-1/k_values[i]))
    t_E_temp.append(np.pi/J*(1+1/k_values[i]))
    t_E_temp.append(np.pi/J*(2-1/k_values[i]))
    t_fluc_temp.append(np.pi/J*(0+1/2*(1-1/k_values[i])))
    t_fluc_temp.append(np.pi/J*(0+1/2*(1+1/k_values[i])))
    t_fluc_temp.append(np.pi/J*(1+1/2*(1-1/k_values[i])))
    t_fluc_temp.append(np.pi/J*(1+1/2*(1+1/k_values[i])))
    t_E_eq[i]=np.array(t_E_temp)
    t_fluc_eq[i]=np.array(t_fluc_temp)

for i in range(len(k_values)):
    plt.figure(1)
    plt.plot(t,E_eq[i],label='$E_b(t)$')
    plt.plot(t_E_eq[i],E_J_eq_g(t_E_eq[i], k_values[i], J, omega_b),'rx',markersize=8,label='$t^E$')
    plt.plot(t_fluc_eq[i],E_J_eq_g(t_fluc_eq[i], k_values[i], J, omega_b),'gx',markersize=8,label='$t^\Sigma$')
    plt.plot(t, np.full(len(t),2*omega_b/4), linestyle='--', label='$E_{max}$')
    plt.xlabel('$t\;\;[1/\omega]$', fontsize=ftsize)
    plt.ylabel('$E_b(t)\;\; [\hbar \omega]$', fontsize=ftsize)
    plt.legend(loc='lower left',bbox_to_anchor=(1/6,0),fontsize=9)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.grid()
    plt.tight_layout()
    plt.savefig('E_J=g=h_k={}.pdf'.format(k_values[i]))
    plt.show()

for i in range(len(k_values)):
    plt.figure(1)
    plt.plot(t,fluc_eq[i],label='$\Sigma^2_b(t)$')
    plt.plot(t_E_eq[i],fluc_J_eq_g(t_E_eq[i], k_values[i], J, omega_b),'rx',markersize=8,label='$t^E$')
    plt.plot(t_fluc_eq[i],fluc_J_eq_g(t_fluc_eq[i], k_values[i], J, omega_b),'gx',markersize=8,label='$t^\Sigma$')
    plt.xlabel('$t\;\;[1/\omega]$', fontsize=ftsize)
    plt.ylabel('$\Sigma^2_b(t)\;\; [\hbar \omega]$', fontsize=ftsize)
    plt.legend(loc='lower left',bbox_to_anchor=(1/6,0),fontsize=10.71)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.grid()
    plt.tight_layout()
    #plt.savefig('fluc_J=g=h_k={}.pdf'.format(k_values[i]))
    plt.show()

'''
#
#Used to make the E_max and fluc_max plots
#
'''
k_values=np.array([2*i+1 for i in range(1,32)])
t_values=np.linspace(3,63,1000)
plt.figure(1)
plt.plot(t_values, omega_b/2*np.cos(np.pi/t_values), color='red', label='Continuous')
plt.plot(k_values, omega_b/2*np.cos(np.pi/k_values), '.',markersize=8, label='Allowed values')
plt.plot(t_values,np.full(len(t_values),omega_b*2/4),linestyle='--',label='$E_{max}$')
plt.xlabel('$k$', fontsize=ftsize)
plt.ylabel('$E_{b,max}\;\;[\hbar \omega]$', fontsize=ftsize)
plt.grid()
plt.legend(fontsize=14)

# Increase font size of tick labels on both axes
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.tight_layout()

plt.savefig('E_max_J=g_Nc=m=2.pdf')
plt.show()

# Plot 2
plt.figure(2)
plt.plot(t_values, omega_b**2/8*(1+np.cos(np.pi/t_values)), color='red',  label='Continuous')
plt.plot(k_values, omega_b**2/8*(1+np.cos(np.pi/k_values)), '.', markersize=8,label='Allowed values')
plt.xlabel('$k$', fontsize=ftsize)
plt.ylabel('$\Sigma^2_{b,max}\;\;[(\hbar \omega)^2]$', fontsize=ftsize)
plt.grid()
plt.legend(fontsize=14)

# Increase font size of tick labels on both axes
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.tight_layout()

plt.savefig('fluc_max_J=g_Nc=m=2.pdf')
plt.show()
'''
'''
#
#Used to check Nb=2 for general 
#
t_J=np.array([np.pi/(2*J)*(2*1-1),np.pi/(2*J)*(2*2-1),np.pi/(2*J)*(2*3-1)])
t_A=np.array([0,np.pi*1/A**(1/2),np.pi*2/A**(1/2),np.pi*3/A**(1/2),np.pi*4/A**(1/2),np.pi*5/A**(1/2),np.pi*6/A**(1/2)])
plt.figure(1)
plt.plot(t, E_t_Nb2_ex(t, Nb, Nc, m, J, omega_b, omega_c, g))
plt.plot(t_J,E_t_Nb2_ex(t_J, Nb, Nc, m, J, omega_b, omega_c, g),'.',markersize=10,label='t_J')
plt.plot(t_A,E_t_Nb2_ex(t_A, Nb, Nc, m, J, omega_b, omega_c, g),'.',markersize=10,label='t_A')
plt.legend()
plt.show()

plt.figure(1)
plt.plot(t, fluc_t_Nb2_ex(t, Nb, Nc, m, J, omega_b, omega_c, g))
plt.plot(t_J,fluc_t_Nb2_ex(t_J, Nb, Nc, m, J, omega_b, omega_c, g),'.',label='t_J')
plt.plot(t_A,fluc_t_Nb2_ex(t_A, Nb, Nc, m, J, omega_b, omega_c, g),'.',label='t_A')
plt.legend()
plt.show()

'''
E_num = np.empty(len(t))
fluc_num=np.empty(len(t))
i=0
for t_step in t:
    E_num[i] = E_t(t_step, Nb, Nc, m, J, omega_b, omega_c, g)
    fluc_num[i]=fluctuations(t_step, Nb, Nc, m, J, omega_b, omega_c, g)
    i+=1

index_E=[]
index_fluc=[]
index_E.append(np.argmax(E_num[0:int(len(t)/5)]))
index_E.append(np.argmax(E_num[int(len(t)/5):int(len(t)/5*3)])+int(len(t)/5))
index_E.append(np.argmax(E_num[int(3*len(t)/5):len(t)])+int(3*len(t)/5))
index_fluc.append(np.argmax(fluc_num[int(0*len(t)/5):int(2*len(t)/5)]))
index_fluc.append(np.argmax(fluc_num[int(2*len(t)/5):int(4*len(t)/5)])+int(2*len(t)/5))
index_fluc.append(np.argmax(fluc_num[int(4*len(t)/5):int(len(t)+1)])+int(4*len(t)/5))

plt.figure(1)
plt.plot(t, E_num,label='$E_b(t)$')
plt.plot([t[index_E[0]],t[index_E[1]],t[index_E[2]]],[E_num[index_E[0]],E_num[index_E[1]],E_num[index_E[2]]],'rx',markersize=8,label='$t^E$')
plt.plot([t[index_fluc[0]],t[index_fluc[1]],t[index_fluc[2]]],[E_num[index_fluc[0]],E_num[index_fluc[1]],E_num[index_fluc[2]]],'gx',markersize=8,label='$t^\Sigma$')
plt.plot(t, np.full(len(t),Nb*omega_b/4), linestyle='--', label='$E_{max}$')
plt.xlabel('$t\;\;[1/\omega]$', fontsize=ftsize)
plt.ylabel('$E_b(t)\;\; [\hbar \omega]$', fontsize=ftsize)
#plt.legend(loc='lower left',bbox_to_anchor=(1/6,0),fontsize=10.71)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.grid()
plt.tight_layout()
plt.savefig('E_Nb={}_m={}_Nc={}_h={}.pdf'.format(Nb,m,Nc,h))
plt.show()

plt.figure(2)
plt.plot(t, fluc_num,label='$\Sigma^2_b(t)$')
plt.xlabel('$t\;\;[1/\omega]$', fontsize=ftsize)
plt.ylabel('$\Sigma^2_b(t)\;\; [(\hbar \omega)^2]$', fontsize=ftsize)
plt.plot([t[index_E[0]],t[index_E[1]],t[index_E[2]]],[fluc_num[index_E[0]],fluc_num[index_E[1]],fluc_num[index_E[2]]],'rx',label='$t^E$')
plt.plot([t[index_fluc[0]],t[index_fluc[1]],t[index_fluc[2]]],[fluc_num[index_fluc[0]],fluc_num[index_fluc[1]],fluc_num[index_fluc[2]]],'gx',markersize=8,label='$t^\Sigma$')
#plt.legend(loc='lower left',bbox_to_anchor=(1/6,0),fontsize=10.71)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.grid()
plt.tight_layout()
plt.savefig('fluc_Nb={}_m={}_Nc={}_h={}.pdf'.format(Nb,m,Nc,h))
plt.show()


'''
E_num = np.empty(len(t))
fluc_num=np.empty(len(t))
i=0
for t_step in t:
    E_num[i] = E_t(t_step, Nb, Nc, m, J, omega_b, omega_c, g)
    fluc_num[i]=fluctuations(t_step, Nb, Nc, m, J, omega_b, omega_c, g)
    i+=1

plt.figure(1)
plt.plot(t, E_num,label='$E_b(t)$')
#plt.plot([t[index_E[0]],t[index_E[1]],t[index_E[2]]],[E_num[index_E[0]],E_num[index_E[1]],E_num[index_E[2]]],'rx',label='$t^E$')
#plt.plot([t[index_fluc[0]],t[index_fluc[1]],t[index_fluc[2]]],[E_num[index_fluc[0]],E_num[index_fluc[1]],E_num[index_fluc[2]]],'gx',label='$t^\Sigma$')
plt.xlabel('$t\;\;[1/\omega]$', fontsize=ftsize)
plt.ylabel('$E_b(t)\;\; [\hbar \omega]$', fontsize=ftsize)
#plt.legend(loc='lower left',bbox_to_anchor=(1/6,0),fontsize=10.71)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.grid()
plt.tight_layout()
plt.savefig('E_Nb={}_m={}_Nc={}_h={}'.format(Nb,m,Nc,h))
plt.show() 
'''
