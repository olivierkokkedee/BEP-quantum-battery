# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 22:44:31 2023

@author: olivi
"""

import pandas as pd
from itertools import combinations
import math

N_max = 24
#generate all possible states with a given number of battery spins N and spin ups m
def generate_spin_states(N, m):
    states = []
    for indices in combinations(range(N), m):
        state = [1 if i in indices else 0 for i in range(N)]
        states.append(state)
    return states

#nJ is the number of contributions n J. If two consecutive spins are opposite oriented,
#nJ will be incremented with 1.
def calculate_nJ(N, m):
    nJ = 0
    states = generate_spin_states(N, m)
    for lst in states:
        nJ_temp = 0
        for i in range(N-1):
            if lst[i] != lst[i+1]:
                nJ_temp += 1
        nJ += nJ_temp
    return nJ

#the remaining code converts the data in an excel sheet
results = pd.DataFrame(columns=['N'] + [f'm_{m}' for m in range(1, N_max//2 + 1)])


for N in range(2, N_max + 1):
    m_values = list(range(1, math.floor(N / 2) + 1))
    nJ_values = [calculate_nJ(N, m) for m in m_values]
    row_data = [N] + nJ_values + [None] * (N_max//2 - len(nJ_values))  # Fill remaining columns with None
    results = results.append(pd.Series(row_data, index=results.columns), ignore_index=True)


results.to_excel('nJ_contribution_per_state.xlsx', index=False)

