#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 09:21:04 2025

@author: davidrosamolina
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# parameters 
Omega_L = 30
Omega_p = 1
gamma_21 = 2
gamma_31 = 2
Delta_L = 0
Delta_p= np.linspace(-20, 20, 500)


y13= []
x13= []


#we assign index in order rhos,x,y with 1,2,3
#(rho11,rho22,rho33,x12,y12,x13,y13,x23,y23)

for p in Delta_p:
    A = np.zeros((9, 9))
    b = np.zeros(9)

    #eq1
    A[0, 0] = 1.0
    A[0, 1] = 1.0
    A[0, 2] = 1.0
    b[0] = 1.0

    #eq2
    A[1, 1] = -gamma_21
    A[1, 4] = -Omega_L

    #eq3
    A[2, 2] = -gamma_31
    A[2, 6] = -Omega_p

    #eq4
    A[3, 3] = -gamma_21/2
    A[3, 4] = -Delta_L
    A[3, 8] = Omega_p/2

    #eq5
    A[4, 0] = -Omega_L/2
    A[4, 1] = Omega_L/2
    A[4, 3] = Delta_L
    A[4, 7] = Omega_p/2
    A[4, 4] = -gamma_21/2

    #eq6
    A[5, 5] = -gamma_31/2
    A[5, 6] = -p
    A[5, 8] = -Omega_L/2

    #eq7
    A[6, 0] = -Omega_p/2
    A[6, 2] = Omega_p/2
    A[6, 5] = p
    A[6, 7] = Omega_L/2
    A[6, 6] = -gamma_31/2

    #eq8
    A[7, 7] = -(gamma_21+gamma_31)/2
    A[7, 8] = Delta_L-p
    A[7, 6] = -Omega_L/2
    A[7, 4] = -Omega_p/2

    #eq9
    A[8, 5] = Omega_L/2
    A[8, 3] = -Omega_p/2
    A[8, 7] = p-Delta_L
    A[8, 8] = -(gamma_21+gamma_31)/2

    
    try:
        sol = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        sol, *_ = np.linalg.lstsq(A, b, rcond=None)

    #save x_13 and y_13
    x13.append(sol[5])
    y13.append(sol[6])
    
# Plot
plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.plot(Delta_p, y13, color='red', label=r'$Y_{13}$ Absorption')
plt.xlabel(r'$\Delta_p$')
plt.ylabel(r'$Y_{13}$')
plt.title('Autler–Townes Doublet - Absorption')
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(Delta_p, x13, color='green', label=r'$X_{13}$ Dispersion')
plt.xlabel(r'$\Delta_p$')
plt.ylabel(r'$X_{13}$')
plt.hlines(0,-20,20, color='gray',linestyle='--', alpha=0.5)
plt.vlines(0,-0.04,0.04, color='gray',linestyle='--', alpha=0.5)
plt.title('Autler–Townes Doublet - Dispersion')
plt.legend()

plt.tight_layout()
plt.show()