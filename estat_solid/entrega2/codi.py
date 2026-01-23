#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt


#dades de l'enunciat
T = np.array([3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 50, 70, 90, 100, 300, 500, 700, 900, 1100])
lambda_ = np.array([129, 271, 494, 816, 1220, 1675, 2150, 2623, 5446, 4968, 2797, 1618, 1060, 845, 136, 70, 48, 36, 30.5])
c_v = np.array([0.0011, 0.0027, 0.0053, 0.0092, 0.0145, 0.0217, 0.0309, 0.0424, 0.3389, 1.1354, 4.5175, 8.8614, 12.6303, 14.1650, 23.2528, 24.3145, 24.6191, 24.7461, 24.8107])

rho =  2330         # kg/m3
c = 6400            # m/s
kB = 1.38*10**(-23) # J/K
h = 6.63*10**(-34)  # Js
M= 28               #g/mol
N_A=6.022*10**(23)    #avogadro

c_v=c_v*1000*rho/M

#b)---------------------------------------------------------------------------
plt.figure(figsize=(12,4))
#Conductivitat vs T
plt.subplot(1, 2, 1)
plt.grid()
plt.plot(T, lambda_, marker='o', linestyle='-', color= "purple", linewidth=1, markersize=3)
plt.xscale("log")
plt.yscale("log")
plt.xlabel('T (K)')
plt.ylabel(r'$\lambda$ (W/mK)')
plt.title('Conductivitat tèrmica per la mostra de silici')
#plt.grid(True, which="both", ls="--")

#Calor específica molar vs T
plt.subplot(1, 2, 2)
plt.grid()
plt.plot(T, c_v, marker='o', linestyle='-', color= "olive", linewidth=1, markersize=3)
plt.xscale("log")
plt.yscale("log")
plt.xlabel('T (K)')
plt.ylabel(r'$c_v$ (J/K$m^3$)')
plt.title('Calor específica molar per la mostra de silici')
#plt.grid(True, which="both", ls="--")

plt.tight_layout()
plt.show()

#c)----------------------------------------------------------------------------
# Ajust a baixes temperatures: c_v = A T^3
temp_b = T <= 40
T_b = T[temp_b]
cv_b = c_v[temp_b]

log_T = np.log(T_b)
log_cv = np.log(cv_b)

coeffs = np.polyfit(log_T, log_cv, 1)
A = np.exp(coeffs[1])

# Densitat atòmica
n = rho / M * N_A

theta_D = ((12*np.pi**4/5) * n * kB / A)**(1/3)

print(f"Temperatura de Debye estimada ≈ {theta_D:.1f} K")

# Representació de l'ajust T^3
plt.figure()
plt.plot(T,  c_v, marker='o', linestyle='-', color= "purple", linewidth=1, markersize=3)
plt.plot(T, A*T**3, '--', label=r'Ajust $c_v \propto T^3$', color="olive")
plt.xscale('log')
plt.yscale('log')
plt.xlabel('T (K)')
plt.ylabel(r'$c_v$ (J m$^{-3}$ K$^{-1}$)')
plt.legend()
plt.grid()
plt.show()


#d)----------------------------------------------------------------------------


l = 3*lambda_/(c_v*c)
plt.plot(T, l, marker='o', linestyle='-', color= "purple", linewidth=1, markersize=3)
plt.title("Camí lliure mitja en funció de T")
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('T (K)')
plt.ylabel('l (m)')
plt.show()

#e)----------------------------------------------------------------------------

#quan se satura l vol dir que l comparable amb L de l'objecte aprox a T=9
l_max=l[6]
T_max = T[6]
print(f"l_saturació = {l_max:.2e} m a T = {T_max} K")


#f)----------------------------------------------------------------------------
#podem aproximar L=l si dominan boundary effects (tau approx tau_b=L/v_s)
#llavors l=v_s*tau approx L

L= l_max
print(rf"Longitud del sistema L ≈ {L:.2e} m")

#g)----------------------------------------------------------------------------

#efectes fonó-fonó a altes temperatures
temp_a=T>=500
T_a=T[temp_a]
l_a=l[temp_a]
log_l=np.log(l_a)
log_T=np.log(T_a)
coeffs = np.polyfit(log_T, log_l, 1) 
alpha = coeffs[0]
print(f"Exponent a={alpha:.2f}")

#h)----------------------------------------------------------------------------

lb= (3 * lambda_) / (c_v * c)
L_vals = [1e-6, 100e-9, 10e-9] # 1um, 100nm, 10nm
etiquetes = [r"1 $\mu$m", "100 nm", "10 nm"]
colors = ['tab:purple', 'tab:olive', 'tab:cyan']

l_L={}
lambda_L={}


for L in L_vals:
    l_eff=L*lb/(L+lb)
    lamb_L=c_v*c*l_eff/3
    l_L[L]=l_eff
    lambda_L[L]=lamb_L

plt.figure(figsize=(12,4))
plt.subplot(1,2,1)
plt.plot(T, lb , 'ko-', label="original")
for i,L in enumerate(L_vals):
    plt.plot(T,l_L[L],'o-', color=colors[i], label=f"L = {etiquetes[i]}")
    
plt.title("Camí lliure mitja en funció de T per a diferents L")
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('T (K)')
plt.ylabel('l (m)')
plt.legend()

plt.subplot(1,2,2)
plt.plot(T, lambda_ , 'ko-', label="original")
for i,L in enumerate(L_vals):
    plt.plot(T,lambda_L[L],'o-', color=colors[i], label=f"L = {etiquetes[i]}")
    
plt.title("Conductivitat per a diferents L")
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('T (K)')
plt.ylabel(r'$\lambda$ (W/Km)')
plt.legend()
plt.show()

#i)----------------------------------------------------------------------------
#mirem al grafic els valors 
l_L_list = list(l_L.values())
L_list = list(l_L.keys())

l_sat1 = l_L_list[0][13]
l_sat2 = l_L_list[1][10]
l_sat3 = l_L_list[2][10]

print(f"{l_sat1:.2e}m, {T[13]}K per a L = {L_list[0]}m")
print(f"{l_sat2:.2e}m, {T[10]}K per a L = {L_list[1]}m")
print(f"{l_sat3:.2e}m, {T[10]}K per a L = {L_list[2]}m")


