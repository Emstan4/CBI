# -*- coding: utf-8 -*-
"""
Created on Mon May 30 15:33:57 2016

@author: StudyCentre
"""

from __future__ import division
from numpy.linalg import solve
import numpy as np
from matplotlib import pyplot as plot

a = 0.1
b = 0.1
g = 1.8
mu_max = 0.25
theta_max = 0.32
csa_m = 35/30
Qf = Q = 0
Cxf = Csf = Claf = Cff = Cef = 0
Cxo = 0.02/24.6
Cso = 100/30
V = 2000
Co = np.array([3.33,  0.000813  ,    0,      1.0])     #[X, S, P, V]             
No = Co*V
Km = 0.006
Kp = 8/30
Csa_k = 35/30
def r_prime(C, mu_max, theta_max):
    Cs, Cx, Csa = C
    mu = mu_max*(1 - (Csa/Csa_k))*(Cs/(Km+Cs))
    theta = (theta_max*(1 + (Csa/Kp))**-1)*(Cs/(Km+Cs))
    if mu < 0:
        mu = 0
    #               vo      v1       v2       v3      v4      v5      v6      v7      v8      v9      v10
    #               vo       v1       v2     v3      v4       v5      v6      v7      v8    v9   v10
    S = np.matrix([[-1,    1+a,       1,      0,      0,      0,      0,      0,      0,    0,    0],
                   [0,       0,      -1,    1.5,   0.75,      0,      0,      0,      0,    0,    0],
                   [0,       0,       0,      0,     -1,      1,      1,      0,      0,    0,    0],
                   [0,       0,       0,      1,      0,      0,      1,     -1,      0,    0,    0],
                   [0,       0,       0,      0,      0,      1,      0,      0,      1,   -1,    0],
                   [0,       0,       0,      2,      0,      0,     -1,      0,      0,    0,    0],
                   [0,       0,       0,      0,      0,      0,      0,    2/3,     -1,    0,    0],
                   [0,       b,     1/3,    0.5,      0,   -0.5,      0,      0,    0.5,    0,   -2],
                   [0,      -g,     1/3,      0,    1/4,      0,      0,      0,    1/4,    0,    3],
                   [0,       0,       0,      0,      0,      0,      0,      0,      0,    0,    1],
                   [0,       1,       0,      0,      0,      0,      0,      0,      0,    0,    0]])
                   
    C = np.matrix( [0,       0,       0,      0,      0,      0,      0,      0,  theta,    0,   mu]).T
    Y = solve(S, C)
    return [-Y[0,0], Y[1,0],  Y[9,0]]

def dNdt_fun(N,t):
    Cs, Cx, Cla, V = N[0]/N[3],N[1]/N[3],N[2]/N[3], N[3]  #calculating and naming concentration and volume 
    r=r_prime([Cs, Cx, Cla], mu_max, theta_max)
    
    return [
            Qf*Csf-Q*Cs+(r[0])*Cx*V, 
            Qf*Cxf-Q*Cx+(r[1])*Cx*V,              
            Qf*Claf-Q*Cla+(r[2])*Cx*V, 
            Qf-Q]
           
def solver(tmax):
    concentrations = []
    rates = []
    yields = []
    z = []
    dt = 0.1          
    tspan=np.arange(0,tmax,dt)            
    No = Co*V
    for t in tspan:        
        dndt = dNdt_fun(No,t)
        N = No + np.array(dndt)*dt
        No = N
        Cs, Cx, Cla = N[0]/N[3],N[1]/N[3],  N[2]/N[3]
        r = r_prime([Cs, Cx, Cla], mu_max, theta_max) # overall rates
        r_g = r_prime([Cs, Cx, Cla], mu_max, 0) # growth rates
        r_m = r_prime([Cs, Cx, Cla], 0, theta_max) #maintenance rates
        z.append(r_g[0]/r[0])
        y_x, y_lac = r[1]/abs(r[0]), r[2]/abs(r[0])
        yields.append([y_x, y_lac])
        rates.append(r)
        concentrations.append([Cs, Cx, Cla]) 
        
    concentrations = np.array(concentrations)
    yields = np.array(yields)
    plot.plot(tspan, concentrations)
    plot.show()

solver(300)
