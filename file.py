# -*- coding: utf-8 -*-
"""
Created on Mon May 30 15:33:57 2016

@author: StudyCentre
"""

from __future__ import division
from numpy.linalg import solve
import numpy as np
from matplotlib import pyplot as plot

a = 0.12
b = 0.135
g = 1.85
cla = 0
mu_max = 0.55
theta_max = 0.2
Kp = 4/30

Qf = Q = 0
Cxf = Csf = Claf = Cff = Cef = 0
Cxo = 0.04/24.6
Cso = 100/30
V = 500
Co=np.array([Cso, Cxo, 0, 0, 0, 1])     #[X, S, P, V]             
No = Co*V
Km = 0.008
def r_prime(C, mu_max, theta_max):
    Cs, Cx, Cla, Cf, Ce = C
    mu = mu_max*(Cs/(Km+Cs))*(1 + (Cla/Kp))**-1
    theta = theta_max*(Cs/(Km+Cs))
    
    S = np.matrix([[-1,    1+a,       1,      0,      0,      0,      0,      0,      0],
                   [0,       0,      -1,      1,      1,      0,      0,      0,      0],
                   [0,       0,       0,      0,     -1,      1,      1,      0,      0],
                   [0,       0,       0,      0,      0,      0,     -1,      1,      1],
                   [0,       0,       0,      0,      0,      0,   -1/3,      1,      0],
                   [0,       b,       0,   -1/3,    1/3,   -1/3,      0,      0,     -1],
                   [0,      -g,    -1/3,      0,    2/3,      0,      0,      0,      0],
                   [0,       0,       0,      1,      0,      0,      0,      0,      0],              
                   [0,       1,       0,      0,      0,      0,      0,      0,      0]])
                   
    C = np.matrix( [0,       0,       0,      0,      0,      0,  theta,      0,     mu]).T
    Y = solve(S, C)
    return [-Y[0,0], Y[1,0], Y[5,0], Y[7,0], Y[8,0]]

def dNdt_fun(N,t):
    Cs, Cx, Cla, Cf, Ce, V = N[0]/N[5],N[1]/N[5],N[2]/N[5],N[3]/N[5],N[4]/N[5],N[5]  #calculating and naming concentration and volume 
    r=r_prime([Cs, Cx, Cla, Cf, Ce], mu_max, theta_max)
    
    return [
            Qf*Csf-Q*Cs+(r[0])*Cx*V, 
            Qf*Cxf-Q*Cx+(r[1])*Cx*V,              
            Qf*Claf-Q*Cla+(r[2])*Cx*V, 
            Qf*Cff-Q*Cf+(r[3])*Cx*V,                        
            Qf*Cef-Q*Ce+(r[4])*Cx*V,
            Qf-Q]
           
def solver(tmax):
    concentrations = []
    rates = []
    yields = []
    z = []
    dt = 0.01          
    tspan=np.arange(0,tmax,dt)            
    No = Co*V
    for t in tspan:        
        dndt = dNdt_fun(No,t)
        N = No + np.array(dndt)*dt
        No = N
        Cs, Cx, Cla, Cf, Ce = N[0]/N[5],N[1]/N[5],  N[2]/N[5],  N[3]/N[5],  N[4]/N[5]
        r = r_prime([Cs, Cx, Cla, Cf, Ce], mu_max, theta_max) # overall rates
        r_g = r_prime([Cs, Cx, Cla, Cf, Ce], mu_max, 0) # growth rates
        r_m = r_prime([Cs, Cx, Cla, Cf, Ce], 0, theta_max) #maintenance rates
        z.append(r_g[0]/r[0])
        y_x, y_lac, y_f, y_e = r[1]/abs(r[0]), r[2]/abs(r[0]),r[3]/abs(r[0]),r[4]/abs(r[0])
        yields.append([y_x, y_lac, y_f, y_e])
        rates.append(r)
        concentrations.append([Cs, Cx, Cla, Cf, Ce]) 
        
    concentrations = np.array(concentrations)
    yields = np.array(yields)
    plot.plot(tspan, z)
    plot.show()

solver(35)
