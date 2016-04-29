# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 10:13:19 2016

@author: StudyCentre
"""

from __future__ import division
import numpy as np
from matplotlib import pyplot as plot

    # X     S     P       ATP
Yg = [0.25, 1 ,  0.52,    0.23]
Ym = [0,    1 ,   0.8,    0.33]
#            Xo So Po Qo Q 
init_cond = [0, 0, 0, 0, 0]
Vo = 5000       
#                 X          S      P     Vo/Vo             #1000 liter initial volume
Co = np.array([0.0006,    4.667,    0,      1.0])             #[X, Gluc, Gly, Et and Vo/Vo] at t=0 in cmol/L
No = Co*Vo
mumax, thethamax, Km = 0.22 , 0.12, 0.0004
Ks = 0.5
#Kp = 8/30
#Csa_k = 35/30
def r_prime(C):
    
    Cx, Cs, Csa, V = C
    C = np.zeros(len(Yg))
    
    for i in range(len(C)):
        C[i] = N[i]/N[-1]
    C[-1] = N[-1] 
    
    mu = mumax*((1 + (C[1]/Ks))**-2)*(C[1]/(Km+C[1]))
    thetha = thethamax*(C[1]/(Km+C[1]))
    if mu < 0:                 # there may be more ethanol production via maintenance, don't want negative mu's
        mu=0
    r = [mu]
    for i in range(len(Yg)-2):
        r.append((Yg[i+1]/Yg[0])*mu + (Ym[i+1]/Ym[-1])*thetha)        
    r[1] = -r[1]
    return r
rplot = []
def dNdt_fun(N,t):
    C = np.zeros(len(Yg))
    for i in range(len(C)):
        C[i] = N[i]/N[-1]
    C[-1] = N[-1] 
    r=r_prime(C)
    rplot.append(r[0])
    var = []
    for i in range(len(r)):        
        var.append(init_cond[-2]*init_cond[i] - init_cond[-1]*C[i] + r[i]*C[0]*C[-1])
    var.append(init_cond[-1]-init_cond[-2])
    return var
tspan = np.arange(0,1443.7,0.1)
dt = tspan[1]

concentration = []
Cx_list =  []
Cs_list =  []
Csa_list =  []
r_list = []
for i, t in enumerate(tspan):
     mat = dNdt_fun(No, t)
     mat = np.array(mat)
     N = No + dt*mat
     No = N
     C = np.zeros(len(Yg))
     for i in range(len(C)):
         C[i] = N[i]/N[-1]
     C[-1] = N[-1]
     r=r_prime(C)
     concentration.append(C)
#     Cx_list.append(C[0])
#     Cs_list.append(C[1])
#     Csa_list.append(C[2])

concentration = np.array(concentration)
#plot.plot(tspan, Csa_list ,color='red',label='$C_{SA}$')
#plot.plot(tspan, Cs_list,color='black',label='$C_{S}$')
#plot.plot(tspan, Cx_list,color='blue',label='$C_{X}$')
plot.plot(tspan, concentration) 
#plot.plot(tspan, rplot)    
#plot.legend(loc='best')

plot.ylabel('Concentration cmol/L') 
plot.xlabel('time (h)') 
plot.ylim([-5,Co[1]+10])
plot.show()
     
#ans = r_prime([0,0,0])
#print ans    
