# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 10:13:19 2016

@author: StudyCentre
"""

from __future__ import division
import numpy as np
from matplotlib import pyplot as plot
from scipy.integrate import odeint


kq = 0.005
    # X     S     P       ATP
Yg = [0.25, 1 ,  0.52,    0.23]
Ym = [0,    1 ,   0.8,    0.333]
#            Xo So Po Qo Q 
init_cond = [0, 12, 0, kq, 0]
Vo = 5000       
#                 X          S      P     Vo/Vo             #1000 liter initial volume
Co = np.array([0.0006,      0.02,    0,      1.0])             #[X, Gluc, Gly, Et and Vo/Vo] at t=0 in cmol/L
No = Co*Vo
mumax, thethamax, Km = 0.22 , 0.12, 0.0004
Ks = 0.5

def r_prime(C):
        
    mu = mumax*C[1]/(Km+C[1])*((1 + C[1]/Ks)**(-2))
    thetha = thethamax*(C[1]/(Km+C[1]))
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
        Qf = (C[-1]*(r[1]*C[0]-kq))/(C[1]-init_cond[1])
        
        var.append(Qf*init_cond[i] - init_cond[-1]*C[i] + r[i]*C[0]*C[-1])
    var.append(Qf - init_cond[-1])
    return var
    
tspan = np.arange(0,500,0.01)
dt = tspan[1]

N = odeint(dNdt_fun, No, tspan)
Ci=(N[:, :3].T/N[:, 3]).T
r = np.asarray([r_prime(C) for C in Ci])


rx = r[:,0]
rs = -r[:,1]
rp = r[:,2]
#
##instantaneous yield
Y_ins = rp/rs
#
##accumulated yield
Y_acc = N[:,2]/(No[1] - N[:,1])
#
#print Y[-1]
#print Y_acc[-1]
z = (Yg[1]/Yg[0])*rx/rs

y_obs = z*Yg[2] + (1-z)*Ym[2]
Cx=N[:,0]/N[:,3]                           #divide cmol amount by the volume to get concentration 
Cs=N[:,1]/N[:,3]
Cp=N[:,2]/N[:,3]
V = N[:,3]
rprod = Cx*rp
time = np.interp(3.548, Cp, tspan)
overall_y = np.interp(3.548, Cp, y_obs)
prod = 3.548/time
#sol = np.interp(0.00, C, tspan[::10])
#x = ((Co[1] - Cs_1430)/Co[1])*100
#Y_acc_1430= np.interp(1430,tspan , Y_acc)
#vol_prod = max(Cp)/t_end
#print "Maximum rate:", overall_y
print "Production rate(0.3548 cmol/L) = ", prod, "cmol/h L"
#print "Volumetric production rate:", prod, "cmol/L.h"
#print "Y accumulated(1430) = ", Y_acc_1430  

plot.plot(tspan, Cp ,color='red',label='$C_{P}$')
plot.plot(tspan, Cs,color='black',label='$C_{S}$')
plot.plot(tspan, Cx,color='blue',label='$C_{X}$')
plot.legend(loc='best')
plot.ylabel('Concentration cmol/L') 
plot.xlabel('time (h)') 
plot.show()
#plot.plot(Cp, rprod)

  
