# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 10:13:19 2016

@author: StudyCentre
"""

from __future__ import division
import numpy as np
from matplotlib import pyplot as plot
from scipy.integrate import odeint
from scipy.optimize import fsolve


Csex=np.array([ 4.        ,  3.99809878,  3.99452932,  3.98782783,  3.97524613,
        3.95162487,  3.90727832,  3.82402459,  3.66773752,  3.37438079,
        2.82385119])

tex=np.array([  0. ,   4.2,   8.4,  12.6,  16.8,  21. ,  25.2,  29.4,  33.6,
        37.8,  42. ])

mumax,thethamax, Km,km2, Kp= 0.15, 0.1, 0.0005,0.00005, 1000
#Cxo, Cso = 0.0003, 3.5
Yg=[0.18, 1, 0.5 , 0.3] 
Ym=[0, 1, 0.65, 0.4]                     
Vo=10000                                             
Cxf=Csf=Cpf=0 

Cso=2.82385119
Cxo = 0.0003
Co=np.array([Cxo, Cso, 0, 1])     #[X, S, P, V]             
No=Co*Vo
Qf=Q=50                            # define throughflow 
Csf=4                             # substrate in feed 
Cxf=Cpf=0
init_cond = [0, 4, 0, 0, 0]
def r_prime(C):

    mu = mumax*C[1]/(Km+C[1])*((1 + C[2]/Kp)**(-1))
    thetha = thethamax*(C[1]/(km2+C[1]))
    r = [mu]
    for i in range(len(Yg)-2):
        r.append((Yg[i+1]/Yg[0])*mu + (Ym[i+1]/Ym[-1])*thetha)        
    r[1] = -r[1]
    return r
    
def steady_state(C,Q):               
    r=r_prime(C)
    vari = []
    Qf = Q
    for i in range(len(r)):  
        vari.append(Qf*init_cond[i] - Q*C[i] + r[i]*C[0]*Vo)
    return vari
Qspan = np.arange(10,1500,1)
D_list = Qspan/Vo
Csteady = []
Co = [0.4,0,2.2]
for qi in Qspan:           
    Co=fsolve(steady_state,Co,args=qi)  
    Csteady.append(Co)#Get estimate from graph above
#print C_steady
Csteady = np.array(Csteady)
plot.plot(D_list, Csteady[:,1]) 
#Qspan = np.arange(200,2500,100)
#Csteady = []
#Co = [0.4,0,2.2]
##for qi in Qspan:           
#Co=fsolve(steady_state,Co,args=Q)  
#Csteady.append(Co)
#function = steady_state(Co,Q)
#print function 
#plot.plot(tex,Csex,'o')
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
        var.append(Qf*init_cond[i] - Q*C[i] + r[i]*C[0]*C[-1])
    var.append(Qf - Q)
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
#print "Production rate(0.3548 cmol/L) = ", prod, "cmol/h L"
#print "Volumetric production rate:", prod, "cmol/L.h"
#print "Y accumulated(1430) = ", Y_acc_1430  
#
#plot.plot(tspan, Cp ,color='red',label='$C_{P}$')
#plot.plot(tspan, Cs,color='black',label='$C_{S}$')
#plot.plot(tspan, Cx,color='blue',label='$C_{X}$')
#plot.legend(loc='best')
#plot.ylabel('Concentration cmol/L') 
#plot.xlabel('time (h)') 
#plot.show()
#plot.plot(Cp, rprod)

  
