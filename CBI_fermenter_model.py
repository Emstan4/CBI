# -*- coding: utf-8 -*-
"""
Created on Wed Apr 06 19:37:44 2016

@author: StudyCentre
"""

from __future__ import division
import numpy as np
from matplotlib import pyplot as plot
from scipy.integrate import odeint





#   gluc      X         SA      Atp
Yg=[1, 0.15, 0.65 , 0.25] 
Ym=[1, 0, 0.85, 0.33] 

Qf = Q = Cxf = Csf  =  Csaf = 0                    #For now it is a batch fermenter
Vo = 4000       
#                 X          Glu    SA      Vo/Vo             #1000 liter initial volume
Co = np.array([0.0007,      3,    0,      1.0])             #[X, Gluc, Gly, Et and Vo/Vo] at t=0 in cmol/L
No = Co*Vo                                  #total initial cmol in the fermenter and volume at the end (see below why)
mumax,thethamax, Km, Ks= 0.4, 0.05, 0.0005, 2



def r_prime(C):
    
    Cx, Cs, Csa = C
    
    mu = mumax*((1 + (Cs/Ks))**-1)*(Cs/(Km+Cs))
    thetha = thethamax*(Cs/(Km+Cs))
    r_x = mu
    r_s = -((Yg[0]/Yg[1])*mu + (Ym[0]/Ym[3])*thetha)    
    r_sa = (Yg[2]/Yg[1])*mu + (Ym[2]/Ym[3])*thetha
       
    return [r_x, r_s, r_sa]
            
            
def dNdt_fun(N,t):
    Cx, Cs, Csa, V = N[0]/N[3],    N[1]/N[3],  N[2]/N[3],  N[3]  #calculating and naming concentration and volume 
    r=r_prime([Cx, Cs, Csa])
    
    return [Qf*Cxf - Q*Cx + (r[0])*Cx*V,
            Qf*Csf - Q*Cs + (r[1])*Cx*V,    
            Qf*Csaf - Q*Csa + (r[2])*Cx*V,
            Qf - Q]  
            
tspan = np.arange(0,40,0.01)
dt = tspan[1]
Cx_list =  []
Cs_list =  []
Csa_list =  []
r_list = []
#for i, t in enumerate(tspan):
#     mat = dNdt_fun(No, t)
#     mat = np.array(mat)
#     N = No + dt*mat
#     No = N
#     
#     Cx, Cs, Csa, V = N[0]/N[3],    N[1]/N[3],  N[2]/N[3],  N[3]
#     
#     r = r_prime([Cx, Cs, Csa])
#     rsa = Cx*r[2]*30.03
#     r_list.append(rsa)
#     Cx_list.append(Cx)
#     Cs_list.append(Cs)
#     Csa_list.append(Csa)

N = odeint(dNdt_fun, No, tspan)
 
Cx=N[:,0]/N[:,3]                           #devide cmol amount by the volume to get concentration 
Cs=N[:,1]/N[:,3]
Cp=N[:,2]/N[:,3]
    
plot.plot(tspan, Cp ,color='red',label='$C_{SA}$')
plot.plot(tspan, Cs,color='black',label='$C_{S}$')
plot.plot(tspan, Cx,color='blue',label='$C_{X}$')
     
#plot.legend(loc='best')
#print np.average(r_list)
#plot.plot(tspan, r_list)
plot.ylabel('Concentration cmol/L') 
plot.xlabel('time (h)') 
plot.show()