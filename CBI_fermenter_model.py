# -*- coding: utf-8 -*-
"""
Created on Wed Apr 06 19:37:44 2016

@author: StudyCentre
"""

from __future__ import division
import numpy as np
from matplotlib import pyplot as plot






#   gluc      X         SA      Atp
Yg = [1,  0.17057569,  0.93816631,    0.30703625]
Ym = [1,           0,         8/7,    0.38095238]

Qf = Q = Cxf = Csf  =  Csaf = 0                    #For now it is a batch fermenter
Vo = 2000       
#                 X          Glu    SA      Vo/Vo             #1000 liter initial volume
Co = np.array([0.000813,    3.33,    0,      1.0])             #[X, Gluc, Gly, Et and Vo/Vo] at t=0 in cmol/L
No = Co*Vo                                  #total initial cmol in the fermenter and volume at the end (see below why)
mumax, thethamax, Km = 0.25 , 0.32, 0.006
Kp = 8/30
Csa_k = 35/30
def r_prime(C):
    
    Cx, Cs, Csa = C
    
    mu = mumax*(1 - (Csa/Csa_k))*(Cs/(Km+Cs))
    thetha = (thethamax*(1 + (Csa/Kp))**-1)*(Cs/(Km+Cs))
    if mu < 0:                 # there may be more ethanol production via maintenance, don't want negative mu's
        mu=0
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
            
tspan = np.arange(0,250,0.01)
dt = tspan[1]
Cx_list =  []
Cs_list =  []
Csa_list =  []
r_list = []
for i, t in enumerate(tspan):
     mat = dNdt_fun(No, t)
     mat = np.array(mat)
     N = No + dt*mat
     No = N
     Cx, Cs, Csa, V = N[0]/N[3],    N[1]/N[3],  N[2]/N[3],  N[3]
     r = r_prime([Cx, Cs, Csa])
     rsa = Cx*r[2]*30.03
     r_list.append(rsa)
     Cx_list.append(Cx)
     Cs_list.append(Cs)
     Csa_list.append(Csa)

     
#plot.plot(tspan, Csa_list ,color='red',label='$r_{sa}$')
#plot.plot(tspan, Cs_list,color='black',label='$r_{S}$')
#plot.plot(tspan, Cx_list,color='blue',label='$r_{X}$')
     
#plot.legend(loc='best')
print np.average(r_list)
plot.plot(tspan, r_list)
plot.ylabel('Concentration cmol/L') 
plot.xlabel('time (h)') 
plot.show()