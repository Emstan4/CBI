# -*- coding: utf-8 -*-
"""
Created on Wed Apr 06 19:37:44 2016

@author: StudyCentre
"""

from __future__ import division
import numpy as np
from matplotlib import pyplot as plot





#   gluc      X         Gly     Lac       Form        Eth        Atp
Yg = [1,    0.14992504,     0,   0.77136432,   0.02023988,     0.04047976,     0.27736132] # [0 Gluc, 1 X, 2 Gly, 3 Et, 4 CO2, 5 ATP]
Ym = [1,         0,     0,         1,        0,          0,     1/3] 

Qf = Q = Cxf = Csf = Cgf = Cef = Clacf = Cff = 0                    #For now it is a batch fermenter
Vo = 500       
#              X     Glu   Gly  Lac Form etha Vo/Vo             #1000 liter initial volume
Co = np.array([0.00159, 3.33,  0,   0,   0,   0,   1.0])             #[X, Gluc, Gly, Et and Vo/Vo] at t=0 in cmol/L
No = Co*Vo                                  #total initial cmol in the fermenter and volume at the end (see below why)
mumax, thethamax, Km = 0.55 , 0.2, 0.008
Kp = 4/30

def r_prime(C):
    
    Cx, Cs, Cg, Clac, Cf, Ce = C
    
    mu = mumax*((1 + (Clac/Kp))**(-1))*(Cs/(Km+Cs))
    thetha = thethamax#*(Cs/(Km+Cs))
    r_x = mu
    r_s = -((Yg[0]/Yg[1])*mu + (Ym[0]/Ym[6])*thetha)
    r_g = (Yg[2]/Yg[1])*mu + (Ym[2]/Ym[6])*thetha
    r_l = (Yg[3]/Yg[1])*mu + (Ym[3]/Ym[6])*thetha
    r_f = (Yg[4]/Yg[1])*mu + (Ym[4]/Ym[6])*thetha
    r_e = (Yg[5]/Yg[1])*mu + (Ym[5]/Ym[6])*thetha
    
    return [r_x, r_s, r_g, r_l, r_f, r_e]
             
def dNdt_fun(N,t):
    Cx, Cs, Cg, Clac, Cf, Ce, V = N[0]/N[6],    N[1]/N[6],  N[2]/N[6],  N[3]/N[6],  N[4]/N[6],  N[5]/N[6],  N[6]  #calculating and naming concentration and volume 
    r=r_prime([Cx, Cs, Cg,Clac, Cf, Ce])
    
    return [Qf*Cxf - Q*Cx + (r[0])*Cx*V,
            Qf*Csf - Q*Cs + (r[1])*Cx*V,    
            Qf*Cgf - Q*Cg + (r[2])*Cx*V,
            Qf*Clacf - Q*Clac + (r[3])*Cx*V,
            Qf*Cff - Q*Cf + (r[4])*Cx*V,            
            Qf*Cef - Q*Ce + (r[5])*Cx*V,
            Qf - Q]  
            
tspan = np.arange(0,30,0.001)
dt = tspan[1]
Cx_list =  []
Cs_list =  []
Cg_list =  []

Cl_list =  []
Cf_list =  []

Ce_list =  []

acc_yld_list = []
ins_yld_list = []
rate_lac = []
rate_glu = []
for i, t in enumerate(tspan):
     mat = dNdt_fun(No, t)
     mat = np.array(mat)
     N = No + dt*mat
     No = N
     Cx, Cs, Cg, Clac, Cf, Ce, V = N[0]/N[6],    N[1]/N[6],  N[2]/N[6],  N[3]/N[6],  N[4]/N[6],  N[5]/N[6],  N[6]
     r = r_prime([Cx, Cs, Cg,Clac, Cf, Ce])
     acc_yld = Clac/(Co[1] - Cs)
     acc_yld_list.append(acc_yld)
     rate_lac.append(r[3])
     rate_glu.append(r[1])
     ins_yld_list.append(r[3]/abs(r[1]))
     
     #print N[3]
     Cx_list.append(N[0]/N[6])
     Cs_list.append(N[1]/N[6])
     Cg_list.append(N[2]/N[6])
     Cl_list.append(N[3]/N[6])
     Cf_list.append(N[4]/N[6])
     Ce_list.append(N[5]/N[6])
     


plot.plot(Cl_list, ins_yld_list, color='red',label='instant')
plot.plot(Cl_list, acc_yld_list, color='black',label='accumu')

plot.ylim([0.8,0.92])
plot.legend(loc='best')
plot.ylabel('Concentration cmol/L') 
plot.xlabel('time (h)') 
plot.show()