# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 14:56:45 2016

@author: user
"""

from __future__ import division
from numpy.linalg import solve
import numpy as np
from matplotlib import pyplot as plot
a = 0.1
b = 0.1
g = 1.8
#cla = np.arange(0,2.5,.01)
#npoints = len(cla)
#
#mumax = 0.55
#Kp = 4/30
#mu = mumax*(1 + (cla/Kp))**-1
#thetha = 0.2
#               v1       v2     v3      v4       v5      v6      v7      v8    v9   v10
S = np.matrix([[1+a,     1,      0,      0,      0,      0,      0,      0,    0,    0],
               [0,      -1,    1.5,   0.75,      0,      0,      0,      0,    0,    0],
               [0,       0,      0,     -1,      1,      1,      0,      0,    0,    0],
               [0,       0,      1,      0,      0,      1,     -1,      0,    0,    0],
               [0,       0,      0,      0,      1,      0,      0,      1,   -1,    0],
               [0,       0,      2,      0,      0,     -1,      0,      0,    0,    0],
               [0,       0,      0,      0,      0,      0,    2/3,     -1,    0,    0],
               [b,     1/3,    0.5,      0,   -0.5,      0,      0,    0.5,    0,   -2],
               [-g,    1/3,      0,    1/4,      0,      0,      0,    1/4,    0,    3],
               [0,       0,      0,      0,      0,      0,      0,      0,    0,    1]])
#              
C = np.matrix( [1,       0,      0,      0,      0,      0,      0,      0,    0,    0]).T
#S[8] = [1,0,0,0,0,0,0,0,0,0]
Yg = solve(S, C)
co2 = 0.25*Yg[3] - 2*0.25*Yg[7] - 0.5*Yg[2]

yatp = (1/6)*Yg[7] -0.25*Yg[3]+(1/3)*Yg[1]

print yatp
print
print co2
print
print Yg

#Ysatp_g=-1/3*Yg[1]+2/3*Yg[3]


#Ym = solve(S, C)
#Ysatp_m=-1/3*Ym[1]+2/3*Ym[3]
##rates
#
#r_la = (Yg[4]/Yg[0])*mu + (Ym[4]/Ysatp_m)*thetha
#
#r_et = (Yg[-1]/Yg[0])*mu + (Ym[-1]/Ysatp_m)*thetha
#r_s = (1/Yg[0])*mu + (1/Ysatp_m)*thetha
#
##observed yields
#Yla_o = np.asarray(r_la/r_s)
#Yla_o = np.reshape(Yla_o, npoints)
#Yet_o = r_et/r_s
#
#plot.plot(cla, Yla_o)
#plot.ylim([0.8, 0.92])
#plot.show()


