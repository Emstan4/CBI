# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 14:56:45 2016

@author: user
"""

from __future__ import division
from numpy.linalg import solve
import numpy as np
from matplotlib import pyplot as plot
a = 0.12
b = 0.135
g = 1.85
cla = np.arange(0,2.5,.01)
npoints = len(cla)

mumax = 0.55
Kp = 4/30
mu = mumax*(1 + (cla/Kp))**-1
thetha = 0.2

S = np.matrix([[1+a,     1,      0,      0,      0,      0,      0,      0],
               [0,      -1,      1,      1,      0,      0,      0,      0],
               [0,       0,      0,     -1,      1,      1,      0,      0],
               [0,       0,      0,      0,      0,     -1,      1,      1],
               [0,       0,      0,      0,      0,   -1/3,      1,      0],
               [b,       0,   -1/3,    1/3,   -1/3,      0,      0,     -1],
               [-g,   -1/3,      0,    2/3,      0,      0,      0,      0],
               [0,       0,      1,      0,      0,      0,      0,      0]])
#              
C = np.matrix( [1,       0,      0,      0,      0,      0,      0,      0]).T
Yg = solve(S, C)

Ysatp_g=-1/3*Yg[1]+2/3*Yg[3]

S[6] = [1,0,0,0,0,0,0,0]
Ym = solve(S, C)
Ysatp_m=-1/3*Ym[1]+2/3*Ym[3]
#rates

r_la = (Yg[4]/Yg[0])*mu + (Ym[4]/Ysatp_m)*thetha

r_et = (Yg[-1]/Yg[0])*mu + (Ym[-1]/Ysatp_m)*thetha
r_s = (1/Yg[0])*mu + (1/Ysatp_m)*thetha

#observed yields
Yla_o = np.asarray(r_la/r_s)
Yla_o = np.reshape(Yla_o, npoints)
Yet_o = r_et/r_s

plot.plot(cla, Yla_o)
plot.ylim([0.8, 0.92])
plot.show()


