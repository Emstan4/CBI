# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 14:56:45 2016

@author: user
"""

from __future__ import division
from numpy.linalg import solve
import numpy as np

a = 0.1
b = 0.1
g = 1.8

#               v1       v2     v3      v4       v5      v6      v7      v8    v9   v10
S = np.matrix([[1+a,     1,      0,      0,      0,      0,      0,      0,    0,    0],
               [0,      -1,    1.5,   0.75,      0,      0,      0,      0,    0,    0],
               [0,       0,      0,     -1,      1,      1,      0,      0,    0,    0],
               [0,       0,      1,      0,      0,      1,     -1,      0,    0,    0],
               [0,       0,      0,      0,      1,      0,      0,      1,   -1,    0],
               [0,       0,      2,      0,      0,     -1,      0,      0,    0,    0],
               [0,       0,      0,      0,      0,      0,    2/3,     -1,    0,    0],
               [b,     1/3,    0.5,      0,   -0.5,      0,      0,    0.5,    0,   -2],
               [-g,      0,    0.5,    1/4,      0,      0,      0,    1/4,    0,    3],
               [0,       0,      0,      0,      0,      0,      0,      0,    0,    1]])
#              
C = np.matrix( [1,       0,      0,      0,      0,      0,      0,      0,    0,    0]).T
S[8] = [1,0,0,0,0,0,0,0,0,0]
Yg = solve(S, C)
co2 = 0.25*Yg[3] - 2*0.25*Yg[7] - 0.5*Yg[2] - a*Yg[0]

yatp = 0.25*Yg[3]+(1/2)*Yg[2] + 0.25*Yg[7]

print yatp
print
print co2
print
print Yg



