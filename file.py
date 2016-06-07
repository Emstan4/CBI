# -*- coding: utf-8 -*-
"""
Created on Mon May 30 15:33:57 2016

@author: StudyCentre
"""

from __future__ import division
import numpy as np


alpha, beta, gamma, PO =0.1, 0.1,2.2, 1.5
mu_max, theta= 0.25, 0.08
Cp_s = 35
Cp = 14
mu = mu_max*(1 - Cp/Cp_s)
klg = 150 # 1/h
c_o2 = 31/1000 # g
r_mt_max = klg*c_o2*(1/32)*.21
Cx = 12/24.6
#               v0      v1     v2      v3       v4      v5      v6      v7      v8      v9      v10     v11
S1 = np.matrix([[-1,1+alpha,     1,      0,      0,      0,      0,      0,      0,      0,      0,      0],
                
                [0,      0,     -1,    1.5,   0.75,      0,      0,      0,      0,      0,      0,      0],

                [0,      0,      0,      0,     -1,      1,      1,      0,      0,      0,      0,      0],

                [0,      0,      0,      1,      0,      0,      1,     -1,      0,      0,      0,      0],

                [0,      0,      0,      0,      0,      1,      0,      0,      1,     -1,      0,      0],

                [0,    0.1,      0,    0.5,  -0.25,      0,      0,    1/3,      0,      0,      0,     -1],

                [0,      0,      0,      2,      0,      0,     -1,      0,      0,      0,      0,      0],

                [0,      0,      0,      0,      0,      0,      0,    2/3,     -1,      0,      0,      0],

                [0, -gamma,      0,    0.5,   0.25,      0,      0,      0,   0.25,      0,      3,      0], 
             
                [0,   beta,    1/3,    0.5,      0,   -0.5,      0,      0,    0.5,      0,     -2,      0],

                [0,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],

                [0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0],])
                
                   
C = np.matrix( [0,      0,      0,      0,      0,      0,      0,      0,      theta,  0,     mu,      r_mt_max/Cx]).T
r = np.linalg.solve(S1, C)        #solving for rates instead of yields
Cx = 8/24.6


ro2 = r_mt_max/Cx


print r[-3]

