# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 10:13:19 2016

@author: StudyCentre
"""

from __future__ import division
import numpy as np
from matplotlib import pyplot as plot


a = .1
b = .1
g = 2.5
mu = 0.12
theta = 0.123
lambda_max = 0.24245455/3

Sg = np.matrix([[1+a,1,0.75,0],
                [b,2,0,-2],
                [-g,2/3,0,3],
                [1,0,0,0]])   
Cg = np.matrix([[1,0,theta/lambda_max,0]]).T                
Yg = np.linalg.solve(Sg,Cg)

Yxs_g = 1/Yg[0]

#print Yg#[2]*(33.5/30)

#Sg[2] = [1,0, 0, 0]       
Yg = np.linalg.solve(Sg,Cg)
Yatp_sm = (2/3)*Yg[1] + 3*Yg[-1]

#lambda_max = Yxs_g*mu + (1/Yatp_sm)*theta
#print lambda_max
#print Yxs_g, Yatp_sm
#print Yg  
dt = 0.1
tspan = np.arange(0,50,dt)
Cx = []
Cs = []

tex = [0,48]
cx_ex = [0.03/24.9, 9.515/24.9]
cs_ex = [100/30, 76.63/30]
cxi = cx_ex[0]
csi = cs_ex[0]

print Yg[-1][0,0]*lambda_max*cx_ex[1]*1000/3

#maximum possible glycolyctic flux

r_max = (Yxs_g[0,0]*mu + (1/Yatp_sm[0,0])*theta)*(33/24.6)
#print r_max
