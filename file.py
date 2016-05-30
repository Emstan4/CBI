# -*- coding: utf-8 -*-
"""
Created on Mon May 30 15:33:57 2016

@author: StudyCentre
"""

from __future__ import division
import numpy
theta = 0.08
mu = 0.25

Sg = numpy.matrix([[-1,1.1,1,0.75,0],
                   [0,0.1,2,0,-2],
                   [0,-2.2,2/3,0,3],
                   [0,1,0,0,0],
                   [0,0,0,1,0]])          
        
Cg = numpy.matrix([[0,0,theta,mu,0]]).T

Yg1 = numpy.linalg.solve(Sg,Cg)
print (Yg1[1,0]/Yg1[0,0])*(24.6/30)