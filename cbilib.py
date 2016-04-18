# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 19:58:58 2016

@author: StudyCentre
"""

from __future__ import division
import numpy
import pandas

Sg = numpy.matrix(pandas.read_excel('constants.xlsx'))
def constants(M, mode):
    a = numpy.zeros(len(M))
    a[0] = 1
    a[-2] = 0.02
    a = numpy.matrix(a).T
    if mode == 'main':
        M[-1] = a        

    return a
        
C = constants(Sg, 'grow')
Yg = numpy.linalg.solve(Sg,C)        

print Yg