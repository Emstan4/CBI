# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 19:58:58 2016

@author: StudyCentre
"""

from __future__ import division
import numpy
import pandas

Sg = numpy.matrix(pandas.read_excel('constant.xlsx'))
def constants(M, mode):
    if mode == 'main':
        a = numpy.zeros(len(M))
        a[0] = 1
        M[-1] = a
        a = numpy.matrix(a).T
    elif mode == 'grow':
        a = numpy.zeros(len(M))
        a[0] = 1
        a = numpy.matrix(a).T
    return a
        
C = constants(Sg, 'grow')
Yg = numpy.linalg.solve(Sg,C)        

print Yg