# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 14:22:09 2016

@author: user
"""

from __future__ import division
import numpy as np
from matplotlib import pyplot as plot

mu_max = 0.55
Kp = 4

conc = np.linspace(0,100,100)
mu = mu_max*(1 + conc/Kp)**-1

plot.plot(conc, mu)
plot.plot()
