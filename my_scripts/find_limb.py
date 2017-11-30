#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 16:03:51 2017

@author: prabhu
"""

import numpy               as np 
import matplotlib.pyplot   as plt
from scipy.interpolate import UnivariateSpline
import scipy.optimize as so
from   astropy.io            import fits

res = fits.open('/home/prabhu/sunrise_holly/work_IMaX_pole_scripts/imax_find_sun/imax_find_sun_input.fits')

res = res[0].data

icont= res[0,4,:,:]
ind = []
for i in range(icont.shape[0]):
    print(i)
    test = icont[80:850,i]

#fig = plt.figure()
#plt.plot(np.gradient(test),'+')
    spl = UnivariateSpline(np.arange(len(test)), np.gradient(test), k=5)
    spl.set_smoothing_factor(1000)
#plt.plot(spl(np.arange(len(test))), label='Smooth Fct 1e3')
    spl.set_smoothing_factor(10000)
#plt.plot(spl(np.arange(len(test))), label='Smooth Fct 1e4')
#plt.legend(loc='lower left')

#    max_idx = np.argmax(spl(np.arange(len(test))))
#    plt.vlines(max_idx, -300, 300, linewidth=5, alpha=0.3)

    F = lambda x: -spl(x)

    t = np.argmax(spl(np.arange(140)))-5
    ee = so.fmin(F,t)
    ind.append((i,ee[0]))
#plt.vlines(ee[0], -300, 300, linewidth=5, alpha=0.3)