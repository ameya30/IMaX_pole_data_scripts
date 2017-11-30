#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 14:55:38 2017

@author: prabhu
"""

from scipy.io import readsav
import numpy as np
#from matplotlib import pyplot as plt
from astropy.io import fits


unbin = readsav("/home/prabhu/sunrise_holly/movie_data/mk_magneto_tr2_reduc_rnr_300_21.sav")

bind = fits.open("/home/prabhu/sunrise_holly/binned_cycles/binned_tr2_mk_restor_300_21_2.fits")

#for bi=2
y1,y2,x1,x2 = 120,170,125,400
#for bi=3
#y1,y2,x1,x2 = 125,145,130,190
#for bi=4
#y1,y2,x1,x2 = 80,105,110,150

#not restored
qcont_unbin_nr = unbin['iidn'][1,4,:,:]
icont_unbin_nr = unbin['iidn'][0,4,:,:]
qcont_bin_nr = bind[1].data[1,4,:,:]
icont_bin_nr = bind[1].data[0,4,:,:]

#fig = plt.figure(figsize=(12,12))
#ax = plt.axes()
#im = plt.imshow(qcont_unbin_nr,cmap='gray')
#plt.show()
#plt.suptitle('qcont')
#fig1 = plt.figure(figsize=(12,12))
#ax1 = plt.axes()
#im1 = plt.imshow(icont_unbin_nr,cmap='gray')
#plt.show()
#plt.suptitle('icont')


nqcont_unbin_nr = qcont_unbin_nr/icont_unbin_nr
nqcont_bin_nr = qcont_bin_nr/icont_bin_nr

std_unbin_nr = np.std(nqcont_unbin_nr[286:336,136:382])
std_bin_nr = np.std(nqcont_bin_nr[y1:y2,x1:x2])

#restored
qcont_unbin_r = unbin['iid'][1,4,:,:]
icont_unbin_r = unbin['iid'][0,4,:,:]
qcont_bin_r = bind[0].data[1,4,:,:]
icont_bin_r = bind[0].data[0,4,:,:]

nqcont_unbin_r = qcont_unbin_r/icont_unbin_r
nqcont_bin_r = qcont_bin_r/icont_bin_r

std_unbin_r = np.std(nqcont_unbin_r[286:336,136:382])
std_bin_r = np.std(nqcont_bin_r[y1:y2,x1:x2])



print("the ratio of nonrestored unbin/bin")
print(std_unbin_nr/std_bin_nr)
print("the ratio of restored unbin/bin")
print(std_unbin_r/std_bin_r)