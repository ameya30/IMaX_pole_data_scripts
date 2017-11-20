#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 18:26:39 2017

@author: prabhu
"""

#from con_re_scipy import congrid
from scipy.io import readsav
from matplotlib import pyplot as plt
import numpy as np
from astropy.io import fits

def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)


data = readsav("/scratch/prabhu/backup_workstation/sunrise_holly/movie_data/mk_magneto_tr2_reduc_rnr_300_22.sav",python_dict=True)

iid = data['iid'] #demodulated and restored 
iidn = data['iidn'] #demodulated and not restored

bi = 2
stok,wvln,dimy,dimx = iid.shape

binr = np.zeros(shape = (stok,wvln,int(dimx/bi),int(dimy/bi)))
binnr = np.zeros(shape = (stok,wvln,int(dimx/bi),int(dimy/bi)))

#binning bi x bi for restored data
for i in range(stok):
	for j in range(wvln):
		ima = iid[i,j,:,:]
		binr[i,j,:,:] = rebin(ima,(int(dimy/bi),int(dimx/bi)))

#binning bi x bi for non-restored data
for i in range(stok):
	for j in range(wvln):
		iman = iidn[i,j,:,:]
		binnr[i,j,:,:] = rebin(iman,(int(dimy/bi),int(dimx/bi)))



# plt.imshow(binr[0,4,:,:],cmap='gray')
# plt.gca().invert_yaxis()
# plt.show()

#creating primary HDU for binned restored data (primary HDU is mandatory in writing fits file)
hdu1 = fits.PrimaryHDU(data=binr)

#creating Image HDU for binned non-restored data
hdu2 = fits.ImageHDU(data=binnr,name="not restored")

#making hdulist to write into a fits file
hdulist = fits.HDUList([hdu1,hdu2])

hdulist.writeto('/scratch/prabhu/backup_workstation/sunrise_holly/binned_cycles/binned_tr2_mk_restor_300_22_' + str(bi) +'.fits')







