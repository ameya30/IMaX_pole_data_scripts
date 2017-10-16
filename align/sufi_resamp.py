import os
import glob

import matplotlib.pyplot   as plt
import matplotlib.gridspec as gridspec

import numpy               as np 

from   astropy.io import fits
from   matplotlib.ticker import FuncFormatter
from   skimage import transform

sufi_300_ps = 0.020690
sufi_397_ps = 0.019832
imax_525_ps = 0.054458

ratio = imax_525_ps / sufi_397_ps 

files = glob.glob('../Data/im_397_20090612_161239_01000ms_0106786_PD_91m_reclevel2_fltrcls1.fits.gz')
for f in files:
    
    n = f.split('_')
    
    save_name = '../Data/sufi_resamp_' + n[1] + '_' + n[3] +'.fits'
    print save_name
    
    orig = fits.open(f)
    orig = orig[0].data
    [xdim, ydim] = np.shape(orig) 
    
    maxOrig = orig.max()
    
    new = transform.resize(orig/maxOrig, (int(xdim/ratio), int(ydim/ratio)), mode = 'constant') 
    #new = maxOrig * new
    
    hdu = fits.PrimaryHDU(new)
    hdu.writeto(save_name)
