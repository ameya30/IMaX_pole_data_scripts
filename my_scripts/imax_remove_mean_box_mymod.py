import os
import glob

import numpy               as np 
import matplotlib.pyplot   as plt

from   astropy.io            import fits
from   matplotlib.ticker     import FuncFormatter

###################### BEGIN ###################### 

# Get data

files = glob.glob('/home/prabhu/sunrise_holly/normalize_mu_output/*.fits')

for el in files:

    res = fits.open(el)
    res = res[0].data
    
    split_name = el.split('_')
    if len(split_name)==6:
        save_out = '/home/prabhu/sunrise_holly/mean_rem_output/mean_rem_output_' + split_name[-1]
    else:
        save_out = '/home/prabhu/sunrise_holly/mean_rem_output/mean_rem_output_nr_' + split_name[-1]
    print (save_out)

    # Find ave intensity of each continuum and remove from images
    
    for stok in range(1, 4):
        mean = np.mean(res[stok, 4, 240:640, 100:700])
        print (mean)
        for wvln in range (0, 5):
            res[stok, wvln, :, :] -= mean
    
    # Save corrected data
    
    hdu = fits.PrimaryHDU(res)
    
    hdu.writeto(save_out)
