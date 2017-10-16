import os
import glob

import numpy               as np 
import matplotlib.pyplot   as plt

from   astropy.io            import fits
from   matplotlib.ticker     import FuncFormatter

###################### BEGIN ###################### 

# Get data

files = glob.glob('../Data/normalised_output_*.fits')

for el in files:

    res = fits.open(el)
    res = res[0].data
    
    split_name = el.split('_')
    save_out = '../Data/mean_rem_output_' + split_name[2]

    print save_out    

    # Find ave intensity of each continuum and remove from images
    
    for stok in range(1, 4):
        mean = np.mean(res[stok, 4, 240:640, 100:700])
        print mean
        for wvln in range (0, 5):
            res[stok, wvln, :, :] -= mean
    
    # Save corrected data
    
    hdu = fits.PrimaryHDU(res)
    
    hdu.writeto(save_out)
