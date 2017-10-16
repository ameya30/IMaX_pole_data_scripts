import os
import glob

import matplotlib.pyplot   as plt
import matplotlib.gridspec as gridspec

import numpy               as np 

from   astropy.io import fits
from   matplotlib.ticker import FuncFormatter

# Global varibales to change 

sigma = 5 

# Funtion to format axis to arc seconds

def arcsec(x, pos):

    pixScale = 50. / 936
    
    return '%2.1f' % (x * pixScale)

# Get data

files = glob.glob('../Data/mean_rem_output_*')
for f in files:
    split_f = f.split('_')
    
    save_name = '../Data/imax_combined_V_output_' + split_f[3]
    
    print save_name
    
    res = fits.open(f)
    res = res[0].data
    
    combV = (res[3, 0, :, :] + res[3, 1, :, :]) - (res[3, 2, :, :] + res[3, 3, :, :])

    try:

        hdu = fits.PrimaryHDU(combV)
        
        hdu.writeto(save_name)

    except:
        pass
   # # Set up figure
   # 
   # fig, axes = plt.subplots(ncols   = 1,
   #                          nrows   = 1, 
   #                          figsize = (6, 5))
   # 
   # 
   # fig.subplots_adjust(left   = 0.04,
   #                     right  = 1.00, 
   #                     top    = 0.99,
   #                     bottom = 0.09,
   #                     wspace = 0.00,
   #                     hspace = 0.00)
   # 
   # # Plot data
   # 
   # im = axes.imshow(combV, cmap = 'gray', clim = [-0.1, 0.1])
   # 
   # fig.colorbar(im, ax = axes, shrink = 0.99, cmap = 'gray', pad = 0.01)
   # 
   # formatter = FuncFormatter(arcsec)
   # 
   # axes.xaxis.set_major_formatter(formatter)
   # axes.yaxis.set_major_formatter(formatter)
   # 
   # axes.set_xlabel('arcsec')
   # axes.set_ylabel('arcsec')
   # 
   # axes.invert_yaxis()
   # 
   # plt.savefig(save_name, dpi = 500)
   # plt.close(fig)
