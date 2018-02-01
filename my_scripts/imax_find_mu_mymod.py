#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 08:50:55 2017

@author: prabhu
"""

import os
import glob

import numpy               as np 
import matplotlib.pyplot   as plt
from scipy.io import readsav
from   astropy.io            import fits
from   matplotlib.ticker     import FuncFormatter

#################### FUNCTIONS #################### 

# Format axis to arc seconds

def arcsec(x, pos):

    pixScale = 0.054458 
    
    return '%2.1f' % (x * pixScale)

# Find mu of each pixel

def find_mu(y,x):
    
    distToCen = np.sqrt((x - xCen) ** 2 + (y - yCen) ** 2)
    
    if distToCen <= radSun:
    
        theta = np.arcsin(distToCen / radSun)

        return np.cos(theta)

    else:

        return 0

################### GLOBAL VARS ###################

# Plot figure? 1 = yes 0 = no

plot_fig = 1

# Best estimate data in pixels: from imax_find_sun.py and then from playing around with numbers until it stopped being off limb and dividing by zero

radSun = np.around(16769.565720491271)

yCen = np.around(16862.994425535013)
xCen = np.around(-1044.9200117954174)

# number of mu slices

no_slices = 12

###################### BEGIN ###################### 

# Get data

#files = glob.glob('../Data/post_demod_*.fits')
files = ['/home/prabhu/sunrise_holly/movie_data/mk_magneto_tr2_reduc_rnr_300_21.sav']

#for f in files:
for i in range(0, 1):
    f = files[i]
    split_f = f.split('_')

#    save_name = '../Data/normalised_output_test_' + split_f[3][0:2] + '.fits'
    save_name = '/home/prabhu/sunrise_holly/normalize_mu_output/normalize_mu_21.fits'

    print (save_name)

#    res = fits.open(f)
    res = readsav(f,python_dict = True)['iidn']
#    res = res[0].data
    
    I_cont = np.copy(res[0, 4, :, :])
    
#    I_cont = np.flipud(I_cont)
    
    # Create mu angle array
    
    dimRes = np.shape(I_cont)
    
    mu_arr = np.zeros(shape = (dimRes[0], dimRes[1]))
    
#    mu_arr = np.flipud(mu_arr)
    
    for y in range(0, dimRes[0]):
    
        for x in range(0, dimRes[1]):
            mu_arr[y,x] = find_mu(y,x)

    # Crop arrays so edge effects don't affect normalisation
    
    mu_copy = np.copy(mu_arr)
    
    mu_arr = mu_arr[80 : 936-80, 80 : 936-80]
    I_cont = I_cont[80 : 936-80, 80 : 936-80]
    
    # Take average of the mu slices 
    
    range_mu = np.max(mu_arr) - np.min(mu_arr)
    
    step = range_mu / no_slices 
    
    mean_list = np.empty(no_slices)
    
    I_copy = np.copy(I_cont)
    
    for i in range (0, no_slices):
    
        edge1 = i * step
        edge2 = edge1 + step
    
        I_cont[mu_arr <= edge1] = 0
        I_cont[mu_arr >  edge2] = 0
    
        mean_list[i] = np.mean(I_cont[I_cont != 0])
           
        I_cont = np.copy(I_copy)
    
    # Find center of slices
    
    cenAn = np.empty(no_slices)
    
    for i in range(0, no_slices):
       edge1 = i * step
    
       if i == no_slices-1:
           edge2 = np.max(mu_copy) + 0.000001
       else:
           edge2 = edge1 + step
    
       cenAn[i] = np.mean([edge1, edge2])
    
    # Find curve of normalistion
    
    polyFit = np.polyfit(cenAn, mean_list, 6)
    
    fit_fn  = np.poly1d(polyFit)
#    
#    # Normalise the data
#    
    res_copy = np.copy(res)
    
    norm_arr = fit_fn(mu_copy)
    
    norm_arr[mu_copy == 0] = np.min(norm_arr)
    xax = np.linspace(cenAn[0], cenAn[no_slices-1], 100)
#    
    for stok in range (0, 4):
    
        for wvln in range(0, 5):
#            res[stok, wvln, :, :] = np.flipud(res[stok, wvln, :, :])
           
            res[stok, wvln, :, :] = res[stok, wvln, :, :] / norm_arr
    
    if plot_fig == 1:
    
        # Set up figure
        
        fig, axes = plt.subplots(nrows = 1,
                                 ncols = 3,
                                 figsize = (12.5, 3.3))
        
        fig.subplots_adjust(left   = 0.04,
                            right  = 0.99, 
                            top    = 0.985,
                            bottom = 0.15,
                            wspace = 0.13,
                            hspace = 0.05)
        
        # Plot image
        
        old = axes[0].imshow(np.flipud(res_copy[0, 4, 80:936-80, 80:936-80]), cmap = 'gray', clim = [0, 5000])
        nor = axes[1].imshow(norm_arr[80:936-80, 80:936-80],       cmap = 'gray')
        new = axes[2].imshow(res[0, 4, 80:936-80, 80:936-80],      cmap = 'gray', clim = [0, 1.5])
        
        # Sort labels and ticks
        
        fig.colorbar(old, ax = axes[0], cmap = 'gray', shrink = 0.99,  pad = 0.01)
        fig.colorbar(nor, ax = axes[1], cmap = 'gray', shrink = 0.99,  pad = 0.01)
        fig.colorbar(new, ax = axes[2], cmap = 'gray', shrink = 0.99,  pad = 0.01)
        
        formatter = FuncFormatter(arcsec)
        axes[0].set_ylabel('arcseconds') 
        for ax in axes:
        
            ax.xaxis.set_major_formatter(formatter)
            ax.yaxis.set_major_formatter(formatter)
            ax.set_xlabel('arcseconds') 
            ax.invert_yaxis()
        
        plt.savefig('/home/prabhu/sunrise_holly/normalize_mu_output/imax_normalise_mu_22.png', dpi = 500)
        
        # Plot fitting curve
        
        fig, axes = plt.subplots(figsize = (6, 4))
        
        xax = np.linspace(cenAn[0], cenAn[no_slices-1], 100)
        
        axes.plot(cenAn, mean_list, '.')
        axes.plot(xax, fit_fn(xax))
        
        axes.set_xlabel('Mu Angle')
        axes.set_ylabel('Mean Intensity')
        
        plt.savefig('/home/prabhu/sunrise_holly/normalize_mu_output/imax_normalisation_curve_22.png', dpi = 500, bbox_inches = 'tight')
    
    # Save corrected data
    
#    hdu = fits.PrimaryHDU(res)
#    
#    hdu.writeto(save_name)
    
    