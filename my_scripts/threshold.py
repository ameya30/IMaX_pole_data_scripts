#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 16:32:05 2017

@author: prabhu
"""
#%%

from scipy.ndimage import convolve
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
def rotate(ang):
    uNew1 = qMes*np.sin(2*ang) + uMes*np.cos(2*ang)
    qNew1 = qMes*np.cos(2*ang) + uMes*-1*np.sin(2*ang)
    
    return qNew1,uNew1

files = ['/home/prabhu/sunrise_holly/imax_lp_max_/tr_imax_lp_max_norm_nr_21.fits']
fa = '/home/prabhu/sunrise_holly/imax_lp_max_/tr_imax_roat_angle_Q_U_norm_nr_21.fits'
#for f in files:
for i in range(0, 1):
    f = files[i]
 
    res = fits.open(f)
 
    q = res[0].data
    
    u = fits.open(f)[1].data #1 is minimized U
    a = fits.open(fa)[0].data

    combq_p = np.sum(q[0:-1],axis=0)
    combu_p = np.sum(u[0:-1],axis=0)
    dim = combq_p.shape


    combu = np.empty(shape=(dim[0],dim[1]))
    combq = np.empty(shape=(dim[0],dim[1]))
    angA = np.empty(shape=(dim[0],dim[1]))

    
    for j in range(dim[0]):
        for i in range(dim[1]):
            qMes = combq_p[j,i]
            uMes = combu_p[j,i]
            if qMes<0:
                combq[j,i],combu[j,i] = rotate(np.pi/2)
                angA[j,i]=a[j,i]+np.pi/2
            else:
                combq[j,i] = qMes
                combu[j,i] = uMes
                angA[j,i] = a[j,i]
                
    kernel= np.ones((3,3))
    kernel[1,1] = 0.5
    kernel[0,:] = 0.5/8
    kernel[2,:] = 0.5/8
    kernel[1,0] = 0.5/8
    kernel[1,2] = 0.5/8
    #
    smooth = convolve(combq,kernel,mode='constant')        

    
    # Create mu angle array
    
    dimRes = np.shape(smooth)
    
    mu_arr = np.zeros(shape = (dimRes[0], dimRes[1]))
    
#    mu_arr = np.flipud(mu_arr)
    
    for y in range(0, dimRes[0]):
    
        for x in range(0, dimRes[1]):
            mu_arr[y,x] = find_mu(y,x)

    # Crop arrays so edge effects don't affect normalisation
    
    mu_copy = np.copy(mu_arr)
    
    mu_arr = mu_arr[80 : 936-80, 80 : 936-80]
    smooth = smooth[80 : 936-80, 80 : 936-80]
    
    # Take average of the mu slices 
    
    range_mu = np.max(mu_arr) - np.min(mu_arr)
    
    step = range_mu / no_slices 
    
    mean_list = np.empty(no_slices)
    
    smooth_copy = np.copy(smooth)
    
    for i in range (0, no_slices):
    
        edge1 = i * step
        edge2 = edge1 + step
    
        smooth[mu_arr <= edge1] = 0
        smooth[mu_arr >  edge2] = 0

        rr = smooth[smooth != 0]
        m = np.mean(rr)
        s = np.std(rr)
        nr = rr[rr>(m+3*s)]
        while nr.shape[0]:
            m = np.mean(rr)
            s = np.std(rr)
            nr = rr[rr>(m+3*s)]
            rr = rr[rr<(m+3*s)]
        mean_list[i] = s  
        smooth = np.copy(smooth_copy)
    
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
    
    xax = np.linspace(cenAn[0], cenAn[no_slices-1], 100)
    plt.plot(xax, fit_fn(xax))
    
    thres = fit_fn(mu_copy)