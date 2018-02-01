#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 10:06:41 2018

@author: prabhu
"""

import numpy as np

from astropy.io import fits
from matplotlib import pyplot as plt
import matplotlib as mpl
from scipy.ndimage.morphology import binary_dilation as dilate
from scipy.ndimage.morphology import binary_erosion as erode
from scipy.ndimage import measurements 

from scipy.ndimage import convolve


def find_mu(y,x):
    
    distToCen = np.sqrt((x - xCen) ** 2 + (y - yCen) ** 2)
    
    if distToCen <= radSun:
    
        theta = np.arcsin(distToCen / radSun)

        return np.cos(theta)

    else:

        return 0

def rotate(ang):
    uNew1 = qMes*np.sin(2*ang) + uMes*np.cos(2*ang)
    qNew1 = qMes*np.cos(2*ang) + uMes*-1*np.sin(2*ang)
    
    return qNew1,uNew1

sig = 3

fac = 3

if fac ==3:
    kernel = np.ones((3,3))
    kernel[1,1] = 0.5
    kernel[0,:] = 0.5/8
    kernel[2,:] = 0.5/8
    kernel[1,0] = 0.5/8
    kernel[1,2] = 0.5/8    

elif fac == 5:
    kernel = np.ones((5,5))
    kernel[2,2] = 0.5
    kernel[1,1:4] = 0.25/8
    kernel[3,1:4] = 0.25/8
    kernel[2,1] = 0.25/8
    kernel[2,3] = 0.25/8  
    kernel[0,:] = 0.25/16
    kernel[4,:] = 0.25/16
    kernel[1:4,0] = 0.25/16
    kernel[1:4,4] = 0.25/16


cy = input("Choose cycle: ")


re = input("Choose r for restored or nr for nonrestored: ")


if re=='r':
    normal = "/home/prabhu/sunrise_holly/mean_rem_output/mean_rem_output_" + cy + '.fits'
    filenameQ = "/home/prabhu/sunrise_holly/imax_lp_max_/tr_imax_lp_max_norm_" + cy+ '.fits'
    filenamea = "/home/prabhu/sunrise_holly/imax_lp_max_/tr_imax_roat_angle_Q_U_norm_" + cy + '.fits'
	

    n = fits.open(normal)[0].data
    q = fits.open(filenameQ)[0].data #0 is max linear pol map, of maxed Q
    u = fits.open(filenameQ)[1].data #1 is minimized U
    a = fits.open(filenamea)[0].data


elif re=='nr':
    normal = "/home/prabhu/sunrise_holly/mean_rem_output/mean_rem_output_nr_" + cy + '.fits'
    filenameQ = "/home/prabhu/sunrise_holly/imax_lp_max_/tr_imax_lp_max_norm_nr_" + cy+ '.fits'
    filenamea = "/home/prabhu/sunrise_holly/imax_lp_max_/tr_imax_roat_angle_Q_U_norm_nr_" + cy + '.fits'
	
	
    n = fits.open(normal)[0].data
    q = fits.open(filenameQ)[0].data #0 is max linear pol map, of maxed Q
    u = fits.open(filenameQ)[1].data #1 is minimized U
    a = fits.open(filenamea)[0].data

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
            
a_copy = np.copy(angA*180/np.pi)
for j in range(a_copy.shape[0]):
    for i in range(a_copy.shape[0]):
        if a_copy[j,i]>180:
            a_copy[j,i]= a_copy[j,i]-180


            
norm_q = n[1,:,:,:]
ncombq = np.sum(norm_q[0:-1],axis=0)

dim = q.shape
norm_v = n[3,:,:,:]
norm_ic =n[0,4,:,:]
combv = norm_v[0,:,:] + norm_v[1,:,:] - norm_v[2,:,:] - norm_v[3,:,:]


#with inflection from find_limb and that as input to my mod of imax_find_sun
theta_fit = np.linspace(-np.pi/2+np.pi/52.025,-np.pi/2+np.pi/26.5, 1800)
xCen = -1044.9200117954174
yCen = 16862.994425535013
radSun = 16767.565720491271

###########
#
smooth = convolve(combq,kernel,mode='constant')
dimRes = np.shape(smooth)

no_slices = 12

mu_arr = np.zeros(shape = (dimRes[0], dimRes[1]))

#    mu_arr = np.flipud(mu_arr)

for y in range(0, dimRes[0]):

    for x in range(0, dimRes[1]):
        mu_arr[y,x] = find_mu(y,x)

# Crop arrays so edge effects don't affect normalisation

mu_copy = np.copy(mu_arr)
c = 120
mu_arr = mu_arr[c : 936-c, c : 936-c]
smooth = smooth[c : 936-c, c : 936-c]

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


thres = fit_fn(mu_copy)
del smooth
smooth = convolve(combq,kernel,mode='constant') 
smooth[smooth<(sig*thres)]=90000 #for combq

x_fit = xCen + (radSun) * np.cos(theta_fit)
y_fit = yCen + (radSun) * np.sin(theta_fit)
limbc = np.array(list(set(list(zip(np.around(y_fit[0:-6]),np.around(x_fit[0:-6])))))).astype(int)

indc = np.meshgrid(*map(np.arange,combq.shape),indexing='ij')
indcc = np.copy(indc)
for i in range(limbc.shape[0]):
    indc[0][0:limbc[i,0]+20,limbc[i,1]]=0
limbmask = np.ma.masked_equal(indc[0],0).mask

for i in range(limbc.shape[0]):
    indcc[0][0:limbc[i,0]-10,limbc[i,1]]=0
limbmask_map = np.ma.masked_equal(indc[0],0).mask
nlimbmask_map =np.logical_not(np.ma.masked_equal(indc[0],0).mask)



#masking values greater than 3 sigma since all are positive after rotation
x = np.ma.greater(smooth,sig*thres)
#removing pixels off limb from array
x_wolimb= np.ma.array(x,mask=limbmask)
#changing the filled value in masked array from True to Nan
np.ma.set_fill_value(x_wolimb,np.nan)
#obtaining this array with nan values in place of True
filledx = x_wolimb.filled()

#masking values lesser than 3 sigma
#y = np.ma.greater(smooth,sig*std_smooth).astype(int)
y = np.logical_not(np.ma.equal(smooth,90000)).astype(int)

##to make it work on le_wolimb, a masked array, had to fill the masked values
# and extract the filled array

le1 = erode(y,structure=np.ones((fac,fac))).astype(x.dtype)

le = dilate(le1,structure=np.ones((fac,fac)).astype(bool)).astype(le1.dtype)
    
le_wolimb = np.ma.array(le,mask=limbmask)
np.ma.set_fill_value(le_wolimb,0.0)

#labeling islands leftover from erode and dilate without off limb pixels
island_lab = measurements.label(le_wolimb.filled())
labels = island_lab[0]


indx,indq = [[] for i in range(2)]
for j in range(dim[1]):
	for i in range(dim[2]):
		dum = q[:,j,i]
		dumdum = filledx[j,i]
		temp = np.dot(dum[0],dum[1:-1])
		np.abs(dum)
		if (temp[temp<0].size) and (not(np.isnan(dumdum))):
			indq.append((j,i))
		if not(np.isnan(dumdum)):
			indx.append((j,i))
print('Number of pixels above 3sigma in LP with mixed sign')
print(len(indq))
# print(indq[10000:10010])
print('Total number of pixels above 3sigma in LP')
print(len(indx))

combq_wo= np.ma.array(combq,mask=limbmask_map)
combv_wo= np.ma.array(combv,mask=limbmask_map)
np.ma.set_fill_value(combq_wo,-0.01)
np.ma.set_fill_value(combv_wo,-0.03)
fcombq = combq_wo.filled()
fcombv = combv_wo.filled()

cdict = {'red': ((0.0,  0.50, 0.50),
                  (1.0,  1.0, 1.0)),

        'green':  ((0.0,  0.50, 0.50),
                   (1.0,  1.0, 1.0)),
          
        'blue':   ((0.0,  0.50, 0.50),
                  (1.0,  1.0, 1.0))}
our_gray_cmap = mpl.colors.LinearSegmentedColormap('our_gray',cdict,256)
mpl.cm.register_cmap(cmap=our_gray_cmap)

levels = 1
fig1 = plt.figure(figsize=(12,12))
ax1 = plt.axes()
im1 = plt.imshow(fcombq[0:895,58:901],cmap='our_gray', vmin=0, vmax=.01)
im1_1 = plt.contour(le_wolimb[0:895,58:901],levels,origin='lower',
                    colors = 'firebrick')
fig1.colorbar(im1)
fig1.tight_layout(pad=1.8)
plt.suptitle('Linear Polarization')

x = np.arange(-100,900,100)*0.05
ticklabx = [str(i) for i in x.astype(int)]
y = np.arange(100,1100,100)*0.05
ticklaby = [str(i) for i in y.astype(int)]

levels = 1
fig2 = plt.figure(figsize=(12,12))
ax2 = plt.axes()
im2 = plt.imshow(fcombv[0:895,58:901],vmax = 0.03,vmin=-0.03,cmap='gray')
im2_1 = plt.contour(le_wolimb[0:895,58:901],levels,origin='lower',colors = 'firebrick')
r = fig2.colorbar(im2)
plt.xlabel('[arcsec]')
plt.ylabel('[arcsec]')
ax2.set_xticklabels(ticklabx)
ax2.set_yticklabels(ticklaby[::-1])
r.set_label(r'$\sum_{\lambda}\ \ \frac{I_V\ (\lambda)}{I_{cont}}$')
#plt.suptitle(r'$\sum_{\lambda}\ \ \frac{I_V\ (\lambda)}{I_{cont}}$')
fig2.tight_layout(pad=1.8)
