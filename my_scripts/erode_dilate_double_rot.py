#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 11:31:16 2017

@author: prabhu
"""

import numpy as np
import cmocean
from astropy.io import fits
from matplotlib import pyplot as plt
from scipy.io import readsav
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
            

norm_q = n[1,:,:,:]
ncombq = np.sum(norm_q[0:-1],axis=0)

dim = q.shape
norm_v = n[3,:,:,:]
norm_ic =n[0,4,:,:]
combv = norm_v[0,:,:] + norm_v[1,:,:] - norm_v[2,:,:] - norm_v[3,:,:]

kernel= np.ones((3,3))
kernel[1,1] = 0.5
kernel[0,:] = 0.5/8
kernel[2,:] = 0.5/8
kernel[1,0] = 0.5/8
kernel[1,2] = 0.5/8
#
smooth = convolve(combq,kernel,mode='constant') 
std_smooth = np.std(smooth[310:350,430:480])

#smooth = convolve(ncombq,kernel,mode='constant') 
#std_smooth = np.std(smooth[310:350,430:480])
#std_smooth = np.std(smooth[660:700,560:585])

#with inflection from find_limb and that as input to my mod of imax_find_sun
theta_fit = np.linspace(-np.pi/2+np.pi/52.025,-np.pi/2+np.pi/26.5, 1800)
xCen = -1044.9200117954174
yCen = 16862.994425535013
radSun = 16767.565720491271
###########
dimRes = np.shape(smooth)
no_slices = 12   
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
    s = np.std(smooth)
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
del smooth
smooth = convolve(combq,kernel,mode='constant') 
smooth[smooth<(3*thres)]=90000 
############
x_fit = xCen + (radSun) * np.cos(theta_fit)
y_fit = yCen + (radSun) * np.sin(theta_fit)
limbc = np.array(list(set(list(zip(np.around(y_fit[0:-6]),np.around(x_fit[0:-6])))))).astype(int)

indc = np.meshgrid(*map(np.arange,combq.shape),indexing='ij')

for i in range(limbc.shape[0]):
    indc[0][0:limbc[i,0]+1,limbc[i,1]]=0
limbmask = np.ma.masked_equal(indc[0],0).mask


#plotting circle fitting the limb on combined LP map or continuum intensity
fig = plt.figure(figsize=(12,12))
ax = plt.axes()
im=plt.imshow(n[0,4,:,:],vmin = 0,vmax= 1.5,cmap = 'gray')
fig.colorbar(im)
fig.tight_layout(pad=1.8)
plt.plot(x_fit, y_fit, 'r-')
plt.gca().invert_yaxis()
plt.suptitle('Continuum Intensity')

#masking values greater than 3 sigma since all are positive after rotation
x = np.ma.greater(smooth,sig*std_smooth)
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

le1 = erode(y,structure=np.ones((3,3))).astype(x.dtype)
le = dilate(le1,structure=[[True,True,True],[True,True,True],[True,True,True]]).astype(le1.dtype)
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


levels = 1
fig2 = plt.figure(figsize=(12,12))
ax2 = plt.axes()
im2 = plt.imshow(combq,vmax=0.01,vmin=-0.01)
im2_1 = plt.contour(le_wolimb,levels,origin='lower',colors = 'r')
fig2.colorbar(im2)
fig2.tight_layout(pad=1.8)
plt.plot(x_fit, y_fit, 'b')
plt.gca().invert_yaxis()
plt.suptitle('Combined Q')

levels = 1
fig2 = plt.figure(figsize=(12,12))
ax2 = plt.axes()
im2 = plt.imshow(combv,vmax = 0.03,vmin=-0.03)
im2_1 = plt.contour(le_wolimb,levels,origin='lower',colors = 'r')
fig2.colorbar(im2)
fig2.tight_layout(pad=1.8)
plt.plot(x_fit, y_fit, 'b')
plt.gca().invert_yaxis()
plt.suptitle('Combined V')

size,psize,nsize,avglp,avglpp,avglpn,avgcv,avgic,ang,tlp,tcv,stda= [[] for i in range(12)]

#labels go from 1 to 967 not zero, 

for i in range(1,island_lab[1]+1):
    mk = np.ma.masked_equal(labels,i).mask
    size.append(labels[labels==i].shape[0])
    avglp.append(np.mean(combq[mk]))
    avgcv.append(np.mean(combv[mk]))
    avgic.append(np.mean(norm_ic[mk]))
    ang.append(np.median(angA[mk]))
    stda.append(np.std(angA[mk]))
    tlp.append(np.sum(combq[mk]))
    tcv.append(np.sum(combv[mk]))
#
for i in range(len(size)):
    if size[i]>700 or size[i]<5:
        size[i] = np.nan
        avglp[i] = np.nan
        avgcv[i] = np.nan
        avgic[i] = np.nan
        ang[i] = np.nan
        stda[i] = np.nan
        tlp[i] = np.nan
        tcv[i] = np.nan
#%%

'''PLots that are important namely :
- size vs angle with mean
- size vs angle with hist
- size vs average linear polarisation 
- size vs total linear polarisation
'''

stdb = np.array(stda)*50
plt.figure()        
plt.scatter(size,np.array(ang)*180/np.pi)#,s=stdb)
plt.xlabel('Size')
plt.ylabel('Angle')
plt.suptitle('Size vs Angle (angle mean)')

plt.figure()        
plt.scatter(size,stda)#,s=stdb)
plt.xlabel('Size')
plt.ylabel('Std')
plt.suptitle('Size vs std of Angle')


angh = []
a_deg = angA*180/np.pi
for i in range(1,island_lab[1]+1):
    mk = np.ma.masked_equal(labels,i).mask
    no,bin_edge = np.histogram(a_deg[mk],bins=10)
    bincen = 0.5*(bin_edge[1:]+bin_edge[:-1])
    angh.append(np.mean(bincen[no==np.max(no)]))

   
plt.figure()        
plt.scatter(size,angh,s=stdb)
plt.xlabel('Size')
plt.ylabel('Angle')
plt.suptitle('Size vs Angle (angle with hist)')

#fitting trend to absolute value of LP,
#already attempted a trend line separately
#+ve and -ve LP separately, the trend seems to be
#mirrored, hence will proceed to fit line to abs LP value
size_won = list(filter(lambda x:not(np.isnan(x)),size)) #removing nans, polyfit doesnt like them much
avglp_won = list(filter(lambda x:not(np.isnan(x)),avglp))
tlp_won = list(filter(lambda x:not(np.isnan(x)),tlp))
polyp1 = np.polyfit(size_won,np.abs(avglp_won),1)
polyp2 = np.polyfit(size_won,np.abs(tlp_won),1)
xaxis = np.unique(np.sort(size_won))


plt.figure()
plt.scatter(size_won,np.abs(avglp_won),marker='.')
plt.ylabel('Avg Linear Pol')
plt.xlabel('Size in pixels')
plt.plot(xaxis,np.polyval(polyp1,xaxis),'r')
plt.ylim([0.0029,0.008])
plt.suptitle('size vs ALP')

plt.figure()
plt.scatter(size_won,np.abs(tlp_won),marker='.')
plt.ylabel('Total Linear Pol')
plt.xlabel('Size in pixels')
plt.plot(xaxis,np.polyval(polyp2,xaxis),'r')
plt.suptitle('size vs TLP')
#%%
'''
Additional plots
- Size vs ACV
- ALP vs ACV
- Map of angle(with double rotation)
- Tlp vs Tcv
- Size vs Tcv

'''

plt.figure()
plt.scatter(size,avgcv,marker='.')
plt.ylabel('Avg Combined V')
plt.xlabel('Size in pixels')
plt.suptitle('Size vs ACV')


plt.figure()
plt.scatter(avglp,avgcv,marker='.')
plt.ylabel('Avg combined V')
plt.xlabel('Average LP')
plt.suptitle('ALP vs ACV')
#


fig3 = plt.figure(figsize=(12,12))
ax3 = plt.axes()
im3 = plt.imshow(angA*180/np.pi,cmap=cmocean.cm.phase)
im3_1 = plt.contour(le_wolimb,levels,origin='lower',colors = 'k')
fig3.colorbar(im3)
fig3.tight_layout(pad=1.8)
plt.plot(x_fit, y_fit, 'b')
plt.gca().invert_yaxis()
plt.suptitle('Angle')



plt.figure()
plt.scatter(tlp,tcv,marker='.')
plt.ylabel('total combined V')
plt.xlabel('total LP')
plt.suptitle('Tlp vs Tcv')


plt.figure()
plt.scatter(size,tcv,marker='.')
plt.ylabel('total combined V')
plt.xlabel('size')
plt.suptitle('Size vs TCV')


