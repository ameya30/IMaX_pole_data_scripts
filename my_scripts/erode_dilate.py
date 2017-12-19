#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 11:01:02 2017

@author: prabhu
"""

import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
from scipy.io import readsav
from scipy.ndimage.morphology import binary_dilation as dilate

from scipy.ndimage.morphology import binary_erosion as erode
from scipy.ndimage import measurements 

from scipy.ndimage import convolve





cy = input("Choose cycle: ")


re = input("Choose r for restored or nr for nonrestored: ")

#for limb
res = fits.open('/home/prabhu/sunrise_holly/work_IMaX_pole_scripts/imax_find_sun/imax_find_sun_input.fits')

res = res[0].data

dum = res[0,4,:,:]
mask = np.ma.masked_inside(dum,0,700).mask #remove was an attempyt at crude limb removal

    

#if re=='r':
#    normal = "/home/prabhu/sunrise_holly/movie_data/mk_magneto_tr2_reduc_rnr_300_" + cy + '.sav'
#    filenameQ = "/home/prabhu/sunrise_holly/imax_lp_max_/imax_lp_max_" + cy+ '.fits'
#    filenamea = "/home/prabhu/sunrise_holly/imax_lp_max_/imax_roat_angle_Q_U_" + cy + '.fits'
#	
#
#    n = readsav(normal,python_dict=True)['iid'] #iid for restored, iidn for nonrestored
#    q = fits.open(filenameQ)[0].data #0 is max linear pol map, of maxed Q
#    u = fits.open(filenameQ)[1].data #1 is minimized U
#    a = fits.open(filenamea)[0].data
#
#
#elif re=='nr':
#    normal = "/home/prabhu/sunrise_holly/movie_data/mk_magneto_tr2_reduc_rnr_300_"+cy+'.sav'
#    filenameQ = "/home/prabhu/sunrise_holly/imax_lp_max_/imax_lp_max_nr_" + cy +'.fits'
#    filenamea = "/home/prabhu/sunrise_holly/imax_lp_max_/imax_roat_angle_Q_U_nr_" + cy +'.fits'
#	
#
#    n = readsav(normal,python_dict=True)['iidn'] #iid for restored, iidn for nonrestored
#    q = fits.open(filenameQ)[0].data #0 is max linear pol map, of maxed Q
#    u = fits.open(filenameQ)[1].data #1 is minimized U
#    a = fits.open(filenamea)[0].data

if re=='r':
    normal = "/home/prabhu/sunrise_holly/normalize_mu_output/normalize_mu_" + cy + '.fits'
    filenameQ = "/home/prabhu/sunrise_holly/imax_lp_max_/imax_lp_max_norm_" + cy+ '.fits'
    filenamea = "/home/prabhu/sunrise_holly/imax_lp_max_/imax_roat_angle_Q_U_norm_" + cy + '.fits'
	

    n = fits.open(normal)[0].data
    q = fits.open(filenameQ)[0].data #0 is max linear pol map, of maxed Q
    u = fits.open(filenameQ)[1].data #1 is minimized U
    a = fits.open(filenamea)[0].data


elif re=='nr':
    normal = "/home/prabhu/sunrise_holly/normalize_mu_output/normalize_mu_nr_" + cy + '.fits'
    filenameQ = "/home/prabhu/sunrise_holly/imax_lp_max_/imax_lp_max_norm_nr_" + cy+ '.fits'
    filenamea = "/home/prabhu/sunrise_holly/imax_lp_max_/imax_roat_angle_Q_U_norm_nr_" + cy + '.fits'
	
	
    n = fits.open(normal)[0].data
    q = fits.open(filenameQ)[0].data #0 is max linear pol map, of maxed Q
    u = fits.open(filenameQ)[1].data #1 is minimized U
    a = fits.open(filenamea)[0].data


#stokv = n[3,:,:,:]
#norm_q = q/n[0,4,:,:]
#norm_v = stokv/n[0,4,:,:]
#norm_ic = n[0,4,:,:]/np.mean(n[0,4,320:400,580:700])
norm_q = q
norm_v = n[3,:,:,:]
norm_ic =n[0,4,:,:]
    
combq = np.sum(norm_q[0:-1],axis=0)
combv = norm_v[0,:,:] + norm_v[1,:,:] - norm_v[2,:,:] - norm_v[3,:,:]

kernel= np.ones((3,3))
kernel[1,1] = 0.3
kernel[0,:] = 0.7/8
kernel[2,:] = 0.7/8
kernel[1,0] = 0.7/8
kernel[1,2] = 0.7/8


#310:350,390:454
smooth = convolve(combq,kernel,mode='constant') 
std_smooth = np.std(smooth[310:350,390:454])

#with holly's imax_find sun, with contours as input
#theta_fit = np.linspace(-np.pi/2+np.pi/50.875,-np.pi/2+np.pi/25,1800)
#xcen = -844.0947807999782-65
#ycen = 14822.481703988178-7
#R_2 = 14716.969579136388

#with inflection from find_limb and that as input to my mod of imax_find_sun
theta_fit = np.linspace(-np.pi/2+np.pi/52.025,-np.pi/2+np.pi/26.5, 1800)
xcen = -1044.9200117954174
ycen = 16862.994425535013
R_2 = 16767.565720491271

x_fit = xcen + (R_2) * np.cos(theta_fit)
y_fit = ycen + (R_2) * np.sin(theta_fit)
limbc = np.array(list(set(list(zip(np.around(y_fit[0:-6]),np.around(x_fit[0:-6])))))).astype(int)

#making a limb mask for all points below the limb
#a grid of indices
indc = np.meshgrid(*map(np.arange,combq.shape),indexing='ij')
#all indices below the limb get a value zero
for i in range(limbc.shape[0]):
    indc[0][0:limbc[i,0]+1,limbc[i,1]]=0
#masking a mask by masking all values zero, gives true for masked values and false for unmasked values
limbmask = np.ma.masked_equal(indc[0],0).mask


#plotting circle fitting the limb on combined LP map or continuum intensity
fig = plt.figure(figsize=(12,12))
ax = plt.axes()
#im=plt.imshow(smooth,cmap = 'gray',vmin=-0.01,vmax=0.01)
im=plt.imshow(n[0,4,:,:],cmap = 'gray')#,vmin=-0.01,vmax=0.01)
fig.colorbar(im)
fig.tight_layout(pad=1.8)
plt.plot(x_fit, y_fit, 'r-')
plt.gca().invert_yaxis()
fig.suptitle('smooth')


#masking smoothed LP map so that all values betweeen -2sigma and 2sigma are masked (including points off limb)
x = np.ma.masked_inside(smooth,-2*std_smooth,2*std_smooth)
#removing pixels off limb from array
x_wolimb= np.ma.array(x,mask=limbmask)
#changing the filled value in masked array from True to Nan
np.ma.set_fill_value(x_wolimb,np.nan)
#obtaining this array with nan values in place of True
filledx = x_wolimb.filled()

#y is the  binary representation of x. ma.masked_inside places True/1 for values which are masked
# and False/0 for values which are not, that is values that remain post masking. Therefore applying logical not on it, to get true for values
#that remain after masking and false for values that are masked and finally converting them to 0 and 1
###includes pixels that are off the limb
y = np.logical_not(np.ma.masked_inside(smooth,-2*std_smooth,2*std_smooth).mask).astype(np.int)




#eroding and dilating Y to keep structures larger than 3x3
###the erode dilate only work on numpy arrays not on masked array
##therefore it worked on le directly because dilate gives an array as output
##to make it work on le_wolimb, a masked array, had to fill the masked values
# and extract the filled array
le1 = erode(y,structure=np.ones((3,3))).astype(x.dtype)
le = dilate(le1,structure=[[True,True,True],[True,True,True],[True,True,True]]).astype(le1.dtype)
le_wolimb = np.ma.array(le,mask=limbmask)
np.ma.set_fill_value(le_wolimb,0.0)
#labeling islands leftover from erode and dilate without off limb pixels
island_lab = measurements.label(le_wolimb.filled())
labels = island_lab[0]


#plotting islands leftover after erode and dilate and pixels not on sun removed
#fig1 = plt.figure(figsize=(12,12))
#ax1 = plt.axes()
##im1 = plt.imshow(x_wolimb,vmax=0.01,vmin=-0.01)
#plt.imshow(labels)
##fig2.colorbar(im1)
#fig1.tight_layout(pad=1.8)
#plt.plot(x_fit, y_fit, 'k')
#plt.gca().invert_yaxis()



levels = 1
fig2 = plt.figure(figsize=(12,12))
ax2 = plt.axes()
im2 = plt.imshow(combq,vmax=0.01,vmin=-0.01)
im2_1 = plt.contour(le_wolimb,levels,origin='lower',colors = 'r')
fig2.colorbar(im2)
fig2.tight_layout(pad=1.8)
plt.plot(x_fit, y_fit, 'b')
plt.gca().invert_yaxis()

size,posi_size,nega_size,avglp,avglpp,avglpn,avgcv,avgcv_plp,avgcv_nlp,avgic,ang = [[] for i in range(11)]

#labels go from 1 to 967 not zero, 

for i in range(1,island_lab[1]+1):
    print(i)
    mk = np.ma.masked_equal(labels,i).mask
    size.append(labels[labels==i].shape[0])
    avglp.append(np.mean(combq[mk]))
    avgcv.append(np.mean(combv[mk]))
    avgic.append(np.mean(norm_ic[mk]))
    ang.append(np.mean(a[mk]))

for i in range(len(size)):
    if size[i]>900 or size[i]<5:
        size[i]=np.nan
        avglp[i]=np.nan
        avgcv[i]=np.nan
        avgic[i]=np.nan

#remove nan, since hist can't handle it
#size_wonan = list(filter(lambda x:not(np.isnan(x)),size))
#
#weights = np.ones_like(size_wonan)/float(len(size_wonan))
#plt.figure()
#plt.hist(size_wonan,50,weights=weights)



#hist for positive and negative polaritites

for i in range(len(size)):
    if avglp[i]>0:
        posi_size.append(size[i])
        avglpp.append(avglp[i])
        avgcv_plp.append(avgcv[i])
    elif avglp[i]<0:
        nega_size.append(size[i])
        avglpn.append(avglp[i])
        avgcv_nlp.append(avgcv[i])

pweights = np.ones_like(posi_size)/len(posi_size)
nweights = np.ones_like(nega_size)/len(nega_size)


plt.figure()
plt.hist([posi_size,nega_size],100,color = ['r','k'])
plt.yscale('log',nonposy='clip')

#polynomials
pp = np.polyfit(posi_size,avglpp,1)
pn = np.polyfit(nega_size,avglpn,1)

xaxis_p = np.unique(np.sort(posi_size))
xaxis_n = np.unique(np.sort(nega_size))

#plt.figure()
#plt.scatter(size,avglp,marker='.')
#plt.ylabel('Avg Linear Pol')
#plt.xlabel('Size in pixels')
#plt.plot(xaxis_p,np.polyval(pp,xaxis_p),'r')
#plt.plot(xaxis_n,np.polyval(pn,xaxis_n),'r')
#
#plt.figure()
#plt.scatter(size,avgcv,marker='.')
#plt.ylabel('Avg Combined V')
#plt.xlabel('Size in pixels')
#
#plt.figure()
#plt.scatter(size,avgic,marker='.')
#plt.ylabel('Avg Continuum intensity')
#plt.xlabel('Size in pixels')
#
#plt.figure()
#plt.scatter(avglp,avgcv,marker='.')
#plt.ylabel('Avg combined V')
#plt.xlabel('Average LP')
#
#plt.figure()
#plt.scatter(avgic,avglp,marker='.')
#plt.ylabel('Avg LP')
#plt.xlabel('Avg I cont')

plt.figure()
plt.scatter(avglp,ang,marker='.')
plt.ylabel('angle')
plt.xlabel('Average LP')

#n,p = [[] for i in range(2)]
#for i in range(len(avglpp)):
#    if avgcv_plp[i]<0:
#        n.append(avglpp[i])
#    else:
#        p.append(avglpp[i])
#    
#plt.figure()
#plt.hist([p,n],200,color=['r','k'])
#plt.xlabel('avg LP of islands, red for patches with +ve Combined V and black for -ve CV ')
#
#del p,n
#
#n,p = [[] for i in range(2)]
#for i in range(len(avglpn)):
#    if avgcv_nlp[i]<0:
#        n.append(avglpn[i])
#    else:
#        p.append(avglpn[i])
#plt.figure()
#plt.hist([p,n],200,color=['r','k'])        
#plt.xlabel('avg LP of islands, red for patches with +ve Combined V and black for -ve CV ')

#%%

#pr[0][0] and pr[0][1] are hists for posi and nega arrays and pr[1] is the bins for them
pr= plt.hist([posi_size,nega_size],bins = 100)
#centre points of bins
bcp = 0.5*(pr[1][1:]+pr[1][:-1])

plt.figure()
plt.plot(bcp,pr[0][0],'r-')
plt.plot(bcp,pr[0][1],'b-')


#separating islands for positive and negative LP signal

yp = np.ma.greater(smooth,2*std_smooth).astype(int)

yn = np.ma.less(smooth,-2*std_smooth).astype(int)


le_p = erode(yp,structure=np.ones((3,3))).astype(x.dtype)
le_n = erode(yn,structure=np.ones((3,3))).astype(x.dtype)
dil_p = dilate(le_p,structure=[[True,True,True],[True,True,True],[True,True,True]]).astype(le_p.dtype)
dil_n = dilate(le_n,structure=[[True,True,True],[True,True,True],[True,True,True]]).astype(le_n.dtype)

le_p_wolimb = np.ma.array(dil_p,mask=limbmask)
le_n_wolimb = np.ma.array(dil_n,mask=limbmask)
np.ma.set_fill_value(le_p_wolimb,0.0)
np.ma.set_fill_value(le_n_wolimb,0.0)

fig3 = plt.figure(figsize=(12,12))
ax3 = plt.axes()
im3 = plt.imshow(combq,vmax=0.01,vmin=-0.01)
im3_1 = plt.contour(le_p_wolimb,levels,origin='lower',colors = 'r')
fig3.colorbar(im3)
fig3.tight_layout(pad=1.8)
plt.plot(x_fit, y_fit, 'b')
plt.gca().invert_yaxis()

#labeling islands leftover from erode and dilate without off limb pixels
island_lab_p = measurements.label(le_p_wolimb.filled())
island_lab_n = measurements.label(le_n_wolimb.filled())
labels_p = island_lab_p[0]
labels_n = island_lab_n[0]

avgq_p1,avgq_p3,avgq_n1,avgq_n3,sizep,sizen,angp,angn = [[] for i in range(8)]

for i in range(1,island_lab_p[1]+1):
    mkp = np.ma.masked_equal(labels_p,i).mask
    sizep.append(labels_p[labels_p==i].shape[0])
    avgq_p1.append(np.mean(n[1,1,:,:][mkp]))
    avgq_p3.append(np.mean(n[1,3,:,:][mkp]))
    angp.append(np.median(a[mkp]))

for i in range(1,island_lab_n[1]+1):
    mkn = np.ma.masked_equal(labels_n,i).mask
    sizen.append(labels_n[labels_n==i].shape[0])
    avgq_n1.append(np.mean(n[1,1,:,:][mkn]))
    avgq_n3.append(np.mean(n[1,3,:,:][mkn]))    
    angn.append(np.median(a[mkn]))

plt.figure()
plt.scatter(sizep,avgq_p1,marker='.')
plt.ylabel('Avg Q, wv 1,posi')
plt.xlabel('Size')


plt.figure()
plt.scatter(sizen,avgq_n1,marker='.')
plt.ylabel('Avg Q, wv 1, nega')
plt.xlabel('Size')




