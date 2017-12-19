#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 17:17:55 2017

@author: prabhu
"""

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


sig = 2

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


        
dim = q.shape
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

smooth = convolve(combq,kernel,mode='constant') 
#std_smooth = np.std(smooth[310:350,390:454])
std_smooth = np.std(smooth[310:350,430:480])

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
im=plt.imshow(n[0,4,:,:],vmin = 0,vmax= 1.5,cmap = 'gray')
fig.colorbar(im)
fig.tight_layout(pad=1.8)
plt.plot(x_fit, y_fit, 'r-')
plt.gca().invert_yaxis()
plt.suptitle('Continuum Intensity')


#masking smoothed LP map so that all values betweeen -2sigma and 2sigma are masked (including points off limb)
x = np.ma.masked_inside(smooth,-sig*std_smooth,sig*std_smooth)
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
y = np.logical_not(np.ma.masked_inside(smooth,-sig*std_smooth,sig*std_smooth).mask).astype(np.int)



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
    ang.append(np.median(a[mk]))
    stda.append(np.std(a[mk]))
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
#separating islands for positive and negative LP signal

yp = np.ma.greater(smooth,sig*std_smooth).astype(int)

yn = np.ma.less(smooth,-sig*std_smooth).astype(int)


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
im3_2 = plt.contour(le_n_wolimb,levels,origin='lower',colors = 'b')
fig3.colorbar(im3)
fig3.tight_layout(pad=1.8)
plt.plot(x_fit, y_fit, 'b')
plt.gca().invert_yaxis()
plt.suptitle('Combined Q')

fig3 = plt.figure(figsize=(12,12))
ax3 = plt.axes()
im3 = plt.imshow(a*180/np.pi)
im3_1 = plt.contour(le_p_wolimb,levels,origin='lower',colors = 'r')
im3_2 = plt.contour(le_n_wolimb,levels,origin='lower',colors = 'b')
fig3.colorbar(im3)
fig3.tight_layout(pad=1.8)
plt.plot(x_fit, y_fit, 'b')
plt.gca().invert_yaxis()
plt.suptitle('Rotation angle (deg)')

#labeling islands leftover from erode and dilate without off limb pixels
island_lab_p = measurements.label(le_p_wolimb.filled())
island_lab_n = measurements.label(le_n_wolimb.filled())
labels_p = island_lab_p[0]
labels_n = island_lab_n[0]

tcv_p,tcv_n,tlp_p,tlp_n,sizep,sizen,angp,angn = [[] for i in range(8)]

for i in range(1,island_lab_p[1]+1):
    mkp = np.ma.masked_equal(labels_p,i).mask
    sizep.append(labels_p[labels_p==i].shape[0])
    angp.append(np.mean(a[mkp]))
    tlp_p.append(np.sum(combq[mkp]))
    tcv_p.append(np.sum(combv[mkp]))
    
for i in range(1,island_lab_n[1]+1):
    mkn = np.ma.masked_equal(labels_n,i).mask
    sizen.append(labels_n[labels_n==i].shape[0])    
    angn.append(np.mean(a[mkn]))
    tlp_n.append(np.sum(combq[mkn]))
    tcv_n.append(np.sum(combv[mkn]))
    
for i in range(len(sizep)):
    if sizep[i]>700 or sizep[i]<5:
        sizep[i] = np.nan
        angp[i] = np.nan
        tlp_p[i] = np.nan
        tcv_p[i] = np.nan
        
for i in range(len(sizen)):
    if sizen[i]>700 or sizen[i]<5:
        sizen[i] = np.nan
        angn[i] = np.nan
        tlp_n[i] = np.nan
        tcv_n[i] = np.nan



#%%
'''PLots that are important
namely :
- size vs angle with mean
- size vs angle with hist
- size vs average linear polarisation 
- size vs total linear polarisation
'''
#####Size vs angle plots
stdb = np.array(stda)*50
plt.figure()        
plt.scatter(size,np.array(ang)*180/np.pi,s=stdb)
plt.xlabel('Size')
plt.ylabel('Angle')
plt.suptitle('Size vs Angle')


angh = []
a_deg = a*180/np.pi
for i in range(1,island_lab[1]+1):
    mk = np.ma.masked_equal(labels,i).mask
    no,bin_edge = np.histogram(a_deg[mk],bins=10)
    bincen = 0.5*(bin_edge[1:]+bin_edge[:-1])
    angh.append(np.mean(bincen[no==np.max(no)]))

   
plt.figure()        
plt.scatter(size,angh,s=stdb)
plt.xlabel('Size')
plt.ylabel('Angle')
plt.suptitle('Size vs Angle with hist')



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

#####Size vs average linear pol
plt.figure()
plt.scatter(size_won,np.abs(avglp_won),marker='.')
plt.ylabel('Avg Linear Pol')
plt.xlabel('Size in pixels')
plt.plot(xaxis,np.polyval(polyp1,xaxis),'r')
#plt.ylim([0.002,0.012])
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
-Size vs ACV
-ALP vs ACV
-Map of angle(without double rotation)
- Angle vs Tlp
- Tlp vs Tcv
- Size vs Tcv
- Angle  vs Tlp with hist

'''
#
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

#plotting histograms, but with lines
#pr[0][0] and pr[0][1] are hists for posi and nega arrays and pr[1] is the bins for them
sizep_won = list(filter(lambda x:not(np.isnan(x)),sizep))
sizen_won = list(filter(lambda x:not(np.isnan(x)),sizen))
plt.figure()
pr= plt.hist([sizep_won,sizen_won],bins = 20)
#centre points of bins
bcp = 0.5*(pr[1][1:]+pr[1][:-1])

plt.clf()
plt.plot(bcp,pr[0][0],'r-')
plt.plot(bcp,pr[0][1],'b-')



fig3 = plt.figure(figsize=(12,12))
ax3 = plt.axes()
im3 = plt.imshow(a*180/np.pi)
im3_1 = plt.contour(le_p_wolimb,levels,origin='lower',colors = 'r')
im3_2 = plt.contour(le_n_wolimb,levels,origin='lower',colors = 'b')
fig3.colorbar(im3)
fig3.tight_layout(pad=1.8)
plt.plot(x_fit, y_fit, 'b')
plt.gca().invert_yaxis()
plt.suptitle('Angle')

plt.figure()        
plt.scatter(np.array(angp)*180/np.pi,tlp_p,color ='r')
plt.scatter(np.array(angn)*180/np.pi,tlp_n,color ='b')
plt.xlabel('angle')
plt.ylabel('total linear pol')
plt.suptitle('Angle vs TLP with mean')

plt.figure()
plt.scatter(tlp_p,tcv_p,marker='.',color='r')
plt.scatter(tlp_n,tcv_n,marker='.',color='b')
plt.ylabel('total combined V')
plt.xlabel('total LP')
plt.suptitle('Tlp vs Tcv')


plt.figure()
plt.scatter(sizep,tcv_p,marker='.',color='r')
plt.scatter(sizen,tcv_n,marker='.',color='b')
plt.ylabel('total combined V')
plt.xlabel('size')
plt.suptitle('Size vs TCV')

#angle with histogram
ang_p,ang_n=[],[]
a_deg = a*180/np.pi
for i in range(1,island_lab_p[1]+1):
    mkp = np.ma.masked_equal(labels_p,i).mask
    no,bin_edge = np.histogram(a_deg[mkp],bins=10)
    bincen = 0.5*(bin_edge[1:]+bin_edge[:-1])
    ang_p.append(np.mean(bincen[no==np.max(no)]))

del no,bin_edge,bincen
    
for i in range(1,island_lab_n[1]+1):
    mkn = np.ma.masked_equal(labels_n,i).mask
    no,bin_edge = np.histogram(a_deg[mkn],bins=10)
    bincen = 0.5*(bin_edge[1:]+bin_edge[:-1])
    ang_n.append(np.mean(bincen[no==np.max(no)]))
   
plt.figure()        
plt.scatter(np.array(ang_p),tlp_p,color ='r')
plt.scatter(np.array(ang_n),tlp_n,color ='b')
plt.xlabel('angle')
plt.ylabel('total linear pol')
plt.suptitle('Angle vs TLP with hist')

#%%

#squared sum linear polarization
def LP_l(l):
    return np.sqrt(np.square(q[l,:,:])+np.square(u[l,:,:]))

lp = LP_l(0) + LP_l(1) +LP_l(2) + LP_l(3) + LP_l(4)

ylp = np.ma.greater(lp,0.007).astype(int)
le_lp = erode(ylp,structure=np.ones((3,3))).astype(x.dtype)
dil_lp = dilate(le_lp,structure=[[True,True,True],[True,True,True],[True,True,True]]).astype(le_p.dtype)

le_lp_wolimb = np.ma.array(dil_lp,mask=limbmask)
np.ma.set_fill_value(le_lp_wolimb,0.0)


levels=1
fig2= plt.figure(figsize=(12,12))
plt.imshow(lp,vmax=0.011,vmin=0)
plt.colorbar()
plt.contour(le_lp_wolimb,levels,origin='lower',colors = 'r')
fig2.tight_layout(pad=1.8)