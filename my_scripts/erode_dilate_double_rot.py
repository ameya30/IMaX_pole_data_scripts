#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 11:31:16 2017

@author: prabhu
"""

from matplotlib import pyplot as plt
import matplotlib as mpl
import cmocean

import numpy as np

from astropy.io import fits
from matplotlib import pyplot as plt
from scipy.io import readsav
from scipy.ndimage.morphology import binary_dilation as dilate
from scipy.ndimage.morphology import binary_erosion as erode
from scipy.ndimage import measurements 
from scipy.ndimage import convolve
from scipy.ndimage.morphology import binary_dilation as dilate





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


cdict = {'red': ((0.0,  0.50, 0.50),
                  (1.0,  1.0, 1.0)),

        'green':  ((0.0,  0.50, 0.50),
                   (1.0,  1.0, 1.0)),
          
        'blue':   ((0.0,  0.50, 0.50),
                  (1.0,  1.0, 1.0))}
our_gray_cmap = mpl.colors.LinearSegmentedColormap('our_gray',cdict,256)
mpl.cm.register_cmap(cmap=our_gray_cmap)


sig = 2.5

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
    rnormal = "/home/prabhu/sunrise_holly/mean_rem_output/mean_rem_output_" + cy + '.fits'
    normal = "/home/prabhu/sunrise_holly/mean_rem_output/mean_rem_output_nr_" + cy + '.fits'
    filenameQ = "/home/prabhu/sunrise_holly/imax_lp_max_/tr_imax_lp_max_norm_nr_" + cy+ '.fits'
    filenamea = "/home/prabhu/sunrise_holly/imax_lp_max_/tr_imax_roat_angle_Q_U_norm_nr_" + cy + '.fits'
	
    r = fits.open(rnormal)[0].data
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



###########
xCen = -1044.9200117954174
yCen = 16862.994425535013
radSun = 16767.565720491271

smooth = combu_p

dimRes = np.shape(smooth)
no_slices = 12
mu_arr = np.zeros(shape = (dimRes[0], dimRes[1]))



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
    nr = rr[(rr>(m+3*s)) | (rr<(m-3*s))]
    while nr.shape[0]:
        m = np.mean(rr)
        s = np.std(rr)
        nr = rr[(rr>(m+3*s)) | (rr<(m-3*s))]
        rr = rr[(rr<(m+3*s)) & (rr>(m-3*s))]
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
plt.figure()
plt.plot(xax, fit_fn(xax))
plt.xlabel('Mu')
plt.ylabel('Noise level')


##use in case of estimating noise level on from squareroot of stokes I
#kyan=fit_fn(xax)
#factor = np.mean(kya/kyan)
#plt.figure()
#plt.plot(xax,kya/factor,'g')
#plt.plot(xax,kyan)
#
#
#polyFit2 = np.polyfit(xax, kya/factor, 6)
#fit_fn2  = np.poly1d(polyFit2)
#plt.figure()
#plt.plot(xax, fit_fn(xax))
#plt.plot(xax,fit_fn2(xax),'r')



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
theta_fit = np.linspace(-np.pi/2+np.pi/50.025,-np.pi/2+np.pi/26.7, 1800)
xCen = -1044.9200117954174
yCen = 16862.994425535013
radSun = 16767.565720491271





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

kernel= np.ones((3,3))
kernel[1,1] = 0.5
kernel[0,:] = 0.5/8
kernel[2,:] = 0.5/8
kernel[1,0] = 0.5/8
kernel[1,2] = 0.5/8
#
smooth = convolve(combq,kernel,mode='constant') 

thres = fit_fn(mu_copy)
del smooth
smooth = convolve(combq,kernel,mode='constant') 
smooth[smooth<(sig*thres)]=90000 #for combq

x_fit = xCen + (radSun) * np.cos(theta_fit)
y_fit = yCen + (radSun) * np.sin(theta_fit)
limbc = np.array(list(set(list(zip(np.around(y_fit[0:-6]),np.around(x_fit[0:-6])))))).astype(int)

indc = np.meshgrid(*map(np.arange,combq.shape),indexing='ij')


for i in range(limbc.shape[0]):
    indc[0][0:limbc[i,0]+20,limbc[i,1]]=0
limbmask = np.ma.masked_equal(indc[0],0).mask


lab = np.arange(-200,900,200)*0.05
ticklab = [str(i) for i in lab.astype(int)]

#plotting circle fitting the limb on combined LP map or continuum intensity
fig = plt.figure(figsize=(12,12))
ax = plt.axes()
im=plt.imshow(n[0,4,:,:],vmin = 0,vmax= 1.5,cmap = 'gray')
rr = fig.colorbar(im)
fig.tight_layout()
plt.plot(x_fit, y_fit, 'k')
plt.gca().invert_yaxis()
rr.set_label('Normalized Continuum Intensity')
ax.set_xticklabels(ticklab)
ax.set_yticklabels(ticklab)
plt.xlabel('[arcsec]')
plt.ylabel('[arcsec]')

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

lad = np.copy(labels)
for i in range(1,island_lab[1]+1):
    mk = np.ma.masked_equal(lad,i).mask
    temp = combq[mk]
    lim = 3.5*np.mean(thres[mk])
    ww = temp[temp>lim]
    if ww.shape[0]:
        lad[lad==i]=1
    else:
        lad[lad==i]=0

del labels, island_lab
island_lab = measurements.label(lad)
labels = island_lab[0]


levels = 1
fig2 = plt.figure(figsize=(12,12))
ax2 = plt.axes()
im2 = plt.imshow(combq,vmax=0.01,vmin=-0.01,cmap='gray')
im2_1 = plt.contour(le_wolimb,levels,origin='lower',colors = 'r')
r2 = fig2.colorbar(im2)
fig2.tight_layout()
plt.plot(x_fit, y_fit, 'k')
plt.gca().invert_yaxis()
r2.set_label('Combined Linear Polarization')
#ax2.set_xticklabels(ticklab)
#ax2.set_yticklabels(ticklab)
plt.xlabel('[arcsec]')
plt.ylabel('[arcsec]')
#plt.savefig('tlp_noiseU_sig2.5.png')

fig = plt.figure(figsize=(12,12))
ax = plt.axes()
im=plt.imshow(r[0,4,:,:],vmin = 0,vmax= 1.5,cmap = 'gray')
im2_1 = plt.contour(le_wolimb,levels,origin='lower',colors = 'r')
rr = fig.colorbar(im)
fig.tight_layout()
plt.plot(x_fit, y_fit, 'k')
plt.gca().invert_yaxis()
rr.set_label('Normalized Continuum Intensity')
#ax.set_xticklabels(ticklab)
#ax.set_yticklabels(ticklab)
plt.xlabel('[arcsec]')
plt.ylabel('[arcsec]')
#plt.savefig('Icont_noiseU_sig2.5.png')


levels = 1
fig2 = plt.figure(figsize=(12,12))
ax2 = plt.axes()
im2 = plt.imshow(combv,vmax = 0.03,vmin=-0.03,cmap='gray')
im2_1 = plt.contour(le_wolimb,levels,origin='lower',colors = 'r')
r2 = fig2.colorbar(im2)
fig2.tight_layout()
plt.plot(x_fit, y_fit, 'k')
plt.gca().invert_yaxis()
r2.set_label('Combined Circular Polarization')
#ax2.set_xticklabels(ticklab)
#ax2.set_yticklabels(ticklab)
plt.xlabel('[arcsec]')
plt.ylabel('[arcsec]')
#plt.savefig('tcv_noiseU_sig2.5.png')


size,psize,nsize,avglp,avgcv,avgic,ang,tlp,tcv,stda,angh= [[] for i in range(11)]

#labels go from 1 to 967 not zero, 

for i in range(1,island_lab[1]+1):
    mk = np.ma.masked_equal(labels,i).mask
    size.append(labels[labels==i].shape[0])
    avglp.append(np.mean(combq[mk]))
    avgcv.append(np.mean(combv[mk]))
    avgic.append(np.mean(norm_ic[mk]))
    ang.append(np.median(a_copy[mk]))
    angh.append(a_copy[mk])
    stda.append(np.std(a_copy[mk]))
    tlp.append(np.sum(combq[mk]))
    tcv.append(np.sum(combv[mk]))
#
#for i in range(len(size)):
#    if size[i]>700 or size[i]<5:
#        size[i] = np.nan
#        avglp[i] = np.nan
#        avgcv[i] = np.nan
#        avgic[i] = np.nan
#        ang[i] = np.nan
#        angh[i]= np.nan
#        stda[i] = np.nan
#        tlp[i] = np.nan
#        tcv[i] = np.nan
#%%

'''PLots that are important namely :
- size vs angle with mean
- size vs angle with hist
- histogram of angles 
    -big
    -small
    -all
- size vs average linear polarisation 
- size vs total linear polarisation
- sqrt(size) vs total linear polarisation
'''
pixel_area = 0.05**2#in arcsec


stdb = np.array(stda)*50
plt.figure()        
plt.scatter(size,ang)#,s=stdb)
plt.xlabel('Size')
plt.ylabel('Angle')
plt.suptitle('Size vs Angle (angle mean)')

plt.figure()        
plt.scatter(size,stda)#,s=stdb)
plt.xlabel('Size')
plt.ylabel('Std')
plt.suptitle('Size vs std of Angle')


#a_hist = list(a_copy.reshape([1,-1]))
a_hist = np.concatenate(np.array(angh))
plt.figure()
plt.hist(a_hist,bins = 100)
plt.xlabel('Angle')
plt.suptitle('Histogram of angles')



#fitting trend to absolute value of LP,
#already attempted a trend line separately
#+ve and -ve LP separately, the trend seems to be
#mirrored, hence will proceed to fit line to abs LP value
size_won = np.log(np.array(list(filter(lambda x:not(np.isnan(x)),size)))) #removing nans, polyfit doesnt like them much
avglp_won = list(filter(lambda x:not(np.isnan(x)),avglp))
tlp_won = list(filter(lambda x:not(np.isnan(x)),tlp))
polyp1 = np.polyfit(size_won,np.abs(avglp_won),3)
polyp2 = np.polyfit(size_won,np.abs(tlp_won),3)
polyp3 = np.polyfit(np.sqrt(size_won),np.abs(tlp_won),3)
xaxis = np.unique(np.sort(size_won))


plt.figure()
plt.scatter(size_won,avglp_won,marker='.')
plt.ylabel('Avg Linear Pol')
plt.xlabel('Size in # of pixels')
plt.plot(xaxis,np.polyval(polyp1,xaxis),'r')
plt.tight_layout()
#plt.ylim([0.0029,0.008])
#plt.savefig('avg_size_noiseU_sig2.5.png')


plt.figure()
plt.scatter(size_won,np.abs(tlp_won),marker='.')
plt.ylabel('Total Linear Pol')
plt.xlabel('Size in # of pixels')
plt.plot(xaxis,np.polyval(polyp2,xaxis),'r')
plt.tight_layout()
#plt.savefig('tlp_size_noiseU_sig2.5.png')

plt.figure()
plt.scatter(np.sqrt(size_won),np.abs(tlp_won),marker='.')
plt.plot(np.sqrt(xaxis),np.polyval(polyp3,np.sqrt(xaxis)),'r')
plt.ylabel('Total Linear Pol')
plt.xlabel('sqrt of size in pixels')
plt.tight_layout()
#plt.savefig('tlp_sqsize_noiseU_sig2.5.png')

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
im3 = plt.imshow(a_copy,cmap = cmocean.cm.phase)
im3_1 = plt.contour(le_wolimb,levels,origin='lower',colors = 'k')
r3 = fig3.colorbar(im3)
fig3.tight_layout(pad=1.8)
plt.plot(x_fit, y_fit, 'b')
plt.gca().invert_yaxis()
r3.set_label('Azimuth Angle [deg]')
ax3.set_xticklabels(ticklab)
ax3.set_yticklabels(ticklab)
plt.xlabel('[arcsec]')
plt.ylabel('[arcsec]')
#plt.savefig('anglemap_noiseU_sig2.5.png')


plt.figure()
plt.scatter(tlp,tcv,marker='.')
plt.ylabel('total combined V')
plt.xlabel('total LP')
plt.suptitle('Tlp vs Tcv')


plt.figure()
plt.scatter(size,np.abs(tcv),marker='.')
plt.ylabel('total combined V')
plt.xlabel('size')
plt.suptitle('Size vs abs(TCV)')


#%%


limits = list(range(0,240,40)) + [300,400,4500]
lab_c1 = np.copy(labels)
lab_c2 = np.copy(labels)
lab_c3 = np.copy(labels)
lab_c4 = np.copy(labels)
lab_c5 = np.copy(labels)
lab_c6 = np.copy(labels)
lab_c7 = np.copy(labels)
lab_c8 = np.copy(labels)

a_1,a_2,a_3,a_4,a_5,a_6,a_7,a_8=[[] for i in range(8)]

        
for i in range(1,island_lab[1]+1):
    ss = labels[labels==i].shape[0]
    mk = np.ma.masked_equal(labels,i).mask
    if ss<=60 and ss>0:
        lab_c1[labels==i]=101010
        a_1.append(a_copy[mk])
    elif ss<=400 and ss>60:
        lab_c2[labels==i]=101010
        a_2.append(a_copy[mk])
    elif ss<=4500 and ss>400:
        lab_c3[labels==i]=101010
        a_3.append(a_copy[mk])
#    elif ss<=5400 and ss>500:
#        lab_c4[labels==i]=101010
#        a_4.append(a_copy[mk])
#    elif ss<=4500 and ss>500:
#        lab_c5[labels==i]=101010
#        a_5.append(a_copy[mk])
#    elif ss<=300 and ss>200:
#        lab_c6[labels==i]=101010
#        a_6.append(a_copy[mk])
#    elif ss<=400 and ss>300:
#        lab_c7[labels==i]=101010
#        a_7.append(a_copy[mk])
#    elif ss<=4500 and ss>400:
#        lab_c8[labels==i]=101010
#        a_8.append(a_copy[mk])

lab_c1[lab_c1 != 101010] = False
lab_c1[lab_c1 == 101010] = True
lab_c2[lab_c2 != 101010] = False
lab_c2[lab_c2 == 101010] = True
lab_c3[lab_c3 != 101010] = False
lab_c3[lab_c3 == 101010] = True
#lab_c4[lab_c4 != 101010] = False
#lab_c4[lab_c4 == 101010] = True
#lab_c5[lab_c5 != 101010] = False
#lab_c5[lab_c5 == 101010] = True
#lab_c6[lab_c6 != 101010] = False
#lab_c6[lab_c6 == 101010] = True
#lab_c7[lab_c7 != 101010] = False
#lab_c7[lab_c7 == 101010] = True
#lab_c8[lab_c8 != 101010] = False
#lab_c8[lab_c8 == 101010] = True

a1_hist = np.concatenate(np.array(a_1))
a2_hist = np.concatenate(np.array(a_2))
a3_hist = np.concatenate(np.array(a_3))
#a4_hist = np.concatenate(np.array(a_4))
#a5_hist = np.concatenate(np.array(a_5))
#a6_hist = np.concatenate(np.array(a_6))
#a7_hist = np.concatenate(np.array(a_7))
#a8_hist = np.concatenate(np.array(a_8))


levels = 1
fig2 = plt.figure(figsize=(12,12))
ax2 = plt.axes()
im2 = plt.imshow(combq,vmax=0.01,vmin=-0.01,cmap='gray')
im2_1 = plt.contour(lab_c1,levels,origin='lower',colors = 'r')
rr = fig2.colorbar(im2)
fig2.tight_layout()
rr.set_label('Combined Linear Polarization')
ax2.set_xticklabels(ticklab)
ax2.set_yticklabels(ticklab)
plt.plot(x_fit, y_fit, 'b')
plt.gca().invert_yaxis()
plt.xlabel('[arcsec]')
plt.ylabel('[arcsec]')
#plt.savefig('lab_c1map_noiseU_sig2.5.png')

plt.figure()
plt.hist(a1_hist,bins=25)
plt.xlabel('azimuth angle [degrees]')
plt.ylabel('Number of pixels')
plt.tight_layout()
#plt.savefig('lab_c1_noiseU_sig2.5.png')

levels = 1
fig2 = plt.figure(figsize=(12,12))
ax2 = plt.axes()
im2 = plt.imshow(combq,vmax=0.01,vmin=-0.01,cmap='gray')
im2_1 = plt.contour(lab_c2,levels,origin='lower',colors = 'r')
rr = fig2.colorbar(im2)
fig2.tight_layout()
rr.set_label('Combined Linear Polarization')
ax2.set_xticklabels(ticklab)
ax2.set_yticklabels(ticklab)
plt.plot(x_fit, y_fit, 'b')
plt.gca().invert_yaxis()
plt.xlabel('[arcsec]')
plt.ylabel('[arcsec]')
#plt.savefig('lab_c2map_noiseU_sig2.5.png')

plt.figure()
plt.hist(a2_hist,bins=25)
plt.xlabel('azimuth angle [degrees]')
plt.ylabel('Number of pixels')
plt.tight_layout()
#plt.savefig('lab_c2_noiseU_sig2.5.png')




levels = 1
fig2 = plt.figure(figsize=(12,12))
ax2 = plt.axes()
im2 = plt.imshow(combq,vmax=0.01,vmin=-0.01,cmap='gray')
im2_1 = plt.contour(lab_c3,levels,origin='lower',colors = 'r')
rr = fig2.colorbar(im2)
fig2.tight_layout()
rr.set_label('Combined Linear Polarization')
ax2.set_xticklabels(ticklab)
ax2.set_yticklabels(ticklab)
plt.plot(x_fit, y_fit, 'b')
plt.gca().invert_yaxis()
plt.xlabel('[arcsec]')
plt.ylabel('[arcsec]')
#plt.savefig('lab_c3map_noiseU_sig2.5.png')
#


plt.figure()
plt.hist(a3_hist,bins=25)
plt.xlabel('azimuth angle [degrees]')
plt.ylabel('Number of pixels')
plt.tight_layout()
#plt.savefig('lab_c3_noiseU_sig2.5.png')




#levels = 1
#fig2 = plt.figure(figsize=(12,12))
#ax2 = plt.axes()
#im2 = plt.imshow(combq,vmax=0.01,vmin=-0.01,cmap='gray')
#im2_1 = plt.contour(lab_c4,levels,origin='lower',colors = 'r')
#rr = fig2.colorbar(im2)
#fig2.tight_layout()
#rr.set_label('Combined Linear Polarization')
#ax2.set_xticklabels(ticklab)
#ax2.set_yticklabels(ticklab)
#plt.plot(x_fit, y_fit, 'b')
#plt.gca().invert_yaxis()
#plt.xlabel('[arcsec]')
#plt.ylabel('[arcsec]')
#plt.savefig('lab_c4map_noiseU_sig2.5.png')
#
#
#
#plt.figure()
#plt.hist(a4_hist,bins=25)
#plt.xlabel('azimuth angle [degrees]')
#plt.ylabel('Number of pixels')
#plt.tight_layout()
#plt.savefig('lab_c4_noiseU_sig2.5.png')



#
#levels = 1
#fig2 = plt.figure(figsize=(12,12))
#ax2 = plt.axes()
#im2 = plt.imshow(combq,vmax=0.01,vmin=-0.01,cmap='gray')
#im2_1 = plt.contour(lab_c5,levels,origin='lower',colors = 'r')
#rr = fig2.colorbar(im2)
#fig2.tight_layout()
#rr.set_label('Combined Linear Polarization')
#ax2.set_xticklabels(ticklab)
#ax2.set_yticklabels(ticklab)
#plt.plot(x_fit, y_fit, 'b')
#plt.gca().invert_yaxis()
#plt.xlabel('[arcsec]')
#plt.ylabel('[arcsec]')
#plt.savefig('lab_c5map_noiseU_sig2.5.png')
#
#
#plt.figure()
#plt.hist(a5_hist,bins=25)
#plt.xlabel('azimuth angle [degrees]')
#plt.ylabel('Number of pixels')
#plt.tight_layout()
#plt.savefig('lab_c5_noiseU_sig2.5.png')




