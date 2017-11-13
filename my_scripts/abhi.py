# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 09:04:01 2017

@author: prabhu
"""

import glob
from natsort import natsorted
import numpy               as np
import matplotlib.pyplot   as plt

from astropy.io     import fits
cy = input("Choose cycle: ")
wv = int(input("Choose wv: "))


normal = "/scratch/prabhu/backup_workstation/sunrise_holly/binned_cycles/binned_tr2_mk_restor_300_"+cy.split("_")[-1]+".fits"
filenameQ = "/scratch/prabhu/backup_workstation/sunrise_holly/binned_cycles/imax_lp_max_" + cy
filenamea = "/scratch/prabhu/backup_workstation/sunrise_holly/binned_cycles/imax_roat_angle_Q_U_" + cy
n= fits.open(normal)[1].data
q = fits.open(filenameQ)[0].data
a = fits.open(filenamea)[0].data
u = fits.open(filenameQ)[1].data

w = q[0,:,:]+q[1,:,:]+q[2,:,:]+q[3,:,:]

print (w.shape)


fig = plt.figure(figsize=(12,12))
ax = plt.axes()
im = ax.imshow(w/n[0,4,:,:],cmap='gray',vmin=-0.01,vmax=0.01)
fig.colorbar(im)
fig.tight_layout(pad=1.8)
plt.gca().invert_yaxis()


fig3 = plt.figure(figsize=(12,12))
ax3 = plt.axes()
im3 = ax3.imshow(u[wv,:,:]/n[0,4,:,:],cmap='gray',vmin=-0.01,vmax=0.01)
fig3.colorbar(im3)
fig3.tight_layout(pad=1.8)
plt.gca().invert_yaxis()

fig2 = plt.figure(figsize=(12,12))
ax2 = plt.axes()
im2 = ax2.imshow(a[:,:],cmap='gray')
fig2.colorbar(im2)
fig2.tight_layout(pad=1.8)
plt.gca().invert_yaxis()

plt.show()
