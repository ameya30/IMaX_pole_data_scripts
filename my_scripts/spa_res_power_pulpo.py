import numpy as np
from scipy.fftpack import fft2,fftshift
from scipy import signal
from astropy.io import fits
from matplotlib import pyplot as plt
import matplotlib as mpl
import itertools

from power import PSD2,pspec

dum = fits.open("/scratch/prabhu/HollyWaller/IMaX_pole_data_scripts/primary_scripts/saves_Oct11/post_demod_tr2_output_22.fits")

data = dum[0].data

i_cont = data[0,4,:,:]

stokv = data[3,2,:,:]

# stokq =data[1,2,:,:]

##figure to plot zoomed in images of stokes I, V and Q for the meeting with sami
# figaa = plt.figure(figsize=(12,12))
# figaa.suptitle('Stokes Q, lambda2')
# axaa = plt.axes()
# imaaa= plt.imshow(i_cont,cmap='gray')#,vmin=100,vmax=8350)
# imaaa=plt.imshow(stokq[715:785,419:590],cmap='gray',vmin=-62.5,vmax=127)
# plt.gca().invert_yaxis()
# plt.colorbar(imaaa)
# plt.savefig('l_stokesQ_lambda2.png',dpi=500)
# plt.show()




dim_i = i_cont.shape

#making a two dimensional Hann window
h = signal.hann(i_cont.shape[0])
hann2 = np.outer(h,h)

#making a frequnxy axis for the power spectrum, d=0.05 is the spacing in the spatial domain, each pixel for
#imax is 0.05 arc seconds, therefore the spatial frequency will be 1/0.05 arcsec
kx =np.fft.fftfreq(i_cont.shape[0],d=0.05)
kx = fftshift(kx)


#doing the power spectrum using a git package, written in power.py, PSD2 gives
#two D power spectra and pspec gives the radial profile, for which I wrote this script
dekh1 = PSD2(i_cont)
dekh2 = PSD2(stokv)
# fwn1,zz1 = pspec(dekh1,view=True,wavenumber=True)
# fwn2,zz2 = pspec(dekh2,view=True,wavenumber=True)

#calculation the fourier transform, shifting frequencies and estimating power spectrum
four = fftshift(fft2(hann2*i_cont))
po = np.abs(four)**2

four2 = fftshift(fft2(hann2*stokv))
po2 = np.abs(four2)**2


#plotting power spectrum
fig = plt.figure(figsize=(12,12))
ax = plt.axes()
img = plt.imshow(po,cmap='gray',norm=mpl.colors.LogNorm())
plt.gca().invert_yaxis()
plt.colorbar(img)
# plt.savefig('powerSpec.png',dpi=500)

#code for plotting power spectrum from the power.py git package
# figd = plt.figure(figsize=(12,12))
# axd = plt.axes()
# imgd = plt.imshow(dekh,cmap='gray',norm=mpl.colors.LogNorm())
# plt.gca().invert_yaxis()
# plt.colorbar(imgd)

# fige = plt.figure()
# axe = plt.axes()
# imge1 = plt.loglog(fwn1,zz1,'r')
# imge2 = plt.loglog(fwn2,zz2,'b')



#making an array of all the xindices of a 936*036 array
x = np.array(list(range(0,936))*936).reshape(936,936)
#array for y indices
y = x.T

#estimating the distance of each element of the x and y array from the centre,
#to get an array of distances, to then estimate equi distant points for the cirlce
distances = np.sqrt((x**2+(y)**2)
# rounding them up as distances
distances2 = np.around(distances)


#obtaining the radial profile at a distance d i.e at a certain spatial frequency and summing up along the circle
#the values at a certain distance/radius are obtained with a boolean mask, distances2==i
#gives a boolean array, and that applied on po or po2 gives an array for which this boolean mask
#holds the values are true
new_po = [np.sum(po[distances2==i]) for i in np.arange(np.min(distances2),np.max(distances2)+1)]
new_po2 = [np.sum(po2[distances2==i]) for i in np.arange(np.min(distances2),np.max(distances2)+1)]


#plotting the radial power spectra
fig2 = plt.figure()
ax2 = plt.axes()
img2 = plt.loglog(kx[468::],new_po[0:468],'r')
img3 = plt.loglog(kx[468::],new_po2[0:468],'b')
plt.xlabel('1/arcsec')
# plt.savefig('radial_profile_powerSpec.png',dpi=500)
plt.show()