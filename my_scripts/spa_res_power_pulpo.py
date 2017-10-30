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

figaa = plt.figure(figsize=(12,12))
axaa = plt.axes()
plt.imshow(i_cont,cmap='gray')
plt.gca().invert_yaxis()
dim_i = i_cont.shape
h = signal.hann(i_cont.shape[0])

hann2 = np.outer(h,h)


kx =np.fft.fftfreq(i_cont.shape[0],d=0.05)
kx = fftshift(kx)





dekh1 = PSD2(i_cont)
dekh2 = PSD2(stokv)
# fwn1,zz1 = pspec(dekh1,view=True,wavenumber=True)
# fwn2,zz2 = pspec(dekh2,view=True,wavenumber=True)

four = fftshift(fft2(hann2*i_cont))

po = np.abs(four)**2

four2 = fftshift(fft2(hann2*stokv))
po2 = np.abs(four2)**2

#fig = plt.figure(figsize=(12,12))
#ax = plt.axes()
#img = plt.imshow(po,cmap='gray',norm=mpl.colors.LogNorm())
#plt.gca().invert_yaxis()
#plt.colorbar(img)

# figd = plt.figure(figsize=(12,12))
# axd = plt.axes()
# imgd = plt.imshow(dekh,cmap='gray',norm=mpl.colors.LogNorm())
# plt.gca().invert_yaxis()
# plt.colorbar(imgd)

# fige = plt.figure()
# axe = plt.axes()
# imge1 = plt.loglog(fwn1,zz1,'r')
# imge2 = plt.loglog(fwn2,zz2,'b')




x = np.array(list(range(0,936))*936).reshape(936,936)
y= x.T

distances = np.sqrt((x-467.5)**2+(y-467.5)**2)

distances2 = np.around(distances)



new_po = [np.sum(po[distances2==i]) for i in np.arange(np.min(distances2),np.max(distances2)+1)]
new_po2 = [np.sum(po2[distances2==i]) for i in np.arange(np.min(distances2),np.max(distances2)+1)]
#fig2 = plt.figure()
#ax2 = plt.axes()
# img2 =plt.loglog(np.arange(np.min(distances2),np.max(distances2)+1),new_po,'r')
# img3 = plt.loglog(np.arange(np.min(distances2),np.max(distances2)+1),new_po2,'b')



dd= list(itertools.zip_longest(kx[468::],new_po))
print(len(kx[468::]))
print(len(new_po[0:468]))
ones = np.ones((936,936))

#ones[distances2==50] *=10^7
#img3 = plt.loglog(kx[468::],new_po[0:468],'r')
#img4 = plt.loglog(kx[468::],new_po2[0:468],'b')

# plt.imshow(ones,cmap='gray')
plt.show()
