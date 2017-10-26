import numpy as np
from astropy.io import fits
from scipy.fftpack import fft2,fftshift
from scipy import signal
from matplotlib import pyplot as plt
import matplotlib as mpl
def four(x,win):
	power = np.abs(fftshift(fft2(x*win)))**2
	return power

datae = fits.open("/home/prabhu/sunrise_holly/test_wt_diff_flats/mera_26.fits")

bus = [datae[i].data for i in range(1,len(datae))]
a0,a15,a19 = bus[0],bus[15],bus[19]

h = signal.hann(a0.shape[0])
hann2 = np.outer(h,h)
plt.imshow(hann2)
plt.imshow(hann2,cmap='gray',norm=mpl.colors.LogNorm())
pow0 = four(a0,hann2)
pow15 = four(a15,hann2)
pow19 = four(a19,hann2)

fig,ax = plt.subplots()
p0 = ax.imshow(pow0,cmap='gray',norm=mpl.colors.LogNorm(vmin = pow0.min(),vmax=pow0.max()))
fig.colorbar(p0, ax=ax, extend='max')
ax.set_title("cycle 3, image0")

fig1,ax1 = plt.subplots()
p15 = ax1.imshow(pow15,cmap='gray',norm=mpl.colors.LogNorm(vmin = pow15.min(),vmax=pow15.max()))
fig1.colorbar(p15, ax=ax1, extend='max')
ax1.set_title("cycle 3, image15")

fig2,ax2 = plt.subplots()
p19 = ax2.imshow(pow19,cmap='gray',norm=mpl.colors.LogNorm(vmin = pow19.min(),vmax=pow19.max()))
fig2.colorbar(p19, ax=ax2, extend='max')
ax2.set_title("cycle 3, image19")

plt.show()