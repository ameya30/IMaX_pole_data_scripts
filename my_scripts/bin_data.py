from con_re_scipy import congrid
from scipy.io import readsav
from matplotlib import pyplot as plt
import numpy as np
from astropy.io import fits

data = readsav("/home/prabhu/sunrise_holly/movie_data/mk_magneto_tr2_reduc_rnr_300_21.sav",python_dict=True)

iid = data['iid'] #demodulated and restored 
iidn = data['iidn'] #demodulated and not restored

bi = 3
stok,wvln,dimy,dimx = iid.shape

binr = np.zeros(shape = (stok,wvln,int(dimx/bi),int(dimy/bi)))
binnr = np.zeros(shape = (stok,wvln,int(dimx/bi),int(dimy/bi)))

#binning bi x bi for restored data
for i in range(stok):
	for j in range(wvln):
		ima = iid[i,j,:,:]
		binr[i,j,:,:] = congrid(ima,(dimy/bi,dimx/bi))

#binning bi x bi for non-restored data
for i in range(stok):
	for j in range(wvln):
		iman = iidn[i,j,:,:]
		binnr[i,j,:,:] = congrid(iman,(dimy/bi,dimx/bi))



# plt.imshow(binr[0,4,:,:],cmap='gray')
# plt.gca().invert_yaxis()
# plt.show()

#creating primary HDU for binned restored data (primary HDU is mandatory in writing fits file)
hdu1 = fits.PrimaryHDU(data=binr)

#creating Image HDU for binned non-restored data
hdu2 = fits.ImageHDU(data=binnr,name="not restored")

#making hdulist to write into a fits file
hdulist = fits.HDUList([hdu1,hdu2])

hdulist.writeto('/home/prabhu/sunrise_holly/binned_cycles/binned_tr2_mk_restor_300_21_' + str(bi) +'.fits')