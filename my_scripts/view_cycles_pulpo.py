import os
import glob


import matplotlib.pyplot as plt
import matplotlib.animation as animation

import numpy as np
from scipy.io import readsav
from astropy.io import fits
from natsort import natsorted


########for pulpo##########

#####reading mk_restor files#####

##mine
#files = natsorted(glob.glob('/scratch/prabhu/HollyWaller/IMaX_pole_data_scripts/primary_scripts/saves_Oct11/mk_magneto_tr2*'))

##holly
#files = natsorted(glob.glob('/scratch/prabhu/HollyWaller/IMaX_pole_data_scripts/mk_magneto_outputreduc_rnr_/*.sav'))

#files = list(filter(lambda x: x[-3:-1]=='sa',files))
#data = [readsav(i,python_dict = True)['iid'] for i in files]


##### to read in .fits from holly's post demod fringe removal#####

##mine
files = natsorted(glob.glob('/scratch/prabhu/HollyWaller/IMaX_pole_data_scripts/primary_scripts/saves_Oct11/*.fits')) 
#data = [fits.open(i)[0].data for i in files]

##holly
#files = natsorted(glob.glob('/scratch/prabhu/HollyWaller/IMaX_pole_data_scripts/post_demod_output_/*.fits')) 
data = [fits.open(i)[0].data for i in files]



lam = 2
pol = 1

print(data[7].shape)
print (data[7][pol,lam,:,:].max())
print (data[7][pol,lam,:,:].min())


fig = plt.figure(figsize=(12,12))
ax = plt.axes()
im = ax.imshow(data[0][pol,lam,:,:],cmap='gray',vmax=175,vmin=-85)
fig.colorbar(im)
fig.tight_layout(pad=1.8)
plt.gca().invert_yaxis()

def init():
	im.set_data(np.ma.array(data[0][3,2,:,:]))#,mask=True),cmap='gray')
	return im,


def animate(i):
	im.set_data(data[i][pol,lam,:,:])#,cmap='gray',vmax=90,vmin=-90)
	ax.set_title("cycle number "+str(i+7))
	return im,



ani = animation.FuncAnimation(fig,animate,np.arange(len(data)),interval = 1000,blit = False)

ani

plt.show()

