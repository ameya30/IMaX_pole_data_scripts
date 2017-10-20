import os
import glob


import matplotlib.pyplot as plt
import matplotlib.animation as animation

import numpy as np
from scipy.io import readsav
from natsort import natsorted

#for workstation
#files = natsorted(glob.glob('/home/prabhu/sunrise_holly/movie_data/*.sav'))

#for pulpo
#mine
#files = natsorted(glob.glob('/scratch/prabhu/HollyWaller/IMaX_pole_data_scripts/primary_scripts/saves_Oct11/mk_magneto_tr2*'))
#holly
files = natsorted(glob.glob('/scratch/prabhu/HollyWaller/IMaX_pole_data_scripts/mk_magneto_outputreduc_rnr_/*.sav'))

files = list(filter(lambda x: x[-3:-1]=='sa',files))

data = [readsav(i,python_dict = True)['iid'] for i in files]

lam = 4
pol = 0

print (data[7][pol,lam,:,:].max())
print (data[7][pol,lam,:,:].min())


fig = plt.figure(figsize=(12,12))
ax = plt.axes()
im = ax.imshow(data[0][pol,lam,:,:],cmap='gray',vmax=8300,vmin=100)
fig.colorbar(im)
fig.tight_layout(pad=1.8)
plt.gca().invert_yaxis()

def init():
	im.set_data(np.ma.array(data[0][3,2,:,:]))#,mask=True),cmap='gray')
	return im,

# im = plt.imshow(data[0][3,2,:,:],cmap='gray',vmax=340,vmin=-140)

# im = plt.imshow(np.ma.array(data[0][3,2,:,:],mask=True),cmap='gray')

def animate(i):
	im.set_data(data[i][pol,lam,:,:])#,cmap='gray',vmax=340,vmin=-140)
	ax.set_title("cycle number "+str(i+7))
	return im,


# ani = animation.FuncAnimation(fig,animate,np.arange(len(data)),init_func = init,interval = 300,blit = True)

ani = animation.FuncAnimation(fig,animate,np.arange(len(data)),interval = 1000,blit = False)

plt.show()

