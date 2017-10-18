import os
import glob


import matplotlib.pyplot as plt
import matplotlib.animation as animation

import numpy as np
from scipy.io import readsav
from natsort import natsorted

files = natsorted(glob.glob('/home/prabhu/sunrise_holly/movie_data/*.sav'))

data = [readsav(i,python_dict = True)['iid'] for i in files]

lam = 4
pol = 0

print (data[7][pol,lam,:,:].max())
print (data[7][pol,lam,:,:].min())


fig = plt.figure(figsize=(12,12))
ax = plt.axes()
im = ax.imshow(data[0][pol,lam,:,:],cmap='gray',vmax=8250,vmin=75)
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

