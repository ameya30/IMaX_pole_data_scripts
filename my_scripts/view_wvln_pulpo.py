import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from natsort import natsorted
from itertools import chain
import glob

files = natsorted(glob.glob("/scratch/prabhu/HollyWaller/IMaX_pole_data_scripts/primary_scripts/saves_Oct11/*.fits"))

cy = int(input("Choose cycle: ") )

st = int(input("Choose stokes: "))

stokes = {0:'I',1:'Q',2:'U',3:'V'}

#opening only the cycle we choose above as an array and putting it into a list fima
fima = [fits.open(i)[0].data for i in files if '_'+str(cy)+'.' in i]

#taking mean of the maximum intensities for various spectral positions for a particular stokes

# cmax = np.mean([np.max(fima[0][st,i,:,:]) for i in range(fima[0].shape[1])])
# cmin = np.mean([np.min(fima[0][st,i,:,:]) for i in range(fima[0].shape[1])])
dim = fima[0].shape
maif = np.zeros(shape=(dim[0],dim[1],dim[2],dim[3]))
if st==0:
	maif[st,:,:,:] = fima[0][st,:,:,:]/np.mean(fima[0][0,4,:,:])
	up,down=1.5,0.5
elif st==1:
	maif[st,:,:,:] = fima[0][st,:,:,:]/fima[0][0,4,:,:]
	up,down=0.04,-0.04
elif st==2:
	maif[st,:,:,:] = fima[0][st,:,:,:]/fima[0][0,4,:,:]
	up,down=0.04,-0.04
else:
	maif[st,:,:,:] = fima[0][st,:,:,:]/fima[0][0,4,:,:]
	up,down=0.08,-0.08


print(fima[0][st,1,0,1],maif[st,1,0,1],fima[0][st,4,3,1],maif[st,4,3,1],fima[0][0,4,0,1],fima[0][0,4,3,1])

frames_1 = [list(range(fima[0].shape[1]))]*20
frames = list(chain.from_iterable(frames_1))

# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=1, metadata=dict(artist='Me'), bitrate=1800)

fig = plt.figure(figsize=(12,12))
ax = plt.axes()
im = ax.imshow(fima[0][st,0,:,:],cmap='gray',vmax=up,vmin=down)
fig.colorbar(im)
fig.tight_layout(pad=1.8)
plt.gca().invert_yaxis()



def animate(i):
	im.set_data(fima[0][st,i,:,:])#,cmap='gray',vmax=340,vmin=-140)
	ax.set_title("cycle number = "+str(cy)+', lambda = ' + str(i) +', stokes = ' +stokes[st])
	return im

ani = animation.FuncAnimation(fig,animate,frames,interval = 1000,blit=False)

# ani.save('cycle_{}_stokes_{}.mp4'.format(str(cy),stokes[st]), writer = writer,dpi=200)



plt.show()
