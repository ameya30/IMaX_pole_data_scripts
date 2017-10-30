import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt




fima = fits.open('/scratch/prabhu/HollyWaller/IMaX_pole_data_scripts/primary_scripts/saves_Oct11/post_demod_tr2_output_21.fits')[0].data 


st = int(input("Choose stokes: "))

stokes = {0:'I',1:'Q',2:'U',3:'V'}
dim = fima.shape

print(dim)

maif = np.zeros(shape=(dim[0],dim[1],dim[2],dim[3]))
print(maif.shape)
if st==0:
	maif[st,:,:,:] = fima[st,:,:,:]/np.mean(fima[0,4,:,:])
	up,down=1.5,0.5
elif st==1:
	maif[st,:,:,:] = fima[st,:,:,:]/fima[0,4,:,:]
	up,down=0.04,-0.04
elif st==2:
	maif[st,:,:,:] = fima[st,:,:,:]/fima[0,4,:,:]
	up,down=0.04,-0.04
else:
	maif[st,:,:,:] = fima[st,:,:,:]/fima[0,4,:,:]
	up,down=0.08,-0.08




fig = plt.figure(figsize=(12,12))
ax = plt.axes()
im = ax.imshow(maif[st,4,:,:],cmap='gray',vmax=up,vmin=down)
fig.colorbar(im)
fig.tight_layout(pad=1.8)
plt.gca().invert_yaxis()


plt.show()

y1 = 200
y2 = 280
x1 = 340
x2 = 500



std = np.std(maif[st,4,y1:y2,x1:x2])
meanie = np.mean(maif[st,4,y1:y2,x1:x2])
rms = std/meanie
print("rms is {}".format(rms))
print("std is {}".format(std))