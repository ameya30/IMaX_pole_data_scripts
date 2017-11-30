import glob
from natsort import natsorted
import numpy               as np
import matplotlib.pyplot   as plt


from astropy.io     import fits

from scipy.ndimage.morphology import binary_dilation as dilate

from scipy.ndimage.morphology import binary_erosion as erode


####y1,y2,x1,x2 = 120,200,225,400 (80*125) dim box

cy = input("Choose cycle: ")

bi = input("Choose binning factor: ")

re = input("Choose r for restored or nr for nonrestored: ")

#opening the chosen restored and non restored data 
if bi=='2':
    y1,y2,x1,x2 = 120,170,125,400
elif bi=='3':
    y1,y2,x1,x2 = 125,145,130,190
elif bi=='4':
    y1,y2,x1,x2 = 80,105,110,150

    

if re=='r':
	normal = "/home/prabhu/sunrise_holly/binned_cycles/binned_tr2_mk_restor_300_"+cy+'_'+bi+'.fits'
	filenameQ = "/home/prabhu/sunrise_holly/imax_lp_max_/imax_lp_max_" + cy+'_'+bi+'.fits'
	filenamea = "/home/prabhu/sunrise_holly/imax_lp_max_/imax_roat_angle_Q_U_" + cy+'_'+bi+'.fits'
	

	n = fits.open(normal)[0].data #0 for restored data 1 for non restored data
	q = fits.open(filenameQ)[0].data #0 is max linear pol map, of maxed Q
	u = fits.open(filenameQ)[1].data #1 is minimized U
	a = fits.open(filenamea)[0].data

elif re=='nr':
	normal = "/home/prabhu/sunrise_holly/binned_cycles/binned_tr2_mk_restor_300_"+cy+'_'+bi+'.fits'
	filenameQ = "/home/prabhu/sunrise_holly/imax_lp_max_/imax_lp_max_nr_" + cy+'_'+bi+'.fits'
	filenamea = "/home/prabhu/sunrise_holly/imax_lp_max_/imax_roat_angle_Q_U_nr_" + cy+'_'+bi+'.fits'
	

	n = fits.open(normal)[1].data #0 for restored data 1 for non restored data
	q = fits.open(filenameQ)[0].data #0 is max linear pol map, of maxed Q
	u = fits.open(filenameQ)[1].data #1 is minimized U
	a = fits.open(filenamea)[0].data

#normalizing the linear polarization array with quiet sun I conti intensity fpr the same box as below for q conti
dim = q.shape
norm_q = np.empty([dim[0],dim[1],dim[2]])
norm_q = q/n[0,4,:,:] #normalizing it with quiet sun intensity instead of pixel by pixel, local normalization

# norm_q = q/n[0,4,:,:]
#stokv = n[3,:,:,:]
#norm_v = stokv/np.mean(n[0,4,y1:y2,x1:x2])
#
#norm_i = n[0,4,:,:]/np.mean(n[0,4,y1:y2,x1:x2])

# plotting Q_max continuum image to chose a box to find the std/noise level, same box is chosen above for I conti to calculate mean to normalize q
#fig = plt.figure(figsize=(12,12))
#ax = plt.axes()
#im = ax.imshow(norm_q[4,:,:],cmap='gray',vmin=-0.01,vmax=0.01)
#fig.colorbar(im)
#fig.tight_layout(pad=1.8)
#plt.gca().invert_yaxis()
# plt.show()


# y1,y2,x1,x2, = map(int,input("Choose y1,y2,x1,x2: ").split(','))

print("Chosen values of pixels for normalization of q with I conti and also for calculating noise level from q conti")
print(y1,y2,x1,x2)

#summing up normalized q over different spectral positions
combq = norm_q[0,:,:]+norm_q[1,:,:]+norm_q[2,:,:]+norm_q[3,:,:]
#combv = norm_v[0,:,:] + norm_v[1,:,:] - norm_v[2,:,:] - norm_v[3,:,:]


#estimating noise level for q
#std = np.std(norm_q[4,y1:y2,x1:x2])
std = np.std(combq[y1:y2,x1:x2])
meanie = np.mean(norm_q[4,y1:y2,x1:x2])
rms = std/meanie
print("rms for q is {}".format(rms))
print("std for q is {}".format(std))

#estimating noise level for v
#stdv = np.std(norm_v[4,y1:y2,x1:x2])
#meaniev = np.mean(norm_v[4,y1:y2,x1:x2])
#rmsv = stdv/meaniev
#print("rms for v is {}".format(rmsv))
#print("std for v is {}".format(stdv))




fig = plt.figure(figsize=(12,12))
ax = plt.axes()
im = ax.imshow(combq,cmap='gray',vmin=-0.01,vmax=0.01)
fig.colorbar(im)
fig.tight_layout(pad=1.8)
plt.gca().invert_yaxis()
plt.suptitle('Combined LP map')



#fig1 = plt.figure(figsize=(12,12))
#ax1 = plt.axes()
#im1 = ax1.imshow(combv,cmap='gray',vmin=-0.08,vmax=0.08)
#fig1.colorbar(im1)
#fig1.tight_layout(pad=1.8)
#plt.gca().invert_yaxis()
#plt.suptitle('Combined V map')


#
#fig1i = plt.figure(figsize=(12,12))
#ax1i = plt.axes()
#im1i = ax1i.imshow(norm_i,cmap='gray',vmin=0.5,vmax=1.5)
#fig1i.colorbar(im1i)
#fig1i.tight_layout(pad=1.8)
#plt.gca().invert_yaxis()
#plt.suptitle('Normalized Continuum intensity map')


 
#to read in segmented data from Johann's scripts
# seg = fits.open("/home/prabhu/sunrise_holly/imax_maxlp_combi_wv_seg_/imax_maxlp_combi_wv_seg_21.fits")
# a = seg[0].data

#massking W to retain values only above and below 3*sigma
#x = np.ma.masked_inside(combq/4,-1.5*std,1.5*std)
x = np.ma.masked_inside(combq,-2*std,2*std)
#x = np.ma.masked_inside(np.amax(np.abs(norm_q[0:-1,:,:]),axis=0),0,2*std)
#x=np.ma.masked_inside(np.sum(np.abs(norm_q[0:-1,:,:])>2*std,axis=0),0,1)*1.
print(x.shape)
np.ma.set_fill_value(x,np.nan)
filledx = x.filled()
y = np.logical_not(np.ma.masked_inside(combq,-2*std,2*std).mask).astype(np.int)
# print(np.abs(filledx)[158,180:190])
# print(y.astype(int)[158,180:190])
# print(x[158,180:190])

#xv = np.ma.masked_inside(combq,-1.5*stdv,1.5*stdv)
#np.ma.set_fill_value(xv,np.nan)
#filledxv = xv.filled()
#yv = np.ma.masked_greater_equal(np.abs(filledxv),0.02).mask


levels = 1

#plotting masked contours of Q over V
fig3 = plt.figure(figsize=(12,12))
ax3 = plt.axes()
im3 = plt.imshow(x,vmax=0.01,vmin=-0.01)
#im3 = plt.imshow(combv,vmax=0.08,vmin=-0.08,cmap='gray')
#im3_1 = plt.contour(y.astype(int),levels,origin='lower',colors = 'r')#,cmap=plt.get_cmap('jet'))
fig3.colorbar(im3)
fig3.tight_layout(pad=1.8)
plt.gca().invert_yaxis()
plt.suptitle('Contours of LP over combV')

le1 = erode(y,structure=np.ones((2,2))).astype(x.dtype)
le = dilate(le1,structure=[[True,True],[True,True]]).astype(le1.dtype)

fig3 = plt.figure(figsize=(12,12))
ax3 = plt.axes()
im3 = plt.imshow(le,vmax=0.01,vmin=-0.01)
#im3 = plt.imshow(combv,vmax=0.08,vmin=-0.08,cmap='gray')
#im3_1 = plt.contour(y.astype(int),levels,origin='lower',colors = 'r')#,cmap=plt.get_cmap('jet'))
fig3.colorbar(im3)
fig3.tight_layout(pad=1.8)
plt.gca().invert_yaxis()
plt.suptitle('Eroded and dilated')




##plotting masked contours of V over Q
#fig4 = plt.figure(figsize=(12,12))
#ax4 = plt.axes()
## im4 = plt.imshow(xv,vmax=0.08,vmin=-0.08)
#im4 = plt.imshow(combq,vmax=0.01,vmin=-0.01,cmap='gray')
#im4_1 = plt.contour(yv.astype(int),levels,origin='lower',colors = 'r')#,cmap=plt.get_cmap('jet'))
#fig4.colorbar(im4)
#fig4.tight_layout(pad=1.8)
#plt.gca().invert_yaxis()
#plt.suptitle('Contours of V over comb LP ')


#
#fig5 = plt.figure(figsize=(12,12))
#ax5 = plt.axes()
## im4 = plt.imshow(xv,vmax=0.08,vmin=-0.08)
#im5 = plt.imshow(norm_i,vmax=1.5,vmin=0.5,cmap='gray')
#im5_1 = plt.contour(y.astype(int),levels,origin='lower',colors = 'r')#,cmap=plt.get_cmap('jet'))
#fig5.colorbar(im5)
#fig5.tight_layout(pad=1.8)
#plt.gca().invert_yaxis()
#plt.suptitle('Contours of Q over I conti')


#fig6 = plt.figure(figsize=(12,12))
#ax6 = plt.axes()
## im4 = plt.imshow(xv,vmax=0.08,vmin=-0.08)
#im6 = plt.imshow(norm_i,vmax=1.5,vmin=0.5,cmap='gray')
#im6_1 = plt.contour(yv.astype(int),levels,origin='lower',colors = 'r')#,cmap=plt.get_cmap('jet'))
#fig6.colorbar(im6)
#fig6.tight_layout(pad=1.8)
#plt.gca().invert_yaxis()
#plt.suptitle('Contours of V over I conti')
#plt.show()





indx,indq = [[] for i in range(2)]

for j in range(dim[1]):
	for i in range(dim[2]):
		dum = norm_q[:,j,i]
		dumdum = filledx[j,i]
		temp = np.dot(dum[0],dum[1:-1])
		np.abs(dum)
		if (temp[temp<0].size) and (not(np.isnan(dumdum))):
			indq.append((j,i))
		if not(np.isnan(dumdum)):
			indx.append((j,i))
print('Number of pixels above 3sigma in LP with mixed sign')
print(len(indq))
# print(indq[10000:10010])
print('Total number of pixels above 3sigma in LP')
print(len(indx))
# print(w[158,219])

#plotting Q_max profile at different wavelength points, to see if all negative or positive
# fig4 = plt.figure(figsize=(12,12))
# ax4 = plt.axes()
# im4 = plt.plot(range(5),norm_q[:,158,183],'x-')
# fig4.colorbar(im4)
# fig4.tight_layout(pad=1.8)
# plt.gca().invert_yaxis()
# plt.show()
