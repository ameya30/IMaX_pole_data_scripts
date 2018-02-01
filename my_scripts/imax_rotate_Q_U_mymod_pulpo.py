import glob
from natsort import natsorted
import numpy               as np
import matplotlib.pyplot   as plt
#import matplotlib.gridspec as gridspec

from scipy.optimize import minimize_scalar
from astropy.io     import fits
#from scipy.io import readsav

def Q_eq(phi):
    # Set inputs for fnction
#    sumU = 0 
    uNew = qMes*np.sin(2*phi) + uMes*np.cos(2*phi)
   # qNew = qMes*np.cos(2*phi) + uMes*-1*np.sin(2*phi)
    
   
    sumU = np.sum(np.square(uNew))
#    for wv in range(0, 5):
#
#        uSing = uMes[wv]
#        qSing = qMes[wv]
#     	
#        uNew = -1 * qSing * np.sin(2*phi) + uSing * np.cos(2*phi)
#        sumU += np.square(uNew)
    return sumU
def rotate(ang):
    uNew1 = qMes*np.sin(2*ang) + uMes*np.cos(2*ang)
    qNew1 = qMes*np.cos(2*ang) + uMes*-1*np.sin(2*ang)
    
    return qNew1,uNew1
plot_figs = 0

# Load in array

#input_list = natsorted(glob.glob('/scratch/prabhu/backup_workstation/sunrise_holly/binned_cycles/*.fits'))
#input_list = ['/scratch/prabhu/backup_workstation/sunrise_holly/movie_data/mk_magneto_tr2_reduc_rnr_300_21.sav','/scratch/prabhu/backup_workstation/sunrise_holly/movie_data/mk_magneto_tr2_reduc_rnr_300_22.sav']
input_list = natsorted(glob.glob('/scratch/prabhu/backup_workstation/sunrise_holly/mean_rem_output/mean_rem_output_nr_21.fits'))

             
for i in range (0, len(input_list)):
    fullArr = fits.open(input_list[i])
    fullArr = fullArr[0].data #0 for restored data 1 for non restored #############
#    fullArr = readsav(input_list[i],python_dict=True)['iid'] #iid for restored, iidn for nonrestored
#    save_name = '/scratch/prabhu/backup_workstation/sunrise_holly/imax_lp_max_/imax_lp_max_nr_' + ("_").join(input_list[i].split('_')[-2::])
#    save_angs = '/scratch/prabhu/backup_workstation/sunrise_holly/imax_lp_max_/imax_roat_angle_Q_U_nr_' + ("_").join(input_list[i].split('_')[-2::]) 
#    save_name = '/scratch/prabhu/backup_workstation/sunrise_holly/imax_lp_max_/imax_lp_max_' + input_list[i].split('_')[-1].split('.')[0]+".fits"
#    save_angs = '/scratch/prabhu/backup_workstation/sunrise_holly/imax_lp_max_/imax_roat_angle_Q_U_' + input_list[i].split('_')[-1].split('.')[0]+".fits"
    if len(input_list[i].split('_'))==8:
        save_name = '/scratch/prabhu/backup_workstation/sunrise_holly/imax_lp_max_/tr_imax_lp_max_norm_' + input_list[i].split('_')[-1]
        save_angs = '/scratch/prabhu/backup_workstation/sunrise_holly/imax_lp_max_/tr_imax_roat_angle_Q_U_norm_' + input_list[i].split('_')[-1]
    else:
        save_name = '/scratch/prabhu/backup_workstation/sunrise_holly/imax_lp_max_/tr_imax_lp_max_norm_' + '_'.join(input_list[i].split('_')[-2::])
        save_angs = '/scratch/prabhu/backup_workstation/sunrise_holly/imax_lp_max_/tr_imax_roat_angle_Q_U_norm_' + '_'.join(input_list[i].split('_')[-2::])
    print (save_name)
    print (save_angs)

    fullDim = np.shape(fullArr)
    
    angArr = np.empty(shape = (fullDim[2], fullDim[3]))
    uNew   = np.empty(shape = (fullDim[1], fullDim[2], fullDim[3]))
    qNew   = np.empty(shape = (fullDim[1], fullDim[2], fullDim[3]))
    
    for x in range (0, fullDim[3]):
        print (x)
        for y in range (0, fullDim[2]):
    	    
            qMes = fullArr[1, :, y, x] 
            uMes = fullArr[2, :, y, x] 
            res = minimize_scalar(Q_eq, bounds=(0, np.pi), method='bounded')
            angle = res['x']
            angArr[y, x] = angle
#            
#            for wv in range (0, 5):
#                uNew[wv, y, x] = -1 * fullArr[1, wv, y, x] * np.sin(angle) + fullArr[2, wv, y, x] * np.cos(angle) 
#                qNew[wv, y, x] =      fullArr[1, wv, y, x] * np.cos(angle) + fullArr[2, wv, y, x] * np.sin(angle)
            qNew[:,y,x],uNew[:, y, x] =  rotate(angle)  
#saving the new q and u in the same file as primary[0] and image[1] hdu's
    hdu_ang = fits.PrimaryHDU(angArr)
    hdu_max = fits.PrimaryHDU(qNew)
    hdu_u = fits.ImageHDU(uNew,name='u')
    hdulist = fits.HDUList([hdu_max,hdu_u])
    hdu_ang.writeto(save_angs)
    hdulist.writeto(save_name)


if plot_figs == 1:
    fig, axes = plt.subplots(ncols   = 3,
                             nrows   = 1, 
                             figsize = (12, 18))
    
    fig.subplots_adjust(left   = 0.04,
                        right  = 0.97, 
                        top    = 0.99,
                        bottom = 0.05,
                        wspace = 0.08,
                        hspace = 0.15)
    
    imAng = axes[0].imshow(angArr[:, :],  cmap = 'gray', clim = (-0.05,0.05))
    minpl = axes[1].imshow(qNew[1, :, :], cmap = 'gray', clim = (-0.05,0.05))
    maxpl = axes[2].imshow(uNew[1, :, :], cmap = 'gray', clim = (-0.05,0.05))
   
    
    for ax in axes:
        ax.invert_yaxis()
    
    plt.savefig('../Figures/imax_angle_arr.png', dpi = 500)
    
    fig, axes = plt.subplots(ncols   = 2,
                             nrows   = 2, 
                             figsize = (12, 6))
    
    axes[1,1].plot(uNew[:, 102, 636])
    axes[0,0].plot(fullArr[1,:, 102, 636])
    axes[0,1].plot(fullArr[2,:, 102, 636])
    axes[1,0].plot(qNew[:, 102, 636])
    axes[0,0].set_title('Mesured Q')
    axes[0,1].set_title('Mesured U')
    axes[1,0].set_title('new Q' )
    axes[1,1].set_title('new U')
    plt.savefig('../Figures/imax_roation_line_change.png', dpi = 500)


