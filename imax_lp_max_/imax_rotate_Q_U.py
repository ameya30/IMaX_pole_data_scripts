import glob
import numpy               as np
import matplotlib.pyplot   as plt
import matplotlib.gridspec as gridspec

from scipy.optimize import minimize_scalar
from astropy.io     import fits

def Q_eq(phi):
    # Set inputs for fnction
    sumU = 0 
    for wv in range(0, 5):

        uSing = uMes[wv]
        qSing = qMes[wv]
     
        uNew = -1 * qSing * np.sin(phi) + uSing * np.cos(phi)
        sumU += np.square(uNew)

    return sumU

plot_figs = 0

# Load in array

input_list = glob.glob('../Data/mean_rem_output_*.fits')

for i in range (0, len(input_list)):
    fullArr = fits.open(input_list[i])
    fullArr = fullArr[0].data
   
    save_name = '../Data/imax_lp_max_' + input_list[i].split('_')[3]
    save_angs = '../Data/imax_roat_angle_Q_U_' + input_list[i].split('_')[3] 

    print save_name
    print save_angs

    fullDim = np.shape(fullArr)
    
    angArr = np.empty(shape = (fullDim[2], fullDim[3]))
    uNew   = np.empty(shape = (fullDim[1], fullDim[2], fullDim[3]))
    qNew   = np.empty(shape = (fullDim[1], fullDim[2], fullDim[3]))
    
    for x in range (0, fullDim[3]):
        for y in range (0, fullDim[2]):
    
            qMes = fullArr[1, :, y, x] 
            uMes = fullArr[2, :, y, x] 
    
            res = minimize_scalar(Q_eq, bounds=(0, np.pi), method='bounded')
            angle = res['x']
            angArr[y, x] = angle
            
            for wv in range (0, 5):
                uNew[wv, y, x] = -1 * fullArr[1, wv, y, x] * np.sin(angle) + fullArr[2, wv, y, x] * np.cos(angle) 
                qNew[wv, y, x] =      fullArr[1, wv, y, x] * np.cos(angle) + fullArr[2, wv, y, x] * np.sin(angle)

    hdu_ang = fits.PrimaryHDU(angArr)
    hdu_max = fits.PrimaryHDU(qNew)

    hdu_ang.writeto(save_angs)
    hdu_max.writeto(save_name)


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


