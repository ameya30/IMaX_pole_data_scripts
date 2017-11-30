import os

import numpy               as np 
import matplotlib.pyplot   as plt

from   pandas                import *
from   pandas.tools.plotting import table
from   astropy.io            import fits
from   scipy                 import optimize

#################### FUNCTIONS #################### 

# Find distance of each point from centre

def to_cen(yc, xc):

    return np.sqrt((x_ind - xc) ** 2 + (y_ind - yc) ** 2)

# Find distance from point to mean circle 

def to_cir(c):

    Ri = to_cen(*c)

    return Ri - Ri.mean()

###################### BEGIN ######################

# Get data

res = fits.open('/home/prabhu/sunrise_holly/work_IMaX_pole_scripts/imax_find_sun/imax_find_sun_input.fits')
#this text file was a result of find_limb.py (inflection point stuff), saving the file was done via spyder and not through the script find_limb, dont get confused
ycor = np.loadtxt('/home/prabhu/sunrise_holly/work_IMaX_pole_scripts/my_scripts/inflection_coord.txt')
#this was then used in place of line 52,53 coords from the contour, 80 was subtracted since we consider a smaller box, where origin is moved to 80,80(matrix y,x)
y_ind = ycor[31:900]-80
xcor = np.arange(936)
x_ind = xcor[31:900]-80
res = res[0].data

# Plot limb image

plt.imshow(res[0, 4, 80:936-730, 80:936-70])

# Find limb, plot contours

limb = plt.contour(res[0, 4, 80:936-730, 80:936-80], 3, linewidths = 2)

path = limb.collections[0].get_paths()[0]

points = path.vertices

#x_ind = points[:, 0]
#y_ind = points[:, 1]

# Find circle

pixScale = 50. / 936

radSunAproxPix = 960 / pixScale

cenEst = (y_ind[0] + radSunAproxPix,x_ind[0] - radSunAproxPix)

cenEst2, ier = optimize.leastsq(to_cir, cenEst)

yCen, xCen = cenEst2

Ri_2 = to_cen(xCen, yCen)

R_2 = Ri_2.mean()

residu_2 = sum((Ri_2 - R_2) ** 2)

residu2_2 = sum((Ri_2 ** 2 - R_2 ** 2) ** 2)

# Plot circle

#with contours as inputs ######################comment 33,35 uncomment 50,51
#theta_fit = np.linspace(-np.pi/2+np.pi/50.025,-np.pi/2+np.pi/27.9, 1800)
#x_fit = xCen + (R_2) * np.cos(theta_fit)
#y_fit = (yCen-330) + (R_2) * np.sin(theta_fit)

#with inflection as input, the values were tweaked to fit the best circle, the final,total x_fit,yfit,theta_fit and R_2 were then
#fed into erode_dilate, on lines 90-94
theta_fit = np.linspace(-np.pi/2+np.pi/50.025,-np.pi/2+np.pi/27, 1800)
x_fit = xCen + (R_2) * np.cos(theta_fit)
y_fit = (yCen-325.7) + (R_2) * np.sin(theta_fit)

plt.plot(x_fit, y_fit, 'k-')

# Create results table

headings = ['Radius', 'xCen', 'yCen']

dataTable = [[round(R_2), round(xCen) + 730, round(yCen) + 80]]

bestEstTable = DataFrame(columns = headings, data = dataTable)

tmpHead = ['BED']

tmpTable = [['Best estimate data']]

tmpDf = DataFrame(columns = tmpHead, data = tmpTable)

bestEstTable.index = tmpDf['BED']

# Set up figure

fig, axes = plt.subplots(figsize = (6, 0.5))

axes.xaxis.set_visible(False)
axes.yaxis.set_visible(False)

axes.set_frame_on(False)

table(axes, bestEstTable, loc = 'upper right')

# Save fig

# plt.savefig('../Figures/imax_find_sun.png', bbox_inches = 'tight', dpi = 500)
