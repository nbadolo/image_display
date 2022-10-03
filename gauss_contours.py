#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 09:00:55 2022

@author: nbadolo
"""



import numpy as np
import matplotlib.pyplot as plt
from astropy.nddata import Cutout2D
from astropy.io import fits
import scipy.optimize as opt
import pylab as plt

#%%



# file path seting
fdir='/home/nbadolo//Bureau/Aymard/Donnees_sph/'
fdir_star=fdir+'star/Mira/'
fdir_psf=fdir+'psf/HD204971_mira/'
fname1='zpl_p23_make_polar_maps-ZPL_SCIENCE_P23_REDUCED'
fname2='-zpl_science_p23_REDUCED'
file_I_star= fdir_star +fname1+'_I'+fname2+'_I.fits'
file_PI_star= fdir_star+fname1+'_PI'+fname2+'_PI.fits'
file_DOLP_star= fdir_star+fname1+'_DOLP'+fname2+'_DOLP.fits'
file_AOLP_star= fdir_star+ fname1+'_AOLP'+fname2+'_AOLP.fits'

file_I_psf= fdir_psf+ fname1+'_I'+fname2+'_I.fits'
file_PI_psf= fdir_psf+fname1+'_PI'+fname2+'_PI.fits'
file_DOLP_psf= fdir_psf+ fname1+'_DOLP'+fname2+'_DOLP.fits'
file_AOLP_psf= fdir_psf+fname1+'_AOLP'+fname2+'_AOLP.fits'

# creating list
file_lst=[file_I_star,file_PI_star,file_DOLP_star,file_AOLP_star]
          #,file_I_psf,file_PI_psf,file_DOLP_psf,file_AOLP_psf]
nFrames=len(file_lst)

#parameters
nDim=1024
nSubDim = 100 # plage de pixels que l'on veut afficher
size = (nSubDim, nSubDim)
nDimfigj=[9,10,11]
nDimfigk=[0,1,2]
vmin0=3.5
vmax0=15
pix2mas=6.8  #en mas/pix
x_min=-pix2mas*nSubDim//2
x_max=pix2mas*(nSubDim//2-1)
y_min=-pix2mas*nSubDim//2
y_max=pix2mas*(nSubDim//2-1)
position = (nDim//2,nDim//2)
size = (nSubDim, nSubDim)

x, y = np.meshgrid(np.arange(nSubDim),np.arange(nSubDim)) #cree un tableau 
sub_v_arr=np.empty((nFrames,nSubDim,nSubDim))
#%%
for i in range (nFrames):
      hdu = fits.open(file_lst[i])   
      data = hdu[0].data   
      i_v = data[0,:,:]
      
      cutout = Cutout2D(i_v, position=position, size=size)
      zoom_hdu = hdu.copy()
      sub_v = cutout.data
      sub_v_arr[i] = sub_v
# gauss function definition
#%%  
"""
# meaning of twoD_Gaussian parameters
xdata_tuple: points grids for the model plot
amplitude : exterior factor of gaussian funtion. amplitude = q/(2pi*sigma**2) where q is an arbitrary constant
 xo and  yo : values of µx and µy in gaussian formula, the coordinate of the peak
 theta : rotaional angle
 offset : minimum of the funtion values

"""

#%%
def twoD_Gaussian(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    '''twoD_Gaussian(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset)
    '''
    (x,  y) = xdata_tuple
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()
#%%
# Create x and y indices
    # x = np.linspace(0, 200, 201)
    # y = np.linspace(0, 200, 201)
    # x, y = np.meshgrid(x, y)
    
x, y = np.meshgrid(np.arange(nSubDim),np.arange(nSubDim)) #cree un tableau 
#%%

plt.figure(0)
plt.clf()
plt.imshow(np.log10(sub_v_arr[1]))
plt.title('PI')
#%%
print(np.max(sub_v_arr[1]))
#%%
#create data
data_intensity = sub_v_arr[1] - np.min(sub_v_arr[1]) + 0.1 # ajout d'un offset
initial_guess = (315, nSubDim//2, nSubDim//2, 40, 40, 30, 0) # variable à ajuster avec les moindre carré
popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), data_intensity.ravel(), p0=initial_guess)# ajustement avec les moindre carres
#%%
data_fitted = twoD_Gaussian((x, y), *popt)
#%%
line_array2 = np.asarray([0.01, 0.05, 0.5])
line_array = line_array2*popt[0]

fig, ax = plt.subplots(1, 1)
#ax.imshow(data_fitted.reshape(nSubDim, nSubDim), cmap=plt.cm.jet, origin='lower',
#   extent=(x.min(), x.max(), y.min(), y.max()))

ax.imshow(np.log10(data_intensity.reshape(nSubDim, nSubDim)), cmap=plt.cm.inferno, origin = 'lower',
   extent=(x.min(), x.max(), y.min(), y.max()))
CS = ax.contour(x, y, data_fitted.reshape(nSubDim, nSubDim), line_array, colors='b') # fait les contours à une hauteur de 8

fmt = {}
strs = ['1%', '5%', '50%']
for l, s in zip(CS.levels, strs):
    fmt[l] = s

ax.clabel(CS, CS.levels, inline=True, fontsize=10, fmt = fmt)
plt.show()

