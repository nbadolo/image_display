#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 14:28:00 2022

@author: nbadolo
"""



#packages
import numpy as np
from matplotlib import pyplot as plt
from math import pi, cos, sin
from astropy.nddata import Cutout2D
from astropy.io import fits
import scipy.optimize as opt
from skimage import io, color, measure, draw, img_as_bool
import pylab as plt
#import matplotlib.pyplot as plt


#%%
#file
fdir='/home/nbadolo//Bureau/Aymard/Donnees_sph/'
fdir_star=fdir+'log/SW_Col/star/both/V_N_R/'
fdir_psf=fdir+'psf/HD204971_mira/'
fname1='zpl_p23_make_polar_maps-ZPL_SCIENCE_P23_REDUCED'
fname2='-zpl_science_p23_REDUCED'
file_I_star= fdir_star + fname1 +'_I'+fname2+'_I.fits'
file_PI_star= fdir_star+fname1+'_PI'+fname2+'_PI.fits'
file_DOLP_star= fdir_star+fname1+'_DOLP'+fname2+'_DOLP.fits'
file_AOLP_star= fdir_star+ fname1+'_AOLP'+fname2+'_AOLP.fits'

file_I_psf= fdir_psf+ fname1+'_I'+fname2+'_I.fits'
file_PI_psf= fdir_psf+fname1+'_PI'+fname2+'_PI.fits'
file_DOLP_psf= fdir_psf+ fname1+'_DOLP'+fname2+'_DOLP.fits'
file_AOLP_psf= fdir_psf+fname1+'_AOLP'+fname2+'_AOLP.fits'
#%%
#creating lists
file_lst=[file_I_star,file_PI_star,file_DOLP_star,file_AOLP_star]
          #,file_I_psf,file_PI_psf,file_DOLP_psf,file_AOLP_psf]
nFrames=len(file_lst)
#%%
#parameters
nDim=1024
nSubDim = 100 # plage de pixels que l'on veut afficher
size = (nSubDim, nSubDim)
nDimfigj=[9,10,11]
nDimfigk=[0,1,2]
vmin0=3.5
vmax0=15
pix2mas = 6.8  #en mas/pix
x_min=-pix2mas*nSubDim//2
x_max=pix2mas*(nSubDim//2-1)
y_min=-pix2mas*nSubDim//2
y_max=pix2mas*(nSubDim//2-1)
position = (nDim//2,nDim//2)
size = (nSubDim, nSubDim)

x, y = np.meshgrid(np.arange(nSubDim),np.arange(nSubDim)) #cree un tableau 
sub_v_arr=np.empty((nFrames,nSubDim,nSubDim))
ellips_arr = np.empty((nFrames,nSubDim,nSubDim))
Vmin = np.empty((nFrames))
Vmax = np.empty((nFrames))
#%%
#opening file
for i in range (nFrames):
      hdu = fits.open(file_lst[i])   
      data = hdu[0].data
      intensity = data[1,:,:]
      
      cutout = Cutout2D(intensity, position=position, size=size)
      zoom_hdu = hdu.copy()
      sub_v = cutout.data
            
      Ellips = np.zeros_like(sub_v) #creation d'un tableau de meme forme que sub_v
      Ellips[sub_v>0.05*np.max(sub_v)] = 1 # on retient les points d'intensité 
      ellips_arr[i] = Ellips                             #égale à 5% de Imax et à tous ces points 
                                                         #on donne 1 comme valeur d'intensité
      Vmin[i] = np.min(ellips_arr[i])
      Vmax[i] = np.max(ellips_arr[i])  
      
      

image = ellips_arr[1]
#%%
plt.figure(7)
plt.clf()
plt.imshow(image, cmap ='inferno', vmin=Vmin[1], vmax=Vmax[1], origin='lower')
plt.title('Intensity at 5%')

#%%
#image = img_as_bool(io.imread('bubble.jpg')[..., 0])
regions = measure.regionprops(measure.label(image))
bubble = regions[0]
#%%
y0, x0 = bubble.centroid
r = bubble.major_axis_length / 2.
#%%
def cost(params):
    x0, y0, r = params
    coords = draw.disk((y0, x0), r, shape=image.shape)
    template = np.zeros_like(image)
    template[coords] = 1
    return -np.sum(template == image)

x0, y0, r = opt.fmin(cost, (x0, y0, r))
#%%
f, ax = plt.subplots()
circle = plt.Circle((x0, y0), r)
ax.imshow(image, cmap='inferno', interpolation='nearest', origin='lower')
ax.add_artist(circle)
plt.show()