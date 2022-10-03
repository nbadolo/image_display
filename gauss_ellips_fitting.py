#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 17:23:56 2022

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
#fdir_psf=fdir+'psf/HD204971_mira/'
fname1='zpl_p23_make_polar_maps-ZPL_SCIENCE_P23_REDUCED'
fname2='-zpl_science_p23_REDUCED'
file_I_star= fdir_star + fname1 +'_I'+fname2+'_I.fits'
file_PI_star= fdir_star+fname1+'_PI'+fname2+'_PI.fits'
file_DOLP_star= fdir_star+fname1+'_DOLP'+fname2+'_DOLP.fits'
file_AOLP_star= fdir_star+ fname1+'_AOLP'+fname2+'_AOLP.fits'

# file_I_psf= fdir_psf+ fname1+'_I'+fname2+'_I.fits'
# file_PI_psf= fdir_psf+fname1+'_PI'+fname2+'_PI.fits'
# file_DOLP_psf= fdir_psf+ fname1+'_DOLP'+fname2+'_DOLP.fits'
# file_AOLP_psf= fdir_psf+fname1+'_AOLP'+fname2+'_AOLP.fits'
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
threshold = 0.05
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
sub_v_arr = np.empty((nFrames,nSubDim,nSubDim))
ellips_arr = np.empty((nFrames,nSubDim,nSubDim))
ellips_im_arr = np.empty((nFrames,nSubDim,nSubDim))
Vmin = np.empty((nFrames))
Vmax = np.empty((nFrames))
Vmin_r = np.empty((nFrames))
Vmax_r = np.empty((nFrames)) 
Vmin_w = np.empty((nFrames))
Vmax_w = np.empty((nFrames)) 
#%%
#opening file
for i in range (nFrames):
      hdu = fits.open(file_lst[i])   
      data = hdu[0].data
      intensity = data[1,:,:]
      
      cutout = Cutout2D(intensity, position=position, size=size)
      zoom_hdu = hdu.copy()
      sub_v = cutout.data
      
      sub_v_arr[i] = sub_v
      
      Ellips_im = np.zeros_like(sub_v_arr[i]) #creation d'un tableau de meme forme que sub_v
      Ellips_im[sub_v > threshold*np.max(sub_v)] = sub_v[sub_v > threshold*np.max(sub_v)]# on retient les points d'intensité égale à 5% de Imax 
      ellips_im_arr[i] = Ellips_im
      
      Ellips = np.zeros_like(sub_v)          #creation d'un tableau de meme forme que sub_v
      Ellips[sub_v> threshold*np.max(sub_v)] = 1   # on retient les points d'intensité 
      ellips_arr[i] = Ellips                 # égale à 5% de Imax et à tous ces points 
                                             # on donne 1 comme valeur d'intensité

        
      Vmin = np.min(sub_v_arr[i])
      Vmax = np.max(sub_v_arr[i])
      
      Vmin_r[i] = np.min(ellips_im_arr[i])
      Vmax_r[i] = np.max(ellips_im_arr[i])  
      
      Vmin_w[i] = np.min(ellips_arr[i])
      Vmax_w[i] = np.max(ellips_arr[i])  
#%%      
      
im_white = ellips_arr[1]
im_real = ellips_im_arr[1]

#%%
# plot of the white image at 5% *Imax
plt.figure('white ellips')
plt.clf()
plt.imshow(im_white, cmap ='inferno', vmin=Vmin_w[1], vmax=Vmax_w[1], origin='lower')
plt.title('Intensity at 5%')
#%%
# plot of the real image at 5% *Imax
plt.figure('real image')
plt.clf()
plt.imshow(im_real, cmap ='inferno', vmin=Vmin_r[1], vmax=Vmax_r[1], origin='lower')
plt.title('Intensity at 5%')

#%%
#image = img_as_bool(io.imread('bubble.jpg')[..., 0])
regions = measure.regionprops(measure.label(im_white))
bubble = regions[0]
#%%
y_i, x_i = bubble.centroid
a_i = bubble.major_axis_length / 2.
b_i = 0.75*bubble.major_axis_length / 2.
theta_i  = pi/4
t = np.linspace(0, 2*pi, nSubDim)
#%%
def cost(params):
    x0, y0, a, b, theta = params
    #coords = draw.disk((y0, x0), r, shape=image.shape)
    coords = draw.ellipse(y0, x0, a, b, shape=None, rotation= theta)
    template = np.zeros_like(im_white)
    template[coords] = 1
    return -np.sum(template == im_white)

x_f, y_f, a_f, b_f, theta_f = opt.fmin(cost, (x_i, y_i, a_i, b_i, theta_i))
#%%
#def ellips(t, x_f, y_f, a_f, bb_f, theta_f):

    
Ell = np.array([a_f*np.cos(t) , b_f*np.sin(t)])  
     #u,v removed to keep the same center location
M_rot = np.array([[cos(theta_f) , -sin(theta_f)],[sin(theta_f) , cos(theta_f)]])  
     #2-D rotation matrix

Ell_rot = np.zeros((2, Ell.shape[1]))
for i in range(Ell.shape[1]):
    Ell_rot[:,i] = np.dot(M_rot,Ell[:,i])
    Ell_rot[0,i] = Ell_rot[0,i] + x_f
    Ell_rot[1,i] = Ell_rot[1,i] + y_f
#return Ell_rot.ravel() # .ravel permet de passer de deux dimension à une seule
#%%
plt.figure('white ellips contour at 5%')
plt.clf()
plt.imshow(ellips_arr[1], cmap ='inferno', vmin=Vmin_w[1], vmax=Vmax_w[1], origin='lower')
#plt.plot( u + Ell_rot[0,:] , v + Ell_rot[1,:],'darkorange' )  #rotated ellipse
plt.plot( Ell_rot[0,:] , Ell_rot[1,:],'darkorange' ) #rotated fit
#plt.grid(color='lightgray',linestyle='--')
plt.show()
#%%
plt.figure('real image contour at 5%')
plt.clf()
plt.imshow(ellips_im_arr[1], cmap ='inferno', vmin = Vmin_r[1], vmax = Vmax_r[1], origin='lower')
#plt.plot( u + Ell_rot[0,:] , v + Ell_rot[1,:],'darkorange' )  #rotated ellipse
plt.plot( Ell_rot[0,:] , Ell_rot[1,:],'darkorange' ) #rotated fit
#plt.grid(color='lightgray',linestyle='--')
plt.show()
#%%
plt.figure('full image and contour at 5%')
plt.clf()
plt.imshow(sub_v_arr[1], cmap ='inferno', vmin=Vmin_w[1], vmax=Vmax_r[1], origin='lower')
#plt.plot( u + Ell_rot[0,:] , v + Ell_rot[1,:],'darkorange' )  #rotated ellipse
plt.plot( Ell_rot[0,:] , Ell_rot[1,:],'darkorange' ) #rotated fit
#plt.grid(color='lightgray',linestyle='--')
plt.show()


