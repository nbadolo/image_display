#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 10:08:35 2022

@author: nbadolo
"""


#packages
import numpy as np
from matplotlib import pyplot as plt
from math import pi, cos, sin
from astropy.nddata import Cutout2D
from astropy.io import fits
import scipy.optimize as opt
import pylab as plt

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
      sub_v_arr[i] = sub_v
      Ellips = np.zeros_like(sub_v_arr[i]) #creation d'un tableau de meme forme que sub_v
      Ellips[sub_v_arr[i]> 0.05*np.max(sub_v_arr[i])] = sub_v_arr[i][sub_v_arr[i] > 0.05*np.max(sub_v_arr[i])]# on retient les points d'intensité égale à 5% de Imax 
      ellips_arr[i] = Ellips
      
      #if np.any(np.min(ellips_arr[i]) != 'nan'): #i==3 or i==7:
      Vmin[i] = np.min(ellips_arr[i])
      Vmax[i] = np.max(ellips_arr[i])  
      # else:
      #    Vmin[i] = np.min(np.log10(ellips_arr[i]))
      #    Vmax[i] = np.max(np.log10(ellips_arr[i]))
#%%
plt.figure(0)
plt.clf()
#if np.any(np.min(ellips_arr[i]) != 'nan'): #i==3 or i==7:
plt.imshow(ellips_arr[1], cmap ='inferno', vmin=Vmin[1], vmax=Vmax[1], origin='lower')
plt.title('Intensity at 5%')

# else:
#     plt.imshow(np.log10(ellips_arr[1]), cmap='inferno', origin='lower',
#     vmin=Vmin[1], vmax = Vmax[1], )
#%%
#ellipse function 
"""
## ellipse function variables
u= x-position of the center
v= y-position of the center
a= radius on the x-axis
b= radius on the y-axis
theta =  #rotation angle
""" 

def ellips(t, u, v, a, b, theta):

    
    Ell = np.array([a*np.cos(t) , b*np.sin(t)])  
         #u,v removed to keep the same center location
    M_rot = np.array([[cos(theta) , -sin(theta)],[sin(theta) , cos(theta)]])  
         #2-D rotation matrix

    Ell_rot = np.zeros((nSubDim, Ell.shape[1]))
    for i in range(Ell.shape[1]):
        Ell_rot[:,i] = np.dot(M_rot,Ell[:,i])
        Ell_rot[0,i] = Ell_rot[0,i] + u
        Ell_rot[1,i] = Ell_rot[1,i] + v
    return Ell_rot.ravel() # .ravel permet de passer de deux dimension à une seule
#%%
t = np.linspace(0, 2*pi, nSubDim)
data_intensity = ellips_arr[1] - np.min(ellips_arr[1]) + 0.1 # ajout d'un offset
data_intensity = data_intensity.ravel()
initial_guess = (1, 0.5, 2, 1.5, 0) # variable initiales à ajuster avec les moindre carré
popt, pcov = opt.curve_fit(ellips, t, data_intensity, p0=initial_guess)# ajustement avec les moindre carres
#%%
data_fitted = ellips(t, *popt)
#%%
#plt.plot( u+Ell[0,:] , v+Ell[1,:] )     #initial ellipse
plt.imshow(ellips_arr[1], cmap ='inferno', vmin=Vmin[1], vmax=Vmax[1], origin='lower')
#plt.plot( u + Ell_rot[0,:] , v + Ell_rot[1,:],'darkorange' )  #rotated ellipse
plt.plot( data_fitted[0,:] , data_fitted[1,:],'darkorange' ) #rotated fit
plt.grid(color='lightgray',linestyle='--')
plt.show()
