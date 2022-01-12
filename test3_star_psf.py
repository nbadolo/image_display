#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 15:15:09 2021

@author: nbadolo
"""

"""
Code simplifié 
"""

import numpy as np
from astropy.io import fits
from scipy import optimize
from astropy.nddata import Cutout2D
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.pyplot import Figure, subplot
#%%
fdir='/home/nbadolo//Bureau/Aymard/Donnees_sph/'
fdir_star=fdir+'star/'
fdir_psf=fdir+'psf/'
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
#%%
file_lst=[file_I_star,file_PI_star,file_DOLP_star,file_AOLP_star,
          file_I_psf,file_PI_psf,file_DOLP_psf,file_AOLP_psf]
nFrames=len(file_lst)
#%%
nDim=1024
nSubDim = 200 # plage de pixels que l'on veut afficher
size = (200, 200)
nDimfigj=[9,10,11]
nDimfigk=[0,1,2]
vmin0=3.5
vmax0=15
x_min=-3.5*nSubDim//2
x_max=3.5*(nSubDim//2-1)
y_min=-3.5*nSubDim//2
y_max=3.5*(nSubDim//2-1)
#%%
mean_sub_v_arr=np.empty((nFrames,nSubDim//2-1))
sub_v_arr=np.empty((nFrames,nSubDim,nSubDim))
im_name_lst =['Mira I','Mira PI','Mira DOLP','Mira AOP',
              'HD204971 I','HD204971 PI','HD204971 DOLP','HD204971 AOLP']
Vmin=np.empty((nFrames))
Vmax=np.empty((nFrames))
#%%
position = (nDim//2,nDim//2)
size = (nSubDim, nSubDim)

x, y = np.meshgrid(np.arange(nSubDim),np.arange(nSubDim)) #cree un tableau 

R = np.sqrt((x-nSubDim/2)**2+(y-nSubDim/2)**2)
r = np.linspace(1,nSubDim//2-1,nSubDim//2-1)
pix2mas=3.5  #en mas/pix
r_mas=pix2mas*r #  où r est en pixels et r_mas en millièmes d'arcseconde


for i in range(nFrames):
    hdu = fits.open(file_lst[i])   
    data= hdu[0].data   
    i_v= data[0,:,:]
   
    cutout = Cutout2D(i_v, position=position, size=size)
    zoom_hdu = hdu.copy()
    sub_v = cutout.data
    
    f = lambda r : sub_v[(R >= r-0.5) & (R < r+0.5)].mean()   
    mean_sub_v = np.vectorize(f)(r) 
    
    mean_sub_v_arr[i]=mean_sub_v 
    sub_v_arr[i]=sub_v
    if i==3 or i==7:
        Vmin[i]=np.min(sub_v_arr[i])
        Vmax[i]=np.max(sub_v_arr[i])  
    else:
        Vmin[i]=np.min(np.log10(sub_v_arr[i]))
        Vmax[i]=np.max(np.log10(sub_v_arr[i]))  
      
#%%
plt.figure(0)
plt.clf()
for i in range (nFrames):   
    plt.subplot(3,4,i+1)
    if i==3 or i==7:
        plt.imshow(sub_v_arr[i], cmap='inferno', origin='lower',vmin=Vmin[i], vmax=Vmax[i], extent = [x_min , x_max, y_min , y_max])
    else:
        plt.imshow(np.log10(sub_v_arr[i]), cmap='inferno', origin='lower',vmin=Vmin[i], vmax=Vmax[i], extent = [x_min , x_max, y_min , y_max])   
    plt.text(size[0]//3., size[1]-size[1]/8., im_name_lst[i], color='cyan',
             fontsize='x-small', ha='center')
    plt.colorbar(label='ADU in log$_{10}$ scale')
    plt.clim(Vmin[i],Vmax[i])
plt.xlabel('Relative R.A.(mas)', size=10)   
    # plt.ylabel('Relative Dec.(mas)', size=10)
   

for j in range(len(nDimfigj)):      
    plt.subplot(3,4,nDimfigj[j])
    plt.plot(r_mas, np.log10(mean_sub_v_arr[j]), color='darkorange',linewidth=2, label='Mira') 
    plt.plot(r_mas, np.log10(mean_sub_v_arr[j+4]),color='purple',linewidth=2, label='HD204971') 
    plt.legend(loc=0) 
    plt.xlabel('r (mas)', size=10) 
    if j==0:
        plt.ylabel(r'Intensity in log$_{10}$ scale', size=10)
