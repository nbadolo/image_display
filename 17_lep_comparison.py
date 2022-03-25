#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 14:41:41 2022

@author: nbadolo
"""

import numpy as np
from astropy.io import fits
from scipy import optimize
from astropy.nddata import Cutout2D
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.pyplot import Figure, subplot
#%%
fdir='/home/nbadolo/Bureau/Aymard/Donnees_sph/'
fdir_star=fdir+'star/17_lep_center/'
fdir_psf=fdir+'psf/psf_17_lep_center/'
fname1='zpl_p23_make_polar_maps-ZPL_SCIENCE_P23_REDUCED'
fname2='-zpl_science_p23_REDUCED'
file_I_star_c = fdir_star +fname1+'_I'+fname2+'_I.fits'
file_PI_star_c = fdir_star+fname1+'_PI'+fname2+'_PI.fits'
file_DOLP_star_c = fdir_star+fname1+'_DOLP'+fname2+'_DOLP.fits'
file_AOLP_star_c = fdir_star+ fname1+'_AOLP'+fname2+'_AOLP.fits'

# file_I_psf= fdir_psf+ fname1+'_I'+fname2+'_I.fits'
# file_PI_psf= fdir_psf+fname1+'_PI'+fname2+'_PI.fits'
# file_DOLP_psf= fdir_psf+ fname1+'_DOLP'+fname2+'_DOLP.fits'
# file_AOLP_psf= fdir_psf+fname1+'_AOLP'+fname2+'_AOLP.fits'
#%%
file_lst_cent = [file_I_star_c,file_PI_star_c,file_DOLP_star_c,file_AOLP_star_c]
nFrames_c = len(file_lst_cent)
#%%
fdir='/home/nbadolo//Bureau/Aymard/Donnees_sph/'
fdir_star=fdir+'star/17_lep_bs/'
fdir_psf=fdir+'psf/psf_17_lep_bs/'
fname1='zpl_p23_make_polar_maps-ZPL_SCIENCE_P23_REDUCED'
fname2='-zpl_science_p23_REDUCED'
file_I_star_bs = fdir_star +fname1 +'_I'+ fname2 +'_I.fits'
file_PI_star_bs = fdir_star+fname1 +'_PI'+ fname2 +'_PI.fits'
file_DOLP_star_bs = fdir_star+fname1 +'_DOLP'+ fname2 +'_DOLP.fits'
file_AOLP_star_bs = fdir_star+ fname1 +'_AOLP'+ fname2 +'_AOLP.fits'
#%%
file_lst_bs = [file_I_star_bs,file_PI_star_bs,file_DOLP_star_bs,file_AOLP_star_bs]
nFrames_bs = len(file_lst_bs)
#%%
nDim = 1024
nSubDim = 200 # plage de pixels que l'on veut afficher
size = (nSubDim, nSubDim)
position = (nDim//2, nDim//2)
nDimfigj = [9,10,11]
nDimfigk = [0,1,2]
vmin0 = 3.5
vmax0 = 15
pix2mas = 6.8  #en mas/pix
x_min = -pix2mas*nSubDim//2
x_max = pix2mas*(nSubDim//2-1)
y_min = -pix2mas*nSubDim//2
y_max = pix2mas*(nSubDim//2-1)
#%%
#mean_sub_v_arr=np.empty((nFrames_c,nSubDim//2-1))
sub_cent_arr = np.empty((nFrames_c,nSubDim,nSubDim))
sub_bs_arr = np.empty((nFrames_bs,nSubDim,nSubDim))
im_name_lst =['17_lep I','17_lep PI','17_lep DOLP','17_lep AOLP',
              'HD19887 I','HD19887 PI','HD19887 DOLP','HD19887 AOLP']
im_name_lst_all = ['17_lep_c I','17_lep_c PI','17_lep_c DOLP','17_lep_c AOLP',
                   '17_lep_bs I','17_lep_bs PI','17_lep_bs DOLP','17_lep_bs AOLP',]
Vmin_c = np.empty((nFrames_c))
Vmax_c = np.empty((nFrames_c))

Vmin_bs = np.empty((nFrames_bs))
Vmax_bs = np.empty((nFrames_bs))
#%%
for i in range(nFrames_c):
    hdu = fits.open(file_lst_cent[i])   
    data= hdu[0].data   
    i_cent = data[0,:,:]
   
    cutout = Cutout2D(i_cent, position = position, size = size)
    zoom_hdu = hdu.copy()
    sub_cent = cutout.data
        
    sub_cent_arr[i] = sub_cent
    
    if i == 3:
        Vmin_c[i] = np.min(sub_cent_arr[i])
        Vmax_c[i] = np.max(sub_cent_arr[i])  
    else:
        Vmin_c[i] = np.min(np.log10(sub_cent_arr[i]))
        Vmax_c[i] = np.max(np.log10(sub_cent_arr[i]))

for j in range(nFrames_bs):
    hdu = fits.open(file_lst_cent[j])   
    data = hdu[0].data   
    i_bs = data[0,:,:]
   
    cutout = Cutout2D(i_bs, position = position, size = size)
    zoom_hdu = hdu.copy()
    sub_bs = cutout.data
        
    sub_bs_arr[j] = sub_bs
    
    if j == 3:
        Vmin_bs[j] = np.min(sub_bs_arr[j])
        Vmax_bs[j] = np.max(sub_bs_arr[j])  
    else:
        Vmin_bs[j] = np.min(np.log10(sub_bs_arr[j]))
        Vmax_bs[j] = np.max(np.log10(sub_bs_arr[j]))
#%%
plt.figure('I_PRIM (les deux modes comparés)')
plt.clf()
for k in range (nFrames_c + nFrames_bs):   
    plt.subplot(2,4,k+1)
    if 0 <= k< nFrames_c:
        if k == 3:
            plt.imshow(sub_cent_arr[k], cmap ='inferno', origin ='lower',
            vmin = Vmin_c[k], vmax = Vmax_c[k], extent = [x_min , x_max, y_min , y_max])
            plt.colorbar(label = 'ADU')
            plt.clim(Vmin_c[k],Vmax_c[k])
        else:
            plt.imshow(np.log10(sub_cent_arr[k]), cmap = 'inferno', origin ='lower',
            vmin = Vmin_c[k], vmax = Vmax_c[k], extent = [x_min , x_max, y_min , y_max])
            plt.colorbar(label ='')
            plt.clim(Vmin_c[k],Vmax_c[k])
    else:
        if k == 7:
            plt.imshow(sub_bs_arr[k-4], cmap = 'inferno', origin='lower',
            vmin = Vmin_bs[k-4], vmax=Vmax_bs[k-4], extent = [x_min , x_max, y_min , y_max])
            plt.colorbar(label='ADU')
            plt.clim(Vmin_bs[k-4],Vmax_bs[k-4])
        else:
            plt.imshow(np.log10(sub_bs_arr[k-4]), cmap ='inferno', origin ='lower',
            vmin = Vmin_bs[k-4], vmax = Vmax_bs[k-4], extent = [x_min , x_max, y_min , y_max]) 
            plt.colorbar(label='')
            plt.clim(Vmin_bs[k-4],Vmax_bs[k-4])
    if k == 6 or k == 7:
        plt.text(-1.1*pix2mas*size[0]//6., 2*pix2mas*size[1]//6., im_name_lst_all[k], color='w',
              fontsize='large', ha='center')
    else:
        plt.text(-1.3*pix2mas*size[0]//6., 2*pix2mas*size[1]//6., im_name_lst_all[k], color='w',
              fontsize='large', ha='center')
    #plt.colorbar(label='ADU in log$_{10}$ scale')
    # plt.clim(Vmin[i],Vmax[i])
# plt.xlabel('Relative R.A.(mas)', size=10)   
# plt.ylabel('Relative Dec.(mas)', size=10)

    if k != 1 and k != 2 and k != 3 :
        if k ==4:
            plt.xlabel('Relative R.A.(mas)', size=10) 
            plt.ylabel('Relative Dec.(mas)', size=10)
        else:
            if k == 0:
                plt.ylabel('Relative Dec.(mas)', size=10)
            else:
                plt.xlabel('Relative R.A.(mas)', size=10) 
    # else:
    #     plt.ylabel('Relative Dec.(mas)', size=10)    
#%%