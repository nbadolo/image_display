#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 15:43:46 2022

@author: nbadolo
"""


"""
Code simplifié pour comparer les profiles radiaux de s-lep; 17_lep et w_pup
et profile radial d'intensité'
"""

import numpy as np
import os
import scipy 
from os.path import exists
from astropy.io import fits
from scipy import optimize
from astropy.nddata import Cutout2D
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.pyplot import Figure, subplot
import webbrowser
#%% 
# =============================================================================
# 1:17_Lep
# 2: S_Lep
# 3: w_pup
# =============================================================================

        
fdir1= '/home/nbadolo/Bureau/Aymard/Donnees_sph/log/17_Lep/star/alone/I_PRIM/'
fdir2= '/home/nbadolo/Bureau/Aymard/Donnees_sph/log/S_Lep/star/alone/I_PRIM/'
fdir3= '/home/nbadolo/Bureau/Aymard/Donnees_sph/w_Pup/star/alone/I_PRIM/'


    
fname1='zpl_p23_make_polar_maps-ZPL_SCIENCE_P23_REDUCED'
fname2='-zpl_science_p23_REDUCED'

#ouverture de 17_Lep
file_I_star1= fdir1+ fname1+'_I'+fname2+'_I.fits'
file_PI_star1= fdir1 +fname1+'_PI'+fname2+'_PI.fits'
file_DOLP_star1= fdir1 +fname1+'_DOLP'+fname2+'_DOLP.fits'
file_AOLP_star1= fdir1 + fname1+'_AOLP'+fname2+'_AOLP.fits'

#ouverture de s_lep
file_I_star2= fdir2 + fname1+'_I'+fname2+'_I.fits'
file_PI_star2= fdir2 +fname1+'_PI'+fname2+'_PI.fits'
file_DOLP_star2= fdir2 +fname1+'_DOLP'+fname2+'_DOLP.fits'
file_AOLP_star2= fdir2 + fname1+'_AOLP'+fname2+'_AOLP.fits'

#ouverture w_pup
file_I_star3= fdir3 + fname1+'_I'+fname2+'_I.fits'
file_PI_star3= fdir3 +fname1+'_PI'+fname2+'_PI.fits'
file_DOLP_star3= fdir3 +fname1+'_DOLP'+fname2+'_DOLP.fits'
file_AOLP_star3= fdir3 + fname1+'_AOLP'+fname2+'_AOLP.fits'

  
# file_lst = [file_I_star,file_PI_star,file_DOLP_star,file_AOLP_star,
#           file_I_psf,file_PI_psf,file_DOLP_psf,file_AOLP_psf]

file_lst1 = [file_I_star1, file_PI_star1, file_DOLP_star1, file_AOLP_star1]
file_lst2 = [file_I_star2, file_PI_star2, file_DOLP_star2, file_AOLP_star2]
file_lst3 = [file_I_star2, file_PI_star2, file_DOLP_star2, file_AOLP_star2]          

curve_name1 = ['17_Lep', 'S_Lep', 'w_Pup']
curve_name2 = ['I', 'PI', 'DOLP', 'ALOP']
nFrames1 = len(file_lst1)
nFrames2 = len(file_lst2)
nFrames3 = len(file_lst3)

nDimfigj = [3, 4, 5]
nDimfigk = [1, 2, 3]
 
nDim = 1024
nSubDim = 200 # plage de pixels que l'on veut afficher
size = (nSubDim, nSubDim)
# nDimfigj = [3, 4, 5]
# nDimfigk = [6, 7, 8]
vmin0 = 3.5
vmax0 = 15
pix2mas = 6.8  #en mas/pix
x_min = -pix2mas*nSubDim//2
x_max = pix2mas*(nSubDim//2-1)
y_min = -pix2mas*nSubDim//2
y_max = pix2mas*(nSubDim//2-1)
X, Y= np.meshgrid(np.linspace(-100,99,200), np.linspace(-100,99,200))
X_, Y_= np.meshgrid(np.linspace(-nDim/2,nDim/2-1,nDim), np.linspace(-nDim/2,nDim/2-1,nDim))

X *= pix2mas
Y *= pix2mas
X_ *= pix2mas
Y_ *= pix2mas

X_step = 10
X_step_ = 50

mean_sub_v_arr1 = np.empty((nFrames1,nSubDim//2-1))
mean_sub_v_arr2 = np.empty((nFrames2,nSubDim//2-1))
mean_sub_v_arr3 = np.empty((nFrames3,nSubDim//2-1))
sub_v_arr1 = np.empty((nFrames1,nSubDim,nSubDim))
sub_v_arr2 = np.empty((nFrames2,nSubDim,nSubDim))
sub_v_arr3 = np.empty((nFrames3,nSubDim,nSubDim))
im_name_lst = ['I','PI','DOLP','AOLP',
                'I','PI','DOLP','AOLP']
Vmin1 = np.empty((nFrames1))
Vmax1 = np.empty((nFrames1))

Vmin2 = np.empty((nFrames2))
Vmax2 = np.empty((nFrames2))

Vmin3 = np.empty((nFrames3))
Vmax3 = np.empty((nFrames3))

position = (nDim//2,nDim//2)
size = (nSubDim, nSubDim)

x, y = np.meshgrid(np.arange(nSubDim), np.arange(nSubDim)) #cree un tableau 

R = np.sqrt((x-nSubDim/2)**2+(y-nSubDim/2)**2)
r = np.linspace(1,nSubDim//2-1,nSubDim//2-1) # creation d'un tableau de distance radiale

r_mas = pix2mas*r #  où r est en pixels et r_mas en millièmes d'arcseconde


for i in range (nFrames1):
      hdu = fits.open(file_lst2[i])[0]   
      data1 = hdu.data   
      i_v1 = data1[0,:,:]
      fltr = hdu.header.get('HIERARCH ESO INS3 OPTI5 NAME')     
      #print(fltr)                   
      cutout1 = Cutout2D(i_v1, position=position, size=size)
      zoom_hdu = hdu.copy()
      sub_v1 = cutout1.data
    
      f = lambda r : sub_v1[(R >= r-0.5) & (R < r+0.5)].mean()   
      mean_sub_v = np.vectorize(f)(r) 
    
      mean_sub_v_arr1[i] = mean_sub_v 
      sub_v_arr1[i] = sub_v1
          
      
      Vmin1[i] = np.min(np.log10(sub_v1+np.abs(np.min(sub_v1))+10))
      Vmax1[i] = np.max(np.log10(sub_v1+np.abs(np.min(sub_v1))+10))
      




   
for i in range (nFrames2):
      hdu = fits.open(file_lst2[i])[0]   
      data2 = hdu.data   
      i_v2 = data2[0,:,:]
      fltr = hdu.header.get('HIERARCH ESO INS3 OPTI5 NAME')     
      #print(fltr)                   
      cutout2 = Cutout2D(i_v2, position=position, size=size)
      zoom_hdu = hdu.copy()
      sub_v2 = cutout2.data
    
      f = lambda r : sub_v2[(R >= r-0.5) & (R < r+0.5)].mean()   
      mean_sub_v = np.vectorize(f)(r) 
    
      mean_sub_v_arr2[i] = mean_sub_v 
      sub_v_arr2[i] = sub_v2
      
      Vmin2[i] = np.min(np.log10(sub_v2+np.abs(np.min(sub_v2))+10))
      Vmax2[i] = np.max(np.log10(sub_v2+np.abs(np.min(sub_v2))+10))
      
      
      
      
for i in range (nFrames3):
      hdu3 = fits.open(file_lst3[i])   
      data3 = hdu3[0].data   
      i_v3 = data3[0,:,:]
   
      cutout3 = Cutout2D(i_v3, position=position, size=size)
      zoom_hdu3 = hdu3.copy()
      sub_v3 = cutout3.data
    
      f = lambda r : sub_v3[(R >= r-0.5) & (R < r+0.5)].mean()   
      mean_sub_v3 = np.vectorize(f)(r) 
    
      mean_sub_v_arr3[i] = mean_sub_v3
      sub_v_arr3[i] = sub_v3
      
      Vmin3[i] = np.min(np.log10(sub_v3+np.abs(np.min(sub_v3))+10))
      Vmax3[i] = np.max(np.log10(sub_v3+np.abs(np.min(sub_v3))+10))
      
plt.clf()
fig = plt.figure('Comparison of radial profile')
fig.set_size_inches(18.5, 10, forward = True)

# for k in range(nFrames1):  
#     if k < 3:    
#         plt.subplot(1,3,k+1)
#         #plt.plot(r_mas, np.log10(mean_sub_v_arr1[k]), color='blue',
#         #        linewidth = 2, label= f'{curve_name1[k]}' +'_'  +f'{curve_name2[k]}') 
#         #plt.plot(r_mas, np.log10(mean_sub_v_arr2[k]), color='darkorange',
#         #        linewidth = 2, label= f'{curve_name1[k]}' +'_' + f'{curve_name2[k]}') 
#         plt.plot(r_mas, np.log10(mean_sub_v_arr3[k]),color='purple',
#                 linewidth = 2, label = f'{curve_name1[k]}' + '_' + f'{curve_name2[k]}') 
#         plt.legend(loc=0) 
#         plt.xlabel('r (mas)', size=10) 
       
#         plt.ylabel(r'Intensity in log$_{10}$ scale', size=10)
    
#         plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/w_Pup/plots/comparison.pdf', 
#                 dpi=100, bbox_inches ='tight')

        
plt.subplot(1,3,1)
plt.plot(r_mas, np.log10(mean_sub_v_arr1[0]), color='blue',
        linewidth = 2, label= f'{curve_name1[0]}' +'_'  +f'{curve_name2[0]}') 
plt.plot(r_mas, np.log10(mean_sub_v_arr2[0]), color='darkorange',
        linewidth = 2, label= f'{curve_name1[0]}' +'_' + f'{curve_name2[0]}') 
# plt.plot(r_mas, np.log10(mean_sub_v_arr3[0]),color='purple',
#         linewidth = 2, label = f'{curve_name1[0]}' + '_' + f'{curve_name2[0]}') 
# plt.legend(loc=0) 
plt.xlabel('r (mas)', size=10) 
   
plt.ylabel(r'Intensity in log$_{10}$ scale', size=10)

plt.subplot(1,3,2)
plt.plot(r_mas, np.log10(mean_sub_v_arr1[1]), color='blue',
        linewidth = 2, label= f'{curve_name1[1]}' + '_'  +f'{curve_name2[1]}') 
plt.plot(r_mas, np.log10(mean_sub_v_arr2[1]), color ='darkorange',
        linewidth = 2, label= f'{curve_name1[1]}' + '_' + f'{curve_name2[1]}') 
plt.plot(r_mas, np.log10(mean_sub_v_arr3[1]),color='purple',
        linewidth = 2, label = f'{curve_name1[1]}' + '_' + f'{curve_name2[1]}') 
plt.legend(loc=0) 
plt.xlabel('r (mas)', size=10) 
   
plt.ylabel(r'Intensity in log$_{10}$ scale', size=10)

plt.subplot(1,3,3)
plt.plot(r_mas, np.log10(mean_sub_v_arr1[2]), color='blue',
        linewidth = 2, label= f'{curve_name1[2]}' +'_'  +f'{curve_name2[2]}') 
plt.plot(r_mas, np.log10(mean_sub_v_arr2[2]), color='darkorange',
        linewidth = 2, label= f'{curve_name1[2]}' +'_' + f'{curve_name2[2]}') 
plt.plot(r_mas, np.log10(mean_sub_v_arr3[2]),color='purple',
        linewidth = 2, label = f'{curve_name1[2]}' + '_' + f'{curve_name2[2]}') 
plt.legend(loc=0) 
plt.xlabel('r (mas)', size=10) 
   
plt.ylabel(r'Intensity in log$_{10}$ scale', size=10)

plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/w_Pup/plots/comparison.pdf', 
        dpi=100, bbox_inches ='tight')
        
# for k in range(len(nDimfigk)):      
#       plt.subplot(3,3,nDimfigk[k] + 1)
#       plt.plot(r_mas, np.log10(mean_sub_v_arr2[k]), color='darkorange',
#               linewidth = 2, label= f'{star_name}') 
#       plt.plot(r_mas, np.log10(mean_sub_v_arr3[k]),color='purple',
#               linewidth = 2, label = f'{star_name}' + '_psf') 
#       plt.legend(loc=0) 
#       plt.xlabel('r (mas)', size=10) 
        
        
plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/w_Pup/plots/comparison.png', 
                dpi=100, bbox_inches ='tight')
plt.tight_layout()
 