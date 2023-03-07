#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 17:51:34 2023

@author: nbadolo
"""

"""
Code simple pour tester rapidement la qualité d'une carte de couleur'
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


star_name = 'AH-SCO'
fdir= '/home/nbadolo/Bureau/Aymard/Donnees_sph/'+star_name+ '/'

fname1='zpl_p23_make_polar_maps-ZPL_SCIENCE_P23_REDUCED'
fname2='-zpl_science_p23_REDUCED'
file_I_star= fdir + fname1+'_I'+fname2+'_I.fits'
file_PI_star= fdir +fname1+'_PI'+fname2+'_PI.fits'
file_DOLP_star= fdir +fname1+'_DOLP'+fname2+'_DOLP.fits'
file_AOLP_star= fdir + fname1+'_AOLP'+fname2+'_AOLP.fits'

file_lst = [file_I_star,file_PI_star, file_DOLP_star, file_AOLP_star]
nFrames = len(file_lst)

nDim = 1024
nSubDim = 200 # plage de pixels que l'on veut afficher
size = (nSubDim, nSubDim)
# nDimfigj = [3, 4, 5]
# nDimfigk = [6, 7, 8]




pix2mas = 3.4  # en mas/pix
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


sub_v_arr = np.empty((nFrames,nSubDim,nSubDim))

im_name_lst = ['I','PI','DOLP', 'AOLP']
Vmin = np.empty((nFrames))
Vmax = np.empty((nFrames))

position = (nDim//2,nDim//2)
size = (nSubDim, nSubDim)

x, y = np.meshgrid(np.arange(nSubDim), np.arange(nSubDim)) #cree un tableau 

R = np.sqrt((x-nSubDim/2)**2+(y-nSubDim/2)**2)
r = np.linspace(1,nSubDim//2-1,nSubDim//2-1) # creation d'un tableau de distance radiale

r_mas=pix2mas*r #  où r est en pixels et r_mas en millièmes d'arcseconde

for i in range (nFrames):
      hdu = fits.open(file_lst[i])[0]   
      data = hdu.data   
      i_v = data[0,:,:] 
                        
      cutout = Cutout2D(i_v, position = position, size=size)
      zoom_hdu = hdu.copy()
      sub_v = cutout.data

      
      sub_v_arr[i] = sub_v
      Vmin[i] = np.min(sub_v_arr[i])
      Vmax[i] = np.max(sub_v_arr[i])
      
        
      fig = plt.figure()
      plt.clf()
      fig.set_size_inches(18.5, 10, forward = True)
      for i in range (nFrames):
            plt.subplot(2,2,i+1)
            plt.imshow(sub_v_arr[i], cmap='inferno', origin='lower',
            vmin=Vmin[i], vmax=Vmax[i], extent = [x_min , x_max, y_min , y_max])
            
            plt.text(size[0]//10, 2*pix2mas*size[1]//6.,
                      f'{star_name}' + '_' + f'{im_name_lst[i]}', color='w',
                  fontsize='large', ha='center')
            plt.colorbar(label='ADU in log$_{10}$ scale', shrink = 0.6)