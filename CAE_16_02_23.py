#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 17:01:07 2023

@author: nbadolo
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 23 14:29:41 2022

@author: nbadolo
"""

"""
Code simplifié pour l'affichage simultané de tous les alone et both  des étoiles 
sans psf:  flux d'intensité
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

def log_image(star_name, obsmod):
#%%       
    fdir= '/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+ '/'
    fdir_star = fdir + 'star/'+obsmod+ '/' 
    fdir_psf = fdir +'psf/'+obsmod+ '/'
    lst_fltr_star1 = os.listdir(fdir_star)
    #print(lst_fltr_star1)
    n_lst_fltr_star1 = len(lst_fltr_star1)
    #print(n_lst_fltr_star1)
    lst_fltr_star2 = []
    nDimfigj = [3, 4, 5]
    nDimfigk = [6, 7, 8]
    for p in range(n_lst_fltr_star1):
        fdir_fltr_data_star = fdir_star + lst_fltr_star1[p]
        lst_fltr_data_star = os.listdir(fdir_fltr_data_star) 
        n_lst_fltr_data_star = len(lst_fltr_data_star)
        if n_lst_fltr_data_star != 0:
            lst_fltr_star2.append(lst_fltr_star1[p])
    n_lst_fltr_star2 = len(lst_fltr_star2)
    #print(lst_fltr_star2)
    
    
    for l in range(n_lst_fltr_star2):
       
        fdir_star_fltr = fdir_star + lst_fltr_star2[l] +'/'
        fdir_psf_fltr = fdir_psf + lst_fltr_star2[l] + '/'
                
        fname1='zpl_p23_make_polar_maps-ZPL_SCIENCE_P23_REDUCED'
        fname2='-zpl_science_p23_REDUCED'
        file_I_star= fdir_star_fltr + fname1+'_I'+fname2+'_I.fits'
        file_PI_star= fdir_star_fltr +fname1+'_PI'+fname2+'_PI.fits'
        file_DOLP_star= fdir_star_fltr +fname1+'_DOLP'+fname2+'_DOLP.fits'
        file_AOLP_star= fdir_star_fltr + fname1+'_AOLP'+fname2+'_AOLP.fits'
        
        file_lst2 = [file_I_star, file_PI_star, file_DOLP_star, file_AOLP_star]
        nFrames2 = len(file_lst2)
        nDim = 1024
        nSubDim = 60 # plage de pixels que l'on veut afficher
        size = (nSubDim, nSubDim)
         # nDimfigj = [3, 4, 5]
         # nDimfigk = [6, 7, 8]
    
        pix2mas = 3.4  # en mas/pix
        x_min = -pix2mas*nSubDim//2
        x_max = pix2mas*(nSubDim//2-1)
        y_min = -pix2mas*nSubDim//2
        y_max = pix2mas*(nSubDim//2-1)
        X, Y= np.meshgrid(np.linspace(-nSubDim/2,nSubDim/2-1,nSubDim), np.linspace(-nSubDim/2,  nSubDim/2-1,nSubDim))
        X_, Y_= np.meshgrid(np.linspace(-nDim/2,nDim/2-1,nDim), np.linspace(-nDim/2,nDim/2-1,nDim))
        
        X *= pix2mas
        Y *= pix2mas
        X_ *= pix2mas
        Y_ *= pix2mas
        
        X_step = 10
        X_step_ = 50
        
         #mean_sub_v_arr = np.empty((nFrames,nSubDim//2-1))
        mean_sub_v_arr2 = np.empty((nFrames2,nSubDim//2-1))
         #mean_sub_v_arr3 = np.empty((nFrames3,nSubDim//2-1))
         #sub_v_arr = np.empty((nFrames,nSubDim,nSubDim))
        sub_v_arr2 = np.empty((nFrames2,nSubDim,nSubDim))
         #sub_v_arr3 = np.empty((nFrames3,nSubDim,nSubDim))
        im_name_lst = ['I','PI','DOLP', 'AOLP']
        Vmin = np.empty((nFrames2))
        Vmax = np.empty((nFrames2))
        
        position = (nDim//2,nDim//2)
        size = (nSubDim, nSubDim)
        
        x, y = np.meshgrid(np.arange(nSubDim), np.arange(nSubDim)) #cree un tableau 
        
        R = np.sqrt((x-nSubDim/2)**2+(y-nSubDim/2)**2)
        r = np.linspace(1,nSubDim//2-1,nSubDim//2-1) # creation d'un tableau de distance radiale
        
        r_mas = pix2mas*r #  où r est en pixels et r_mas en millièmes d'arcseconde
        
        for i in range (nFrames2):
              hdu = fits.open(file_lst2[i])[0]   
              data2 = hdu.data   
              i_v = data2[0,:,:]
              fltr = hdu.header.get('HIERARCH ESO INS3 OPTI5 NAME')     
              #print(fltr)                   
              cutout2 = Cutout2D(i_v, position=position, size=size)
              zoom_hdu = hdu.copy()
              sub_v = cutout2.data
            
              f = lambda r : sub_v[(R >= r-0.5) & (R < r+0.5)].mean()   
              mean_sub_v = np.vectorize(f)(r) 
            
              mean_sub_v_arr2[i] = mean_sub_v 
              sub_v_arr2[i] = sub_v
              
              Vmin[i] = np.min(sub_v+np.abs(np.min(sub_v))+10)
              Vmax[i] = np.max(sub_v+np.abs(np.min(sub_v))+10)
         
        U2 = np.cos(np.pi*sub_v_arr2[3]/180)
        V2 = np.sin(np.pi*sub_v_arr2[3]/180)
              
        plt.clf()
        fig = plt.figure(f'{star_name}' +'(' + f'{fltr}' + '_Cam1' + '_I )')
        fig.set_size_inches(12, 6, forward = True)
        plt.imshow(np.log10(sub_v_arr2[0] + np.abs(np.min(sub_v_arr2[0]))), cmap='inferno', origin='lower',
                        vmin= np.log10(Vmin[0]), vmax= np.log10(Vmax[0]), extent = [x_min , x_max, y_min , y_max])
                       
        plt.text(size[0]//10, 2.3*pix2mas*size[1]//6.,
                                  f'{star_name}' + '_' + f'{im_name_lst[0]}', color='w',
                              fontsize='large', ha='center')
        plt.colorbar(label='ADU in log$_{10}$ scale', shrink = 0.6)
        plt.xlabel('Relative R.A.(mas)', size=14)
        plt.ylabel('Relative Dec.(mas)', size=14)
                  
        plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+
                          '/plots/CAE/'+ star_name +'_' + fltr + '_Cam1' + '_I.pdf', 
                          dpi=100, bbox_inches ='tight')
        
        
        plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+
                          '/plots/CAE/'+ star_name +'_' + fltr + '_Cam1' + '_I.png', 
                          dpi=100, bbox_inches ='tight')
        plt.tight_layout()
        
 #%%       
        plt.clf()
        fig = plt.figure(f'{star_name}' +'(' + f'{fltr}' + '_Cam1' + '_PI )')
        fig.set_size_inches(12, 6, forward = True)          
        plt.imshow(sub_v_arr2[1], cmap='inferno', origin='lower',
                      vmin=Vmin[1], vmax=Vmax[1], extent = [x_min , x_max, y_min , y_max])
                      
        plt.text(size[0]//10, 2.3*pix2mas*size[1]//6.,
                                f'{star_name}' + '_' + f'{im_name_lst[1]}', color='w',
                            fontsize='large', ha='center')
        plt.colorbar(label='ADU', shrink = 0.6)

        plt.xlabel('Relative R.A.(mas)', size=14)
        plt.ylabel('Relative Dec.(mas)', size=14)
                  
        plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+
                          '/plots/CAE/'+ star_name +'_' + fltr + '_Cam1' + '_PI.pdf', 
                          dpi=100, bbox_inches ='tight')
        
        
        plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+
                          '/plots/CAE/'+ star_name +'_' + fltr + '_Cam1' + '_PI.png', 
                          dpi=100, bbox_inches ='tight')
        plt.tight_layout()      
#%%     
        plt.clf()
        fig = plt.figure(f'{star_name}' +'(' + f'{fltr}' + '_Cam1' + '_PI Polvect )')
        fig.set_size_inches(12, 6, forward = True)
        plt.imshow(sub_v_arr2[1], cmap ='inferno', origin='lower',vmin=Vmin[1], 
                    vmax=Vmax[1], extent = [x_min , x_max, y_min , y_max])   
        plt.colorbar(label='ADU in log$_{10}$ scale', shrink = 0.6)       
        #q = plt.quiver(X[::X_step,::X_step],Y[::X_step,::X_step],U2[::X_step,::X_step], V2[::X_step,::X_step])
        q = plt.quiver(X[::X_step,::X_step],Y[::X_step,::X_step],U2[::X_step,::X_step], V2[::X_step,::X_step], color='w', pivot='mid', width=0.004, headwidth=0.8, headlength=1.5)
                 # headlength=1e-10, headaxislength=1e-10)
        #plt.quiverkey(q, X = 0.1, Y = 1.03, U = 0, label='', labelpos='E')                 
        plt.text(size[0]//10, 2.3*pix2mas*size[1]//6.,
               f'{star_name}' + f'{im_name_lst[1]}'+ '_&_Pol. vect', color='c',
                    fontsize='large', ha='center')
        
        plt.xlabel('Relative R.A.(mas)', size=14)
        plt.ylabel('Relative Dec.(mas)', size=14)
        plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+
                          '/plots/CAE/'+ star_name +'_' + fltr + '_Cam1' + '_PI_vect.pdf', 
                          dpi=100, bbox_inches ='tight')
        
        
        plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+
                          '/plots/CAE/'+ star_name +'_' + fltr + '_Cam1' + '_PI_vect.png', 
                          dpi=100, bbox_inches ='tight')
        plt.tight_layout()  
         
        msg1='reduction okay for '+ star_name+'_Cam1'
         #return(msg1)
        print(msg1)
#log_image('17_Lep', 'alone')

star_name = '17_Lep'
obsmod ='alone'
log_image(star_name, obsmod)