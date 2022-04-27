#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 12:58:11 2022

@author: nbadolo
"""

"""
Code simplifié pour l'affichage simultané de tous les alone  et de sa psf: flux 
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
def log_image(star_name, obsmod):   
    

#%%        
    fdir= '/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+ '/'
    fdir_star = fdir + 'star/'+obsmod+ '/' 
    fdir_psf = fdir +'psf/'+obsmod+ '/'
    lst_fltr_star = os.listdir(fdir_star)
    n_lst_fltr_star = len(lst_fltr_star)
    lst_fltr2_star = []
    nDimfigj = [3, 4, 5]
    nDimfigk = [6, 7, 8]
    for p in range(n_lst_fltr_star):
        
        fdir_fltr_data_star = fdir_star + lst_fltr_star[p]
        lst_fltr_data_star = os.listdir(fdir_fltr_data_star) 
        n_lst_fltr_data_star = len(lst_fltr_data_star)
        if n_lst_fltr_data_star != 0:
            lst_fltr2_star.append(lst_fltr_star[p])
    print(lst_fltr2_star)
    
    
    lst_fltr_psf = os.listdir(fdir_psf)
    n_lst_fltr_psf = len(lst_fltr_psf)
    lst_fltr2_psf = []
    for n in range(n_lst_fltr_psf):
        
        fdir_fltr_data_psf = fdir_psf + lst_fltr_psf[n]
        lst_fltr_data_psf = os.listdir(fdir_fltr_data_psf) 
        n_lst_fltr_data_psf = len(lst_fltr_data_psf)
        if n_lst_fltr_data_psf != 0:
            lst_fltr2_psf.append(lst_fltr_psf[n])
    print(lst_fltr2_psf)
    
    lst_fltr3 = list(set(lst_fltr2_star).intersection(lst_fltr2_psf))
    print(lst_fltr3)
    n_lst_fltr3 = len(lst_fltr3)
    for j in range(n_lst_fltr3):
        fdir_star_fltr = fdir_star + lst_fltr3[j] +'/'
        fdir_psf_fltr = fdir_psf + lst_fltr3[j] + '/'
        
        fname1='zpl_p23_make_polar_maps-ZPL_SCIENCE_P23_REDUCED'
        fname2='-zpl_science_p23_REDUCED'
        file_I_star= fdir_star_fltr + fname1+'_I'+fname2+'_I.fits'
        file_PI_star= fdir_star_fltr +fname1+'_PI'+fname2+'_PI.fits'
        file_DOLP_star= fdir_star_fltr +fname1+'_DOLP'+fname2+'_DOLP.fits'
        file_AOLP_star= fdir_star_fltr + fname1+'_AOLP'+fname2+'_AOLP.fits'
    
        file_I_psf = fdir_psf_fltr + fname1+'_I'+fname2+'_I.fits'
        file_PI_psf = fdir_psf_fltr +fname1+'_PI'+fname2+'_PI.fits'
        file_DOLP_psf = fdir_psf_fltr + fname1+'_DOLP'+fname2+'_DOLP.fits'
        file_AOLP_psf = fdir_psf_fltr +fname1+'_AOLP'+fname2+'_AOLP.fits'
      
        file_lst = [file_I_star,file_PI_star,file_DOLP_star,file_AOLP_star,
                  file_I_psf,file_PI_psf,file_DOLP_psf,file_AOLP_psf]
        
        file_lst2 = [file_I_star, file_PI_star, file_DOLP_star, file_AOLP_star]
        file_lst3 = [file_I_psf, file_PI_psf, file_DOLP_psf, file_AOLP_psf]          
        
        nFrames = len(file_lst)
        nFrames2 = len(file_lst2)
        nFrames3 = len(file_lst3)
    
   
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
        
        #mean_sub_v_arr = np.empty((nFrames,nSubDim//2-1))
        mean_sub_v_arr2 = np.empty((nFrames2,nSubDim//2-1))
        mean_sub_v_arr3 = np.empty((nFrames3,nSubDim//2-1))
        #sub_v_arr = np.empty((nFrames,nSubDim,nSubDim))
        sub_v_arr2 = np.empty((nFrames2,nSubDim,nSubDim))
        sub_v_arr3 = np.empty((nFrames3,nSubDim,nSubDim))
        im_name_lst = ['I','PI','DOLP','AOLP',
                        'I','PI','DOLP','AOLP']
        Vmin2 = np.empty((nFrames3))
        Vmax2 = np.empty((nFrames3))
        
        Vmin3 = np.empty((nFrames3))
        Vmax3 = np.empty((nFrames3))
    
        position = (nDim//2,nDim//2)
        size = (nSubDim, nSubDim)
        
        x, y = np.meshgrid(np.arange(nSubDim), np.arange(nSubDim)) #cree un tableau 
        
        R = np.sqrt((x-nSubDim/2)**2+(y-nSubDim/2)**2)
        r = np.linspace(1,nSubDim//2-1,nSubDim//2-1)
        
        r_mas=pix2mas*r #  où r est en pixels et r_mas en millièmes d'arcseconde
    
      # """
      # Filtre utilisé: I_PRIM 
      # """
    
        for i in range (nFrames2):
              hdu = fits.open(file_lst2[i])   
              data2 = hdu[0].data   
              i_v2 = data2[0,:,:]
           
              cutout2 = Cutout2D(i_v2, position=position, size=size)
              zoom_hdu = hdu.copy()
              sub_v2 = cutout2.data
            
              f = lambda r : sub_v2[(R >= r-0.5) & (R < r+0.5)].mean()   
              mean_sub_v = np.vectorize(f)(r) 
            
              mean_sub_v_arr2[i] = mean_sub_v 
              sub_v_arr2[i] = sub_v2
              if np.any(np.min(sub_v_arr2[i])<= 0): #i==3 or i==7:
                  Vmin2[i] = np.min(sub_v_arr2[i])
                  Vmax2[i] = np.max(sub_v_arr2[i])  
              else:
                  Vmin2[i] = np.min(np.log10(sub_v_arr2[i]))
                  Vmax2[i] = np.max(np.log10(sub_v_arr2[i]))  
         
              U2 = sub_v_arr2[2]*np.cos(np.pi*sub_v_arr2[3]/180)
              V2 = sub_v_arr2[2]*np.sin(np.pi*sub_v_arr2[3]/180)
              
              
              
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
              if np.any(np.min(sub_v_arr3[i])<= 0): #i==3 or i==7:
                  Vmin3[i] = np.min(sub_v_arr3[i])
                  Vmax3[i] = np.max(sub_v_arr3[i])  
              else:
                  Vmin3[i] = np.min(np.log10(sub_v_arr3[i]))
                  Vmax3[i] = np.max(np.log10(sub_v_arr3[i]))  
         
              U3 = sub_v_arr3[2]*np.cos(np.pi*sub_v_arr3[3]/180)
              V3 = sub_v_arr3[2]*np.sin(np.pi*sub_v_arr3[3]/180)
              
              shap = np.shape(sub_v_arr2)
              print(shap)
        mean_sub_v_arr =  mean_sub_v_arr2 + mean_sub_v_arr3   
        sub_v_arr = sub_v_arr2 + sub_v_arr3     
        Vmin = Vmin2 + Vmin3
        Vmax = Vmax2 + Vmax3
        
        plt.figure(f'{star_name}' +'(' + f'{ lst_fltr3[j]}' + ')', figsize=(40,20.5))
        plt.clf()    
        for i in range (nFrames2):
              plt.subplot(3,3,i+1)
              if i != 2 and i!= 3:
                  if np.any(np.min(sub_v_arr2[i])<= 0):           
                      plt.imshow(sub_v_arr2[i], cmap='inferno', origin='lower',
                      vmin=Vmin2[i], vmax=Vmax2[i], extent = [x_min , x_max, y_min , y_max])
                      
                      plt.text(size[0]//10, 2*pix2mas*size[1]//6.,
                                f'{star_name}' + '_' + f'{im_name_lst[i]}', color='w',
                            fontsize='large', ha='center')
                      plt.colorbar(label='ADU in log$_{10}$ scale')
                  else:
                       plt.imshow(np.log10(sub_v_arr2[i]), cmap='inferno', origin='lower',
                       vmin=Vmin2[i], vmax=Vmax2[i], extent = [x_min , x_max, y_min , y_max])
                       
                       plt.text(size[0]//10, 2*pix2mas*size[1]//6.,
                                 f'{star_name}' + '_' + f'{im_name_lst[i]}', color='w',
                             fontsize='large', ha='center')
                       plt.colorbar(label='ADU in log$_{10}$ scale')
              else :
                  if i == 2:                  
                      if np.any(np.min(sub_v_arr2[1])<= 0):
                          plt.imshow(sub_v_arr2[1], cmap ='inferno', origin='lower',vmin=Vmin2[1], 
                                        vmax=Vmax2[1], extent = [x_min , x_max, y_min , y_max])   
                          plt.colorbar(label='ADU in log$_{10}$ scale')       
                          q = plt.quiver(X[::X_step,::X_step],Y[::X_step,::X_step],U2[::X_step,::X_step], V2[::X_step,::X_step])
                          plt.quiverkey(q, X = -0.05, Y = 1.03, U = 0.03, label='pol. degree vector norm scale 0.03 ', labelpos='E')
                          # plt.text(-17*size[0]//25., 3*size[1]//2, im_name_lst[3], color='cyan',
                          #            fontsize='x-small', ha='center')
                          
                          plt.text(size[0]//10, 2*pix2mas*size[1]//6.,
                                    f'{star_name}' + '_' + f'{im_name_lst[i]}', color='w',
                                fontsize='large', ha='center')
                          #plt.colorbar(label='ADU in log$_{10}$ scale')
                    
                      else :
                          plt.imshow(np.log10(sub_v_arr2[1]), cmap ='inferno', origin='lower',vmin=Vmin2[1], 
                                         vmax=Vmax2[1], extent = [x_min , x_max, y_min , y_max])   
                          plt.colorbar(label='ADU in log$_{10}$ scale')       
                          q = plt.quiver(X[::X_step,::X_step],Y[::X_step,::X_step],U2[::X_step,::X_step], V2[::X_step,::X_step])
                          plt.quiverkey(q, X = -0.05, Y = 1.03, U = 0.03, label='pol. degree vector norm scale 0.03 ', labelpos='E')
                          # plt.text(-17*size[0]//25., 3*size[1]//2, im_name_lst[3], color='cyan',
                          #         fontsize='x-small', ha='center') 
                          
                          
                          plt.text(size[0]//10, 2*pix2mas*size[1]//6.,
                                    f'{star_name}' + '_' + f'{im_name_lst[i]}'+ 'Pol. vect', color='w',
                                fontsize='large', ha='center')
                         
        for j in range (len(nDimfigj)):   
            plt.subplot(3,3,(nDimfigj[j] + 1))
            if j == 2:
                print(nDimfigj[j])
                if np.any(np.min(sub_v_arr3[1])<= 0):
                    plt.imshow(sub_v_arr3[1], cmap ='inferno', origin='lower',vmin=Vmin3[1], 
                                        vmax=Vmax3[1], extent = [x_min , x_max, y_min , y_max])   
                    plt.colorbar(label='ADU in log$_{10}$ scale')       
                    q_ = plt.quiver(X[::X_step,::X_step],Y[::X_step,::X_step],U3[::X_step,::X_step], V3[::X_step,::X_step])
                    plt.quiverkey(q_, X = -0.05, Y = 1.03, U = 0.03, label='pol. degree vector norm scale 0.03 ', labelpos='E')
                          # plt.text(-17*size[0]//25., 3*size[1]//2, im_name_lst[3], color='cyan',
                          #            fontsize='x-small', ha='center')      
                    plt.text(size[0]//10, 2*pix2mas*size[1]//6., 
                                f'{star_name}' + '_psf_' + f'{im_name_lst[j-1]}', color='w',
                              fontsize='large', ha='center')
                      #plt.colorbar(label='ADU in log$_{10}$ scale')
                    
                else :
                    plt.imshow(np.log10(sub_v_arr3[1]), cmap ='inferno', origin='lower',vmin=Vmin3[1], 
                                         vmax=Vmax3[1], extent = [x_min , x_max, y_min , y_max])   
                    plt.colorbar(label='ADU in log$_{10}$ scale')       
                    q_ = plt.quiver(X[::X_step,::X_step],Y[::X_step,::X_step],U3[::X_step,::X_step], V3[::X_step,::X_step])
                    plt.quiverkey(q_, X = -0.05, Y = 1.03, U = 0.03, label='pol. degree vector norm scale 0.03 ', labelpos='E')
                          # plt.text(-17*size[0]//25., 3*size[1]//2, im_name_lst[3], color='cyan',
                          #         fontsize='x-small', ha='center') 
                          
                          
                    plt.text(size[0]//10, 2*pix2mas*size[1]//6., 
                        f'{star_name}' + '_psf_' + f'{im_name_lst[1]}'+ 'Pol. vect', color='w',
                                  fontsize='large', ha='center')
            else:
                if np.any(np.min(sub_v_arr3[j])<= 0):           
                    plt.imshow(sub_v_arr3[j], cmap='inferno', origin='lower',
                    vmin=Vmin3[j], vmax=Vmax3[j], extent = [x_min , x_max, y_min , y_max])
                      
                    plt.text(size[0]//10, 2*pix2mas*size[1]//6., 
                                f'{star_name}' + '_psf_' + f'{im_name_lst[j]}', color='w',
                              fontsize='large', ha='center')
                    plt.colorbar(label='ADU in log$_{10}$ scale')
                else:
                    plt.imshow(np.log10(sub_v_arr3[j]), cmap='inferno', origin='lower',
                    vmin=Vmin3[j], vmax=Vmax3[j], extent = [x_min , x_max, y_min , y_max])
                       
                    plt.text(size[0]//10, 2*pix2mas*size[1]//6., 
                                  f'{star_name}' + '_psf_' + f'{im_name_lst[j]}', color='w',
                                fontsize='large', ha='center')
                    plt.colorbar(label='ADU in log$_{10}$ scale')
                    
                 
                    
        
        for k in range(len(nDimfigk)):      
              plt.subplot(3,3,nDimfigk[k] + 1)
              plt.plot(r_mas, np.log10(mean_sub_v_arr2[k]), color='darkorange',
                      linewidth = 2, label= f'{star_name}') 
              plt.plot(r_mas, np.log10(mean_sub_v_arr3[k]),color='purple',
                      linewidth = 2, label = f'{star_name}' + '_psf') 
              plt.legend(loc=0) 
              plt.xlabel('r (mas)', size=10) 
              if k == 0:
                  plt.ylabel(r'Intensity in log$_{10}$ scale', size=10)
        
        # plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+
        #                 '/plots/'+star_name+'_' +lst_fltr3[j] + '.pdf', 
        #                 dpi=100, bbox_inches ='tight')
        
        
        # plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+
        #                 '/plots/'+star_name+'_' +lst_fltr3[j] + '.png', 
        #                 dpi=100, bbox_inches ='tight')
        # plt.tight_layout()
    
    msg='reduction okay for '+ star_name
    return(msg)


    
