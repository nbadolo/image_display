#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 09:45:40 2022

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
file_lst=[file_PI_star]
          #,file_I_psf,file_PI_psf,file_DOLP_psf,file_AOLP_psf]
nFrames=len(file_lst)
#%%
#parameters
nDim=1024
nSubDim = 200 # plage de pixels que l'on veut afficher
size = (nSubDim, nSubDim)
nDimfigj=[9,10,11]
nDimfigk=[0,1,2]
lst_threshold = [0.01] #, 0.015, 0.02, 0.03, 0.05, 0.07, 0.1]
n_threshold = len(lst_threshold)
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
sub_v_arr = np.zeros((nFrames,nSubDim,nSubDim))
Ellips_arr = np.zeros((nFrames,nSubDim,nSubDim))
Ellips_im_arr = np.zeros((nFrames,nSubDim,nSubDim))
par_arr = np.zeros((n_threshold, 5)) # les paramètres de l'étoile
Vmin = np.zeros((nFrames)) # pour l'image totale
Vmax = np.zeros((nFrames))
Vmin_r = np.zeros((nFrames))
Vmax_r = np.zeros((nFrames))# pour 'imagge reel tronquee à au pourcentage
Vmin_w = np.zeros((nFrames))
Vmax_w = np.zeros((nFrames)) # pour l'image binaire

fmt = {}
strs = ['1%'] #, '1.5%', '2%', '3%', '5%', '7%','10%']
for l, s in zip(lst_threshold, strs):
    fmt[l] = s

#fig_label = strs
#%%
#opening file
for i in range(nFrames):
      hdu = fits.open(file_lst[i])   
      data = hdu[0].data
      intensity = data[1,:,:]
      
      cutout = Cutout2D(intensity, position=position, size=size)
      zoom_hdu = hdu.copy()
      sub_v = cutout.data
      
      sub_v_arr[i] = sub_v
      
      for j in range(n_threshold):
          Ellips_im = np.zeros_like(sub_v_arr[i]) #creation d'un tableau de meme forme que sub_v
          Ellips_im[sub_v > lst_threshold[j]*np.max(sub_v)] = sub_v[sub_v > lst_threshold[j]*np.max(sub_v)]# on retient les points d'intensité égale à 5% de Imax 
          Ellips_im_arr[i] = Ellips_im
          
          Ellips = np.zeros_like(sub_v)          #creation d'un tableau de meme forme que sub_v
          Ellips[sub_v > lst_threshold[j]*np.max(sub_v)] = 1   # on retient les points d'intensité 
          Ellips_arr[i] = Ellips                # égale à 5% de Imax et à tous ces points 
                                                 # on donne 1 comme valeur d'intensité
    
            
          Vmin[i] = np.min(np.log10(sub_v+np.abs(np.min(sub_v))+10))
          Vmax [i] = np.max(np.log10(sub_v+np.abs(np.min(sub_v))+10))
          
          Vmin_r[i] = np.min(np.log10(Ellips_im+np.abs(np.min(Ellips_im))+10))
          Vmax_r[i] = np.max(np.log10(Ellips_im+np.abs(np.min(Ellips_im))+10))  
          
          Vmin_w[i] = np.min(np.log10(Ellips+np.abs(np.min(Ellips))+10))
          Vmax_w[i] = np.max(np.log10(Ellips+np.abs(np.min(Ellips))+10))  
     
      
          im_white = Ellips_arr[i]
          im_real = Ellips_im_arr[i]

          
           # plot of the white image at 5% *Imax
          # plt.figure('white ellips')
          # plt.clf()
          # plt.imshow(im_white, cmap ='inferno', vmin=Vmin_w[i], vmax=Vmax_w[i], origin='lower')
          # plt.title('Intensity at ' + f'{strs[j]}')
          
         
          
            # plot of the real image at 5% *Imax
          # plt.figure('real image')
          # plt.clf()
          # plt.imshow(im_real, cmap ='inferno', vmin=Vmin_r[i], vmax=Vmax_r[i], origin='lower')
          # plt.title('Intensity at 5%')


          #image = img_as_bool(io.imread('bubble.jpg')[..., 0])
          regions = measure.regionprops(measure.label(im_white))
          bubble = regions[0]
          
          # initial guess (must be to change related on the % considered)
          y_i, x_i = bubble.centroid
          a_i = 0.5*bubble.major_axis_length / 2.
          b_i = 0.5*0.75*bubble.major_axis_length / 2.
          theta_i  = pi/4
          t = np.linspace(0, 2*pi, nSubDim)

          def cost(params):
              x0, y0, a, b, theta = params
              #coords = draw.disk((y0, x0), r, shape=image.shape)
              coords = draw.ellipse(y0, x0, a, b, shape=None, rotation= theta)
              template = np.zeros_like(im_white)
              template[coords] = 1
              return -np.sum(template == im_white)
            
          x_f, y_f, a_f, b_f, theta_f = opt.fmin(cost, (x_i, y_i, a_i, b_i, theta_i))

          #def ellips(t, x_f, y_f, a_f, bb_f, theta_f):

          par_arr[j] = [x_f, y_f, a_f, b_f, theta_f]
          
          Ell = np.array([a_f*np.cos(t) , b_f*np.sin(t)])  
                 #u,v removed to keep the same center location
          M_rot = np.array([[cos(theta_f) , -sin(theta_f)],[sin(theta_f) , cos(theta_f)]])  
                 #2-D rotation matrix
            
          Ell_rot = np.zeros((2, Ell.shape[1]))
          for k in range(Ell.shape[1]):
              Ell_rot[:,k] = np.dot(M_rot,Ell[:,k])
              Ell_rot[0,k] = Ell_rot[0,k] + x_f
              Ell_rot[1,k] = Ell_rot[1,k] + y_f
          #return Ell_rot.ravel() # .ravel permet de passer de deux dimension à une seule

          plt.figure('white ellipse contour at ' + f'{strs[j]}')
          plt.clf()
          plt.imshow(np.log10(Ellips_arr[i]+np.abs(np.min(Ellips_arr[i]))+10), cmap ='inferno', vmin=Vmin_w[i], vmax=Vmax_w[i], origin='lower')
            #plt.plot( u + Ell_rot[0,:] , v + Ell_rot[1,:],'darkorange' )  #rotated ellipse
          plt.plot( Ell_rot[0,:] , Ell_rot[1,:],'darkorange' ) #rotated fit
            #plt.grid(color='lightgray',linestyle='--')
          plt.show()
          plt.title('white ellipse contour at ' + f'{strs[j]}')
          plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/SW_Col/plots/fits/log_scale/' +'white_ellips_contour_at_' + strs[j] + '.pdf', 
                      dpi=100, bbox_inches ='tight')
        
        
          plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/SW_Col/plots/fits/log_scale/' +'white_ellips_contour_at_' + strs[j] +  '.png', 
              dpi=100, bbox_inches ='tight')
          plt.tight_layout()

          plt.figure('real image contour at ' + f'{strs[j]}')
          plt.clf()
          plt.imshow(np.log10(Ellips_im_arr[i]+np.abs(np.min(Ellips_im_arr[i]))+10), cmap ='inferno', vmin = Vmin_r[i], vmax = Vmax_r[i], origin='lower')
        #plt.plot( u + Ell_rot[0,:] , v + Ell_rot[1,:],'darkorange' )  #rotated ellipse
          plt.plot( Ell_rot[0,:] , Ell_rot[1,:],'darkorange' ) #rotated fit
        #plt.grid(color='lightgray',linestyle='--')
          plt.show()
          plt.title('real image contour at ' + f'{strs[j]}')
          plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/SW_Col/plots/fits/log_scale/' +'real_image_contour_at_' + strs[j] + '.pdf', 
                        dpi=100, bbox_inches ='tight')


          plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/SW_Col/plots/fits/log_scale/' +'real_image_contour_at_' + strs[j] +  '.png', 
                dpi=100, bbox_inches ='tight')
          plt.tight_layout()


          plt.figure('full image and contour at ' + f'{strs[j]}')
          plt.clf()
          plt.imshow(np.log10(sub_v_arr[i]+np.abs(np.min(sub_v_arr[i]))+10), cmap ='inferno', vmin=Vmin_w[i], vmax=Vmax_r[i], origin='lower')
          #plt.plot( u + Ell_rot[0,:] , v + Ell_rot[1,:],'darkorange' )  #rotated ellipse
          plt.plot( Ell_rot[0,:] , Ell_rot[1,:],'darkorange' ) #rotated fit
          #plt.grid(color='lightgray',linestyle='--')
          plt.show()
          plt.title('full image and contour at ' + f'{strs[j]}')
          plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/SW_Col/plots/fits/log_scale/' +'full_image_and_contour_at_' + strs[j] + '.pdf', 
                      dpi=100, bbox_inches ='tight')
        
        
          plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/SW_Col/plots/fits/log_scale/' +'full_image_and_contour_at_' + strs[j] +  '.png', 
              dpi=100, bbox_inches ='tight')
          plt.tight_layout()
