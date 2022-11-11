#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 11:10:49 2022

@author: nbadolo
"""

#packages
import numpy as np
from matplotlib import pyplot as plt
from math import pi, cos, sin, atan
from astropy.nddata import Cutout2D
from astropy.io import fits
from scipy.stats import linregress
import scipy.optimize as opt
import scipy.ndimage
from sklearn.linear_model import LinearRegression
from skimage import io, color, measure, draw, img_as_bool
import pylab as plt
from AymardPack import EllRadialProf as erp



#import matplotlib.pyplot as plt

star_name = 'SW_Col'
#file
fdir='/home/nbadolo/Bureau/Aymard/Donnees_sph/'
fdir_star=fdir+'log/'+star_name+'/star/both/V_N_R/'
#fdir_star=fdir+'log/Y_Scl/star/both/V_N_R/'
#fdir_psf=fdir+'psf/HD204971_mira/'
fname1='zpl_p23_make_polar_maps-ZPL_SCIENCE_P23_REDUCED'
fname2='-zpl_science_p23_REDUCED'
file_I_star= fdir_star + fname1 +'_I'+fname2+'_I.fits'
file_PI_star= fdir_star+fname1+'_PI'+fname2+'_PI.fits'
file_DOLP_star= fdir_star+fname1+'_DOLP'+fname2+'_DOLP.fits'
file_AOLP_star= fdir_star+ fname1+'_AOLP'+fname2+'_AOLP.fits'


# file_I_psf= fdir_psf + fname1+'_I'+fname2 + '_I.fits'
# file_PI_psf= fdir_psf+fname1+'_PI'+fname2 +'_PI.fits'
# file_DOLP_psf= fdir_psf+ fname1+'_DOLP'+fname2+'_DOLP.fits'
# file_AOLP_psf= fdir_psf+fname1+'_AOLP'+fname2+'_AOLP.fits'

#--creating lists
file_lst = [file_I_star, file_PI_star]
          #,file_I_psf,file_PI_psf,file_DOLP_psf,file_AOLP_psf]
          
nFrames = len(file_lst)
lst_Frame_name = ['Intensity', 'Pol_Intensity'] 

#--parameters
nDim=1024
nSubDim = 200 # plage de pixels que l'on veut afficher
size = (nSubDim, nSubDim)
nDimfigj=[9,10,11]
nDimfigk=[0,1,2]
lst_threshold = [0.01, 0.015, 0.02, 0.03, 0.05, 0.07, 0.1, 0.5]
n_threshold = len(lst_threshold)
vmin0 = 3.5
vmax0 = 15
pix2mas = 6.8  #en mas/pix
x_min = -pix2mas*nSubDim//2
x_max = pix2mas*(nSubDim//2-1)
y_min = -pix2mas*nSubDim//2
y_max = pix2mas*(nSubDim//2-1)
position = (nDim//2,nDim//2)
size = (nSubDim, nSubDim)

x, y = np.meshgrid(np.arange(nSubDim),np.arange(nSubDim)) #cree un tableau 
sub_v_arr = np.zeros((nFrames, nSubDim, nSubDim))
Ellips_arr = np.zeros((nFrames, n_threshold,nSubDim,nSubDim))
Ellips_im_arr = np.zeros((nFrames, n_threshold,nSubDim,nSubDim))
Ellips_arr2 = np.zeros((n_threshold, nSubDim, nSubDim))
Ellips_ = np.zeros((n_threshold, nSubDim))
Ell_rot_arr = np.zeros((nFrames, n_threshold, 2, nSubDim))
par_arr = np.zeros((nFrames, n_threshold, 5)) # les paramètres de l'étoile
Vmin = np.zeros((nFrames, n_threshold)) # pour l'image totale
Vmax = np.zeros((nFrames, n_threshold))
Vmin_r = np.zeros((nFrames, n_threshold))
Vmax_r = np.zeros((nFrames, n_threshold))# pour l'image reelle tronquee à au pourcentage
Vmin_w = np.zeros((nFrames, n_threshold))
Vmax_w = np.zeros((nFrames, n_threshold)) # pour l'image binaire

strs = [str(x*100) + ' %' for x in lst_threshold]


#fig_label = strs

#opening file
for i in range(nFrames):
      hdu = fits.open(file_lst[i])   
      data = hdu[0].data
      intensity = data[1,:,:]
      
      cutout = Cutout2D(intensity, position = position, size = size)
      zoom_hdu = hdu.copy()
      sub_v = cutout.data
      
      sub_v_arr[i] = sub_v
      
      for j in range(n_threshold):
          Ellips_arr2[i] = sub_v
          Ellips_im = np.zeros_like(sub_v_arr[i]) #creation d'un tableau de meme forme que sub_v
          Ellips_im[sub_v > lst_threshold[j]*np.max(sub_v)] = sub_v[sub_v > lst_threshold[j]*np.max(sub_v)]# on retient les points d'intensité égale à 5% de Imax 
          Ellips_im_arr[i][j] = Ellips_im
          
          Ellips = np.zeros_like(sub_v)          #creation d'un tableau de meme forme que sub_v
          Ellips[sub_v > lst_threshold[j]*np.max(sub_v)] = 1   # on retient les points d'intensité 
          Ellips_arr[i][j] = Ellips                # égale à 5% de Imax et à tous ces points 
                                                 # on donne 1 comme valeur d'intensité
    
            
          Vmin[i][j] = np.min(np.log10(sub_v+np.abs(np.min(sub_v))+10))
          Vmax [i][j] = np.max(np.log10(sub_v+np.abs(np.min(sub_v))+10))
          
          Vmin_r[i][j] = np.min(np.log10(Ellips_im+np.abs(np.min(Ellips_im))+10))
          Vmax_r[i][j] = np.max(np.log10(Ellips_im+np.abs(np.min(Ellips_im))+10))  
          
          Vmin_w[i][j] = np.min(np.log10(Ellips+np.abs(np.min(Ellips))+10))
          Vmax_w[i][j] = np.max(np.log10(Ellips+np.abs(np.min(Ellips))+10))  
     
      
          im_white = Ellips_arr[i][j]
          im_real = Ellips_im_arr[i][j]
          
          
        # =========================================================#
        # Deteminatiion des coordonnées du pixel le plus brillant  #        
        # =========================================================#
        
          index = np.argwhere(sub_v_arr[0])
          index_peak = np.argmax(index)
          index_peak = np.unravel_index(index_peak, 2)
          # X = index[:,1]
          # Y = index[:,0]
          print(index_peak)
          
         
          stop
          
        # =======================#
        #   Calcul de l'ellipse  #
        # =======================#
        
          regions = measure.regionprops(measure.label(im_white))
          bubble = regions[0]
          
          # initial guess (must be to change related on the % considered)
          y_i, x_i = bubble.centroid
          a_i = bubble.major_axis_length / 2.
          b_i = bubble.minor_axis_length / 2.
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
          # x_f = x_f *pix2mas
          # y_f = y_f*pix2mas
          # a_f = a_f*pix2mas 
          # b_f = b_f*pix2mas
          #def ellips(t, x_f, y_f, a_f, bb_f, theta_f):
          theta_f1 = theta_f    
          theta_f = np.pi/2 -theta_f
          par_arr[i][j] = [x_f, y_f, a_f, b_f, theta_f]
          
          theta_f_deg = theta_f*180/pi
          Ell = np.array([a_f*np.cos(t) , b_f*np.sin(t)])  # equation parametrique de l'ellips^
                 #u,v removed to keep the same center location
          M_rot = np.array([[cos(theta_f) , -sin(theta_f)],[sin(theta_f) , cos(theta_f)]])# matrice de rotation dans le sens anti-horaire
          
                 #2-D rotation matrix
            
          Ell_rot_ = np.zeros((2, nSubDim))
          Ell_rot = np.zeros((2, nSubDim))
          for k in range(Ell.shape[1]):
              Ell_rot_[:,k] = np.dot(M_rot,Ell[:,k]) # fait le produit scal de la matrice de rotation par chaq couple parametriq 
              Ell_rot[:,k] = Ell_rot_[:,k]
              Ell_rot[0,k] = Ell_rot[0,k] + x_f
              Ell_rot[1,k] = Ell_rot[1,k] + y_f
          #return Ell_rot.ravel() # .ravel permet de passer de deux dimension à une seule
              Ell_rot_arr[i][j][:,k] = Ell_rot[:,k] 
              #Ell_rot_arr[i][j][1] =  Ell_rot[1,k]
          
          # M_rot2 = np.array([[cos(alpha_rad) , -sin(alpha_rad)],[sin(alpha_rad) , cos(alpha_rad)]])
          # Ell_rot2 = np.zeros((2, Ell.shape[1]))
          # for k in range(Ell.shape[1]):
          #     Ell_rot2[:,k] = np.dot(M_rot2,Ell[:,k])
          #     Ell_rot2[0,k] = Ell_rot2[0,k] + x_f
          #     Ell_rot2[1,k] = Ell_rot2[1,k] + y_f
          
          
              
          plt.figure('white ellipse contour at ' + f'{strs[j]}' +' for ' + f'{lst_Frame_name[i]}')
          plt.clf()
          plt.imshow(np.log10(Ellips_arr[i][j]+np.abs(np.min(Ellips_arr[i][j]))+10), cmap ='inferno', vmin=Vmin_w[i][j], vmax=Vmax_w[i][j], origin='lower')
          #plt.plot( u + Ell_rot[0,:] , v + Ell_rot[1,:],'darkorange' )  #rotated ellipse
          plt.plot( Ell_rot[0,:] , Ell_rot[1,:],'darkorange' ) #rotated fit
          #plt.grid(color='lightgray',linestyle='--')
          #plt.show()
          plt.title('white ellipse contour at ' + f'{strs[j]}' + ' for ' + f'{lst_Frame_name[i]}'+' of ' + star_name, fontsize=10)
          plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+'/plots/fits/log_scale/semi_automatic/' +'white_ellips_contour_at_' + strs[j] + '_for_' + f'{lst_Frame_name[i]}' + '.pdf', 
                      dpi=100, bbox_inches ='tight')
          plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+'/plots/fits/log_scale/semi_automatic/' +'white_ellips_contour_at_' + strs[j] + '_for_' + f'{lst_Frame_name[i]}' + '.png', 
                      dpi=100, bbox_inches ='tight')
        
          
          
          # plt.figure('white ellipse contour for linear reg at ' + f'{strs[j]}')
          # plt.clf()
          # plt.imshow(np.log10(Ellips_arr[i]+np.abs(np.min(Ellips_arr[i]))+10), cmap ='inferno', vmin=Vmin_w[i], vmax=Vmax_w[i], origin='lower')
          #   #plt.plot( u + Ell_rot[0,:] , v + Ell_rot[1,:],'darkorange' )  #rotated ellipse
          # plt.plot( Ell_rot2[0,:] , Ell_rot2[1,:],'darkorange' ) #rotated fit
          #   #plt.grid(color='lightgray',linestyle='--')
          # plt.plot( xx , yy,'r' )
          # plt.show()
          # plt.title('white ellipse contour at ' + f'{strs[j]}')
          # plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/SW_Col/plots/fits/log_scale/' +'white_ellips_contour_at_' + strs[j] + '.pdf', 
          #             dpi=100, bbox_inches ='tight')
          
          
          
          # plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/SW_Col/plots/fits/log_scale/' +'white_ellips_contour_at_' + strs[j] +  '.png', 
          #     dpi=100, bbox_inches ='tight')
          # plt.tight_layout()
          
       
          # stop
          # plt.title('real image contour at ' + f'{strs[j]}' + ' for ' + f'{lst_Frame_name[i]}' +' of ' + star_name, fontsize=10)
          # plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+'/plots/fits/log_scale/semi_automatic/' +'real_image_contour_at_' + strs[j] + '_for_' + f'{lst_Frame_name[i]}' + '.pdf', 
          #               dpi=100, bbox_inches ='tight')
          # plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+'/plots/fits/log_scale/semi_automatic/' +'real_image_contour_at_' + strs[j] + '_for_' + f'{lst_Frame_name[i]}' + '.png', 
          #       dpi=100, bbox_inches ='tight')
          # plt.tight_layout()


         
          
          plt.figure('real image contour at ' + f'{strs[j]}')
          plt.clf()
          plt.imshow(np.log10(Ellips_im_arr[i][j] + np.abs(np.min(Ellips_im_arr[i][j]))+10), cmap ='inferno', vmin = Vmin_r[i][j], vmax = Vmax_r[i][j], origin='lower')
        
          plt.plot(Ell_rot[0,:] , Ell_rot[1,:],'darkorange' ) 
          
          #plt.show()
          plt.title('real image contour at ' + f'{strs[j]}' + ' for ' + f'{lst_Frame_name[i]}' +' of ' + star_name, fontsize=10)
          plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+'/plots/fits/log_scale/semi_automatic/' +'real_image_contour_at_' + strs[j] + '_for_' + f'{lst_Frame_name[i]}' + '.pdf', 
                        dpi=100, bbox_inches ='tight')


          plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+'/plots/fits/log_scale/semi_automatic/' +'real_image_contour_at_' + strs[j] + '_for_' + f'{lst_Frame_name[i]}' + '.png', 
                dpi=100, bbox_inches ='tight')
          plt.tight_layout()

          plt.figure('full image and contour at ' + f'{strs[j]}')
          plt.clf()
          plt.imshow(np.log10(sub_v_arr[i]+np.abs(np.min(sub_v_arr[i]))+10), cmap ='inferno', vmin=Vmin_w[i][j], vmax=Vmax_r[i][j], origin='lower')
          #plt.plot( u + Ell_rot[0,:] , v + Ell_rot[1,:],'darkorange' )  #rotated ellipse
          plt.plot(Ell_rot[0,:], Ell_rot[1,:],'darkorange' ) #rotated fit
          #plt.grid(color='lightgray',linestyle='--')
          #plt.show()
          plt.title('full image and contour at ' + f'{strs[j]}' + ' for ' + f'{lst_Frame_name[i]}' +' of '+ star_name, fontsize=10)
          plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+'/plots/fits/log_scale/semi_automatic/' +'full_image_and_contour_at_' + strs[j] + '_for_' + f'{lst_Frame_name[i]}' + '.pdf', 
                       dpi=100, bbox_inches ='tight')
        
        
          plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+'/plots/fits/log_scale/semi_automatic/' +'full_image_and_contour_at_' + strs[j] + '_for_' + f'{lst_Frame_name[i]}' + '.png', 
               dpi=100, bbox_inches ='tight')
          plt.tight_layout()


      plt.figure('full image and all the  contours' + f'{strs[j]}' + ' for ' f'{lst_Frame_name[i]}')
      plt.clf()
      plt.imshow(np.log10(sub_v_arr[i]+np.abs(np.min(sub_v_arr[i]))+10), cmap ='inferno', vmin = Vmin_w[i][j], vmax = Vmax_r[i][j], origin='lower')
      for j in range(n_threshold): 
          
           plt.plot(Ell_rot_arr[i][j][0,:], Ell_rot_arr[i][j][1,:])
          
      plt.show()
      plt.title('full image and all the contours ' + ' for ' + f'{lst_Frame_name[i]}' +' of '+ star_name, fontsize=10)
      plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+'/plots/fits/log_scale/semi_automatic/' +'full_image_and_all_the_contours' +'_for_' + f'{lst_Frame_name[i]}' + '.pdf', 
                      dpi=100, bbox_inches ='tight')
        
        
      plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+'/plots/fits/log_scale/semi_automatic/' +'full_image_and_all_the_contours'  + '_for_' + f'{lst_Frame_name[i]}' + '.png', 
              dpi=100, bbox_inches ='tight')
      plt.tight_layout()
      
      # ========================================================#
      # For the radial profile at a given orientation, theta_f  #     
      # ========================================================# 
       
     
      # imp = np.log10(sub_v_arr[1]+ np.abs(np.min(sub_v_arr[1])) + 10)
      # for i in range(nFrames):
      #     im = np.log10(sub_v_arr[i]+ np.abs(np.min(sub_v_arr[i])) + 10)
      #     x0, y0, x1, y1, x2, y2,z, zi1, zi2, x1_, y1_, x2_, y2_, zi1_, zi2_ = erp(par_arr[i][0][0], par_arr[i][0][1], par_arr[i][0][2], par_arr[i][0][3], par_arr[i][0][4], im, 100)
          
          
      #     plt.figure('radial profile at a given orientation, theta and theta + pi/2 ('+star_name+ ')' ' for ' f'{lst_Frame_name[i]}')
      #     plt.clf()
      #     plt.subplot(2,2,1)    
      #     plt.imshow(z, cmap ='inferno', vmin=Vmin_w[i][0], vmax = Vmax_r[i][0], origin='lower')# extent = [x_min , x_max, y_min , y_max] )
      #     plt.plot(Ell_rot_arr[i][0][0,:] , Ell_rot_arr[i][0][1,:])
      #     plt.plot([x0, x1], [y0, y1], 'ro-')
      #     plt.plot([x0, x1_], [y0, y1_], 'ro-')
      #     plt.xlabel('Relative R.A.(pix)', size=10)
      #     plt.ylabel('Relative Dec.(pix)', size=10)
      #     plt.title('radial profile at a given orientation, theta ('+star_name+ ')', fontsize=10)
      #     plt.subplot(2,2,2)    
      #     plt.imshow(z, cmap ='inferno', vmin = Vmin_w[i][0], vmax =  Vmax_r[i][0], extent = [x_min , x_max, y_min , y_max],  origin='lower')# extent = [x_min , x_max, y_min , y_max] )
      #     plt.plot(Ell_rot_arr[i][0][0,:] , Ell_rot_arr[i][0][1,:])
      #     plt.plot([x0, x2], [y0, y2], 'ro-')
      #     plt.plot([x0, x2_], [y0, y2_], 'ro-')
      #     plt.xlabel('Relative R.A.(pix)', size=10)
      #     plt.ylabel('Relative Dec.(pix)', size=10)
      #     plt.title('radial profile at a given orientation, theta + pi/2 ('+star_name+ ')', fontsize=10)
      #     plt.subplot(2,2,3)
      #     plt.plot(zi1)
      #     plt.plot(zi1_)
      #     plt.xlabel('r (mas)', size=10)          
      #     plt.ylabel(r'Intensity in log$_{10}$ scale', size=10)
      #     plt.subplot(2,2,4)
      #     plt.plot(zi2) 
      #     plt.plot(zi2_) 
      #     plt.xlabel('r (pix)', size=10) 
      #     plt.ylabel(r'Intensity in log$_{10}$ scale', size=10)
      #     plt.show()
      
      
      
      
      im = np.log10(sub_v_arr[i] + np.abs(np.min(sub_v_arr[i])) + 10)
      imp = np.log10(sub_v_arr[1] + np.abs(np.min(sub_v_arr[1])) + 10)
      x0, y0, x1, y1, x2, y2,z, zi1, zi2, x1_, y1_, x2_, y2_, zi1_, zi2_ = erp(par_arr[0][0][0], par_arr[0][0][1], par_arr[0][0][2], par_arr[0][0][3], par_arr[0][0][4], im, 100)
      x0p, y0p, x1p, y1p, x2p, y2p,zp, zi1p, zi2p, x1p_, y1p_, x2p_, y2p_, zi1p_, zi2p_ = erp(par_arr[1][0][0], par_arr[1][0][1], par_arr[1][0][2], par_arr[1][0][3], par_arr[1][0][4], imp, 100)
      plt.figure('profiles comparison')
      plt.clf()
      plt.subplot(2,2,1)
      plt.imshow(z, cmap ='inferno', vmin = Vmin_w[0][0], vmax = Vmax_r[0][0], origin='lower')
      plt.plot(Ell_rot_arr[0][0][0,:] , Ell_rot_arr[0][0][1,:])
      plt.plot([x0, x1], [y0, y1], 'ro-')
      plt.plot([x0, x1_], [y0, y1_], 'ro-')
      plt.plot([x0, x2], [y0, y2], 'ro-')
      plt.plot([x0, x2_], [y0, y2_], 'ro-')
      plt.xlabel('Relative R.A.(pix)', size=10)
      plt.ylabel('Relative Dec.(pix)', size=10)
      plt.title('radial profile at a given orientation, theta ('+star_name+ ')'  f'{lst_Frame_name[0]}', fontsize=10)
      plt.subplot(2,2,2)    
      plt.imshow(zp, cmap ='inferno', vmin = Vmin_w[1][0], vmax =  Vmax_r[1][0], origin='lower')# extent = [x_min , x_max, y_min , y_max] )
      plt.plot(Ell_rot_arr[1][0][0,:] , Ell_rot_arr[1][0][1,:])
      plt.plot([x0p, x1p], [y0p, y1p], 'ro-')
      plt.plot([x0p, x1p_], [y0p, y1p_], 'ro-')
      plt.plot([x0p, x2p], [y0p, y2p], 'ro-')
      plt.plot([x0p, x2p_], [y0p, y2p_], 'ro-')
      plt.xlabel('Relative R.A.(pix)', size=10)
      plt.ylabel('Relative Dec.(pix)', size=10)
      plt.title('radial profile at a given orientation, theta + pi/2 (' + star_name + ')'  f'{lst_Frame_name[1]}', fontsize=10)
      plt.subplot(2,2,3)
      plt.plot(zi1)
      plt.plot(zi1p)
      plt.legend(["Intensity", "Pol. Intensity"])
      plt.xlabel('r (mas)', size=10)          
      plt.ylabel(r'Intensity in log$_{10}$ scale', size=10)
      plt.subplot(2,2,4)
      plt.plot(zi2) 
      plt.plot(zi2p) 
      plt.legend(["Intensity", "Pol. Intensity"])
      plt.xlabel('r (pix)', size=10) 
      plt.ylabel(r'Intensity in log$_{10}$ scale', size=10)
      plt.show()
      plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+'/plots/radial/radial_profile_at_a_given_orientation.pdf', 
                      dpi=100, bbox_inches ='tight')
      plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+'/plots/radial/radial_profile_at_a_given_orientation.png', 
              dpi=100, bbox_inches ='tight')
      plt.tight_layout()