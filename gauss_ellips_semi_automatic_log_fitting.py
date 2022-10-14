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
#import matplotlib.pyplot as plt
#%%
star_name = 'SW_Col'
#file
fdir='/home/nbadolo//Bureau/Aymard/Donnees_sph/'
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
#%%
#creating lists
file_lst=[file_PI_star]
          #,file_I_psf,file_PI_psf,file_DOLP_psf,file_AOLP_psf]
          
nFrames = len(file_lst)
lst_Frame_name = ['Intensity', 'Pol_Intensity'] 
#%%
#parameters
nDim=1024
nSubDim = 200 # plage de pixels que l'on veut afficher
size = (nSubDim, nSubDim)
nDimfigj=[9,10,11]
nDimfigk=[0,1,2]
lst_threshold = [0.01, 0.015, 0.02, 0.03, 0.05, 0.07, 0.1]
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
#%%
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

          
           # plot of the white image at 5% *Imax
          # plt.figure('white ellips')
          # plt.clf()
          # plt.imshow(im_white, cmap ='inferno', vmin=Vmin_w[i], vmax=Vmax_w[i], origin='lower')
          # plt.title('Intensity at ' + f'{strs[j]}')
#%%                            
          #plot of the real image at 5% *Imax
          # plt.figure(3)
          # plt.clf()
          # plt.imshow(im_real, cmap ='inferno', vmin=Vmin_r[i], vmax=Vmax_r[i], origin='lower')
          # plt.title('Intensity at 5%')
#%%       
# utilisation d'une regression linéaire pour déterminer l'orientation

          # index = np.argwhere(im_white) # recupère les indices  des points dont l'intensité est non nulle 
          # X = index[:,1]
          # Y = index[:,0]
          
          
          # linear_reg = np.polyfit(X, Y, 1, full = False, cov = True)
          # alpha_rad = atan(linear_reg[0][0])    # recupération de la pente de la regression
          # alpha_deg = alpha_rad*180/pi
          # aa = linear_reg[0][0]
          # bb = linear_reg[0][1]
          # xx = np.arange(nSubDim)
          # yy = aa*xx + bb
          
          
          #slope, intercept, r, p, se = linregress(Y, X)
#%%          #image = img_as_bool(io.imread('bubble.jpg')[..., 0])
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
          
          #def ellips(t, x_f, y_f, a_f, bb_f, theta_f):
          theta_f1 = theta_f    
          theta_f = np.pi/2 -theta_f
          par_arr[i][j] = [x_f, y_f, a_f, b_f, theta_f]
          
          theta_f_deg = theta_f*180/pi
          Ell = np.array([a_f*np.cos(t) , b_f*np.sin(t)])  # equation parametrique de l'ellips^
                 #u,v removed to keep the same center location
          M_rot = np.array([[cos(theta_f) , -sin(theta_f)],[sin(theta_f) , cos(theta_f)]]) 
          
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
          
        # =============================================================================
        #  For the radial profile at a given orientation, theta_f       
        # =============================================================================
          
          
          #z  = sub_v_arr[i]
          z = np.log10(sub_v_arr[i]+np.abs(np.min(sub_v_arr[i]))+10)
          #z = sub_v_arr[i] + np.abs(np.min(sub_v_arr[i]) + 0.01
          #-- Extract the line...
          # Make a line with "num" points...

          # x1, y1 = x_f, y_f
          # x2, y2 = 0, y_f - theta_f*x_f
          
          x1, y1 = x_f, y_f      # common point of th two lines
          x2, y2 = x_f - y_f/(np.tan(theta_f)), 0 # the second point of theta line 
          x3, y3 = x_f + y_f/(np.tan(theta_f)), 0  # the second point of (theta + pi/2) line 
         
          
          num = 1000
          # x, y = np.linspace(x0, x1, num), np.linspace(y0, y1, num)
          x, y = np.linspace(x1, x2, num), np.linspace(y1, y2, num)
          x_, y_ = np.linspace(x1, x3, num), np.linspace(y1, y3, num)

          # Extract the values along the line, using cubic interpolation
          zi = scipy.ndimage.map_coordinates(z, np.vstack((x,y))) #  les données pour le profile radial
          zi_ = scipy.ndimage.map_coordinates(z, np.vstack((x_,y_)))
          #-- Plot...
          # fig, axes = plt.subplots(nrows=2)
          # axes[0].imshow(z)
          # axes[0].plot([x1, x2], [y1, y2], 'ro-')
          # axes[0].axis('image')

          # axes[1].plot(zi)
          
          plt.figure('radial profile at a given orientation, theta_f and theta_f + pi/2 ('+star_name+ ')')
          plt.clf()
          plt.subplot(2,2,1)    
          plt.imshow(z, cmap ='inferno', vmin=Vmin_w[i][j], vmax=Vmax_r[i][j], origin='lower')# extent = [x_min , x_max, y_min , y_max] )
          plt.plot(Ell_rot[0,:], Ell_rot[1,:])
          plt.plot([x1, x2], [y1, y2], 'ro-')
          plt.xlabel('Relative R.A.(pix)', size=10)
          plt.ylabel('Relative Dec.(pix)', size=10)
          plt.title('radial profile at a given orientation, theta_f ('+star_name+ ')', fontsize=10)
          plt.subplot(2,2,2)    
          plt.imshow(z, cmap ='inferno', vmin=Vmin_w[i][j], vmax=Vmax_r[i][j], origin='lower')# extent = [x_min , x_max, y_min , y_max] )
          plt.plot(Ell_rot[0,:], Ell_rot[1,:])
          plt.plot([x1, x3], [y1, y3], 'ro-')
          plt.xlabel('Relative R.A.(pix)', size=10)
          plt.ylabel('Relative Dec.(pix)', size=10)
          plt.title('radial profile at a given orientation, theta_f + pi/2 ('+star_name+ ')', fontsize=10)
          plt.subplot(2,2,3)
          plt.plot(zi)
          plt.xlabel('r (mas)', size=10)          
          plt.ylabel(r'Intensity in log$_{10}$ scale', size=10)
          plt.subplot(2,2,4)
          plt.plot(zi_) 
          plt.xlabel('r (pix)', size=10) 
          plt.ylabel(r'Intensity in log$_{10}$ scale', size=10)
          plt.show()
          
          
          stop
          plt.title('real image contour at ' + f'{strs[j]}' + ' for ' + f'{lst_Frame_name[i]}' +' of ' + star_name, fontsize=10)
          plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+'/plots/fits/log_scale/semi_automatic/' +'real_image_contour_at_' + strs[j] + '_for_' + f'{lst_Frame_name[i]}' + '.pdf', 
                        dpi=100, bbox_inches ='tight')
          plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+'/plots/fits/log_scale/semi_automatic/' +'real_image_contour_at_' + strs[j] + '_for_' + f'{lst_Frame_name[i]}' + '.png', 
                dpi=100, bbox_inches ='tight')
          plt.tight_layout()


         
          
          plt.figure('real image contour at ' + f'{strs[j]}')
          plt.clf()
          plt.imshow(np.log10(Ellips_im_arr[i][j] + np.abs(np.min(Ellips_im_arr[i][j]))+10), cmap ='inferno', vmin = Vmin_r[i][j], vmax = Vmax_r[i][j], origin='lower')
        #plt.plot( u + Ell_rot[0,:] , v + Ell_rot[1,:],'darkorange' )  #rotated ellipse
          plt.plot(Ell_rot[0,:] , Ell_rot[1,:],'darkorange' ) #rotated fit
        #plt.grid(color='lightgray',linestyle='--')
          
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


      plt.figure('full image and all the  contours' + f'{strs[j]}')
      plt.clf()
      plt.imshow(np.log10(sub_v_arr[i]+np.abs(np.min(sub_v_arr[i]))+10), cmap ='inferno', vmin=Vmin_w[i][j], vmax=Vmax_r[i][j], origin='lower')
      #plt.plot( Ell_rot_arr[i][0][0,:] , Ell_rot_arr[i][0][1,:],'darkorange' )
      #plt.plot( u + Ell_rot[0,:] , v + Ell_rot[1,:],'darkorange' )  #rotated ellipse
      for j in range(n_threshold): 
          
           plt.plot(Ell_rot_arr[i][j][0,:], Ell_rot_arr[i][j][1,:])
          
           #plt.plot( Ell_rot[0,:] , Ell_rot[1,:],'darkorange' ) #rotated fit
           #plt.grid(color='lightgray',linestyle='--')
      plt.show()
      plt.title('full image and all the contours ' + ' for ' + f'{lst_Frame_name[i]}' +' of '+ star_name, fontsize=10)
      plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+'/plots/fits/log_scale/semi_automatic/' +'full_image_and_all_the_contours' +'_for_' + f'{lst_Frame_name[i]}' + '.pdf', 
                      dpi=100, bbox_inches ='tight')
        
        
      plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+'/plots/fits/log_scale/semi_automatic/' +'full_image_and_all_the_contours'  + '_for_' + f'{lst_Frame_name[i]}' + '.png', 
              dpi=100, bbox_inches ='tight')
      plt.tight_layout()