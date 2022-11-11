#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 10:19:54 2022

@author: nbadolo
"""


#packages
from __future__ import print_function
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
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as mpl
from scipy.optimize import curve_fit


#import matplotlib.pyplot as plt

star_name = 'SW_Col'
#file
fdir='/home/nbadolo/Bureau/Aymard/Donnees_sph/'
fdir_star = fdir+'log/'+star_name+'/star/both/V_N_R/'
#fdir_star = fdir+ 'log/Y_Scl/star/both/V_N_R/'
#fdir_psf = fdir+'psf/HD204971_mira/'
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
nFrames = len(file_lst)


nDim = 1024
nSubDim = 200 # plage de pixels que l'on veut afficher
size = (nSubDim, nSubDim)
pix2mas = 3.6  #en mas/pix
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


mean_sub_v_arr = np.empty((nFrames,nSubDim//2-1))
sub_v_arr = np.empty((nFrames,nSubDim,nSubDim))
im_name_lst = ['I','PI']
Vmin = np.empty((nFrames))
Vmax = np.empty((nFrames))


position = (nDim//2,nDim//2)
size = (nSubDim, nSubDim)

x, y = np.meshgrid(np.arange(nSubDim), np.arange(nSubDim)) #cree un tableau 

R = np.sqrt((x-nSubDim/2)**2+(y-nSubDim/2)**2)
r = np.linspace(1,nSubDim//2-1,nSubDim//2-1)

r_mas=pix2mas*r #  où r est en pixels et r_mas en millièmes d'arcseconde

for i in range (nFrames):
      hdu = fits.open(file_lst[i])   
      data = hdu[0].data   
      i_v = data[0,:,:]
   
      cutout = Cutout2D(i_v, position=position, size=size)
      zoom_hdu = hdu.copy()
      sub_v = cutout.data
    
      f = lambda r : sub_v[(R >= r-0.5) & (R < r + 0.5)].mean()   
      mean_sub_v = np.vectorize(f)(r) 
    
      mean_sub_v_arr[i] = mean_sub_v 
      sub_v_arr[i]=sub_v
      if i==3 or i==7:
          Vmin[i]=np.min(sub_v_arr[i])
          Vmax[i]=np.max(sub_v_arr[i])  
      else:
          Vmin[i]=np.min(np.log10(sub_v_arr[i]))
          Vmax[i]=np.max(np.log10(sub_v_arr[i]))  

      
      mean_sub = np.argwhere(mean_sub_v)
      index = mean_sub
      xdata = index[:,0]
      ydata = index[:,0]
    
      #--Recast xdata and ydata into numpy arrays so we can use their handy features
      xdata = np.asarray(xdata)
      ydata = np.asarray(ydata)
      plt.figure(0)
      plt.clf()
      plt.plot(xdata, ydata, 'o')
      stop
      #--Define the Gaussian function
      def Gauss(x, A, B):
    	   y = A*np.exp(-1*B*x**2)
    	   return y
      parameters, covariance = curve_fit(Gauss, xdata, ydata)
    
      fit_A = parameters[0]
      fit_B = parameters[1]
    
      fit_y = Gauss(xdata, fit_A, fit_B)
      plt.plot(xdata, ydata, 'o', label='data')
      plt.plot(xdata, fit_y, '-', label='fit')
      plt.legend()


plt.figure(f'{star_name}' +'(' + f'{ im_name_lst[i]}' + ')', figsize=(40,20.5))
plt.clf()    
          
for k in range(nFrames):      
      plt.subplot(1, 2, k+1)
      plt.plot(r_mas, mean_sub_v_arr[k], color='m',
              linewidth=2, label='Mira') 
      plt.xlabel('r (mas)', size=10) 
      if k == 0:
          plt.ylabel(r'Intensity in log$_{10}$ scale', size=10)

# plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+
#                 '/plots/'+star_name+'_' +lst_fltr3[j] + '.pdf', 
#                 dpi=100, bbox_inches ='tight')


# plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+
#                 '/plots/'+star_name+'_' +lst_fltr3[j] + '.png', 
#                 dpi=100, bbox_inches ='tight')
plt.tight_layout()
