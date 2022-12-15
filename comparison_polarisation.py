#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 09:21:29 2022

@author: nbadolo
"""



"""
Code simplifié pour la comparaison des deux methodes de calcul du Qphi et l'évaluation de leur impact 
sur les vecteurs de polaristaion. Resultats à envoyer à Miguel :
"""
#packages
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
import Julien_image_tool as imtool


#%%
star_name = 'SW_Col'  
obsmod = 'both'     
fdir= '/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+ '/'
fdir_star = fdir + 'star/'+obsmod+ '/V_N_R/' 

#fdir_star_fltr = fdir + fdir_star
fname1='zpl_p23_make_polar_maps-ZPL_SCIENCE_P23_REDUCED'
fname2='-zpl_science_p23_REDUCED'
file_Q_star= fdir_star + fname1+'_Q'+fname2+'_Q.fits'
file_U_star= fdir_star +fname1+'_U'+fname2+'_U.fits'
file_Q_phi_star= fdir_star +fname1+'_Q_PHI'+fname2+'_Q_PHI.fits'

file_lst = [file_Q_star, file_U_star, file_Q_phi_star]
nFrames = len(file_lst)



nDim = 1024 
nSubDim = 60 # plage de pixels que l'on veut afficher
size = (nSubDim, nSubDim)

pix2mas = 3.4  #en mas/pix
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



i_v_arr = np.empty((nFrames,nDim,nDim))   
sub_v_arr = np.empty((nFrames,nSubDim,nSubDim))


#im_name_lst = ['I','PI','DOLP', 'AOLP']
Vmin = np.empty((nFrames))
Vmax = np.empty((nFrames))

position = (nDim//2,nDim//2)


x, y = np.meshgrid(np.arange(nSubDim), np.arange(nSubDim)) #cree un tableau 

R = np.sqrt((x-nSubDim/2)**2+(y-nSubDim/2)**2)
r = np.linspace(1,nSubDim//2-1,nSubDim//2-1) # creation d'un tableau de distance radiale

r_mas = pix2mas*r #  où r est en pixels et r_mas en millièmes d'arcseconde

ang_arr2 = imtool.angle_array((nSubDim,nSubDim),centerx=None,centery=None)  #calcul of the positional angle selon de Boer
ang_arr1 = imtool.angleN_array((nSubDim,nSubDim),centerx=None,centery=None)  #calcul of the positional angle selon Schmid


# opening Q and U
for i in range (nFrames):
    hdu = fits.open(file_lst[i])[0]   
    Data = hdu.data   
    i_v = Data[0,:,:]
    i_v_arr[i] = i_v.data
    
    #fltr = hdu.header.get('HIERARCH ESO INS3 OPTI5 NAME')     
    #print(fltr)                   
    cutout = Cutout2D(i_v, position=position, size=size)
    #zoom_hdu = hdu.copy()
    sub_v = cutout.data
    sub_v_arr[i] = sub_v
    
    

Q = sub_v_arr[0]
U = sub_v_arr[1]
Q_phi_sphr = sub_v_arr[2]

phi1=ang_arr1 #  le postional angle selon Schmid, calcul inspiré par le code de Julien
phi2=ang_arr2 #  le postional angle selon de Boer et calculer par Julien

Co2 = np.cos(2*phi2) 
Si2 = np.sin(2*phi2)
Q_phi = -Q*Co2 - U*Si2 

Co1 = np.cos(2*phi1) 
Si1 = np.sin(2*phi1)
Q_r = Q*Co1 + U*Si1  # plarisation radiale de  Schmid 2006 

Vmin =[np.min(Q_phi),np.min(Q_phi_sphr), np.min(Q_r)]
Vmax = [np.max(Q_phi),np.max(Q_phi_sphr), np.max(Q_r)] 

plt.figure('Q_phi')
plt.clf()
plt.subplot(1,2,1)
plt.imshow((Q_phi), cmap = 'inferno', origin='lower', vmin=Vmin[0], 
                vmax=Vmax[0], extent = [x_min , x_max, y_min , y_max])
plt.colorbar(label='ADU in log$_{10}$ scale', shrink = 0.6)
plt.text(size[0]//10, 2*pix2mas*size[1]//6.,
           ' Q_phi,' +' De Boer 2020'+ '  ('+f'{star_name}'+ ')', color='w',
       fontsize='large', ha='center')
plt.xlabel('mas', size=14)
plt.ylabel('mas', size=14)
#plt.title('polarisation azimutale de Boer 2020', size = 12)

# plt.subplot(2,2,2)
# plt.imshow((Q_phi_sphr), cmap = 'inferno', origin='lower', vmin=Vmin[1], 
#                 vmax=Vmax[1], extent = [x_min , x_max, y_min , y_max])
# plt.colorbar(label='ADU in log$_{10}$ scale', shrink = 0.6)
# plt.xlabel('mas', size=14)
# plt.ylabel('mas', size=14)
# plt.title('polarisation azimutale de Sphere', size = 12)

plt.subplot(1,2,2)
plt.imshow((Q_r), cmap = 'inferno', origin='lower', vmin=Vmin[2], 
                vmax=Vmax[2], extent = [x_min , x_max, y_min , y_max])
plt.text(size[0]//10, 2*pix2mas*size[1]//6.,
           ' Q_r,' +' Schmid 2006'+ '  (' + f'{star_name}' + ')', color='w',
    fontsize='large', ha='center')
plt.colorbar(label='ADU in log$_{10}$ scale', shrink = 0.6)
plt.xlabel('mas', size=14)
plt.ylabel('mas', size=14)
#plt.title('plarisation radiale de  Schmid 2006', size = 12)

plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/Donnees_Miguel/Q_phi_boer.pdf', dpi=100, bbox_inches ='tight')
plt.tight_layout()


stop 

plt.figure('Q_r')
plt.clf()
plt.imshow((Q_r), cmap = 'inferno', origin='lower')
plt.xlabel('pix', size=20)
plt.ylabel('pix', size=20)
plt.title('plarisation radiale de  Schmid 2006', size = 18)
#%%              
  #calcul de Q_phi

x_pix, y_pix = np.meshgrid(i_v.x, i_v .y) # recupère les coordonnées de l'image en pixel et en fait une grille complete



Q = np.reshape(Q,(nSubDim,nSubDim))
U = np.reshape(U,(nSubDim,nSubDim))

tphi1 = x_pix/y_pix 
#phi1 = np.arctan(tphi1) # positional angle selon de Schmid 2006
Co1 = np.cos(2*phi1) 
Si1 = np.sin(2*phi1)
Q_r = Q*Co1 + U*Si1  # plarisation radiale de  Schmid 2006
#%%
fig, ax = plt.subplots() 
plt.figure('Q_r')
plt.clf()
plt.imshow(np.abs(Q_r), cmap = 'inferno', origin='lower')
plt.xlabel('pix', size=20)
plt.ylabel('pix', size=20)
plt.title('plarisation radiale de  Schmid 2006', size = 18)
plt.savefig('/home/nbadolo/SIM_CODES/RADMC3D/newest_version/radmc3d-2.0/AymardModels/PolData/pol_rad_pyr_mg100.pdf', dpi=100, bbox_inches ='tight')
plt.tight_layout()
#%%
tphi2 = -x_pix/y_pix                      # calcul de la tangente du positional angle selon de Boer2020
phi2 = np.arctan(tphi2)                   # calcul du positional angle selon de Boer2020
Co2 = np.cos(2*phi2) 
Si2 = np.sin(2*phi2)
Q_phi = -Q*Co2 - U*Si2  



#%%        
plt.clf()
fig = plt.figure(f'{star_name}' +'(' + f'{fltr}' + '_Cam1' + ')')
fig.set_size_inches(12, 6, forward = True)










          
plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+
                  '/plots/'+ star_name +'_' + fltr + '_Cam1' + '.pdf', 
                  dpi=100, bbox_inches ='tight')


plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/'+star_name+
                  '/plots/'+ star_name +'_' + fltr + '_Cam1' + '.png', 
                  dpi=100, bbox_inches ='tight')
plt.tight_layout()
msg1='reduction okay for '+ star_name+'_Cam1'
 #return(msg1)
print(msg1)