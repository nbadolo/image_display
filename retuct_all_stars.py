#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 14:34:30 2022

@author: nbadolo
"""
import os 
import log_agb_images
import log_agb_images_wp
import log_irc_images


lst_star = os.listdir('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/')
lst_star.sort()
lst_len = len(lst_star)
print(lst_star)

for i in range(lst_len): # affiche uniquement irc_10420
   if i == 8 :    
       log_irc_images.log_image(lst_star[i], 'alone')
       log_irc_images.log_image(lst_star[i], 'both')   
   
   else :  # Affiche toutes les étoiles qui n'ont pas de psf
       if i == 0 or i == 2 or i == 5 or i == 17 or i == 21 or i == 22:
           log_agb_images_wp.log_image(lst_star[i], 'alone')
           log_agb_images_wp.log_image(lst_star[i], 'both')
       
       else: # Affiche toutes les étoiles normales
           log_agb_images.log_image(lst_star[i], 'alone')
           log_agb_images.log_image(lst_star[i], 'both')
           