#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 14:34:30 2022

@author: nbadolo
"""
import os 
import log_agb_images_alone

lst_star = os.listdir('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/')
lst_star.sort()
lst_len = len(lst_star)


for i in range(lst_len): #range(len(lst_star)):
   # print('Reduct_Numb', i)
    log_agb_images_alone.log_image(lst_star[i])