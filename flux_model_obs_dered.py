#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 13:59:14 2022

@author: nbadolo
"""


import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

"""
Pour la comparaison des spectres des étoiles de mon echantillon issus des fichiers *.sed du code de Iain:
 flux modelisé, flux observé et flux observé dérougi   
"""

lst_star = os.listdir('/home/nbadolo/Bureau/Aymard/Donnees_sph/log/')
print(lst_star)
n_lst_star = len(lst_star)
for i in range (n_lst_star):
    if i != 3 and i != 13:

    #star_name = ['SW_Col', '17_Lep']
    #i= 1
    
        file_path_r = '/home/nbadolo/Bureau/Aymard/Donnees_sph/pyssed_log/magn/'+lst_star[i] + '.ods'
        df_r = pd.read_excel(file_path_r)
        
        # Fixing random state for reproducibility
                  
        fig = plt.figure('Flux as a function of wavelengh (echelle log)')
        plt.clf()
        fig.set_size_inches(18.5, 10, forward = True)
        plt.plot(np.log10(df_r["wavel"]), np.log10(df_r["model"]), 'b+')
        plt.plot(np.log10(df_r["wavel"]),np.log10(df_r["flux"]),  'r+')
        plt.plot(np.log10(df_r["wavel"]),np.log10(df_r["dered"]),  'g+')
        plt.legend(["flux modelisé", "fux observé", "fux observé dérougi"], prop={'size': 22})
        plt.title('Flux observé et modelisé de ' +lst_star[i],  size=20)
        plt.xlabel('lambda', size=20)
        plt.ylabel('Flux', size=20)
        
        plt.savefig('/home/nbadolo/Bureau/Aymard/These/for_papers/plots/log_paper_plots/log_scale/log_scale' + lst_star[i]+ '_flux_model_obs.pdf', 
                        dpi=100, bbox_inches ='tight')
        
        
        plt.savefig('/home/nbadolo/Bureau/Aymard/These/for_papers/plots/log_paper_plots/log_scale/log_scale' + lst_star[i] + '_flux_model_obs.png', 
                        dpi=100, bbox_inches ='tight')
        plt.tight_layout()
        plt.show()    
        fig = plt.figure('Flux as a function of wavelengh (echelle lineaire)')
        plt.clf()
        fig.set_size_inches(18.5, 10, forward = True)
        plt.plot(df_r["wavel"]/10000, (df_r["model"]), 'bo')
        plt.plot(df_r["wavel"]/10000, (df_r["flux"]),  'ro')
        plt.plot(df_r["wavel"]/10000, (df_r["dered"]),  'g^')
        plt.legend(["flux modelisé", "fux observé","fux observé dérougi"], prop={'size': 22})
        plt.title('Flux observé et modelisé de ' + lst_star[i],  size=20)     
        plt.xlabel('lambda(µm)', size=20)
        plt.ylabel('Flux(W/m^2)', size=20)
        plt.xlim(0.7,20,100)
        plt.ylim(0,500,100)
        stop
        plt.savefig('/home/nbadolo/Bureau/Aymard/These/for_papers/plots/log_paper_plots/linear_scale/linear_scale' + lst_star[i]+ '_flux_model_obs.pdf', 
                        dpi=100, bbox_inches ='tight')
        
        
        plt.savefig('/home/nbadolo/Bureau/Aymard/These/for_papers/plots/log_paper_plots/linear_scale/linear_scale' + lst_star[i]+ '_flux_model_obs.png', 
                        dpi=100, bbox_inches ='tight')
        plt.tight_layout()
        plt.show()
