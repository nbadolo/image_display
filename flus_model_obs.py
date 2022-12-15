#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 13:59:14 2022

@author: nbadolo
"""


import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

star_name = 'Y_Scl'
file_path_r = '/home/nbadolo/Bureau/Aymard/Donnees_sph/pyssed_log/magn/'+star_name +'.ods'
df_r = pd.read_excel(file_path_r)

# Fixing random state for reproducibility


fig = plt.figure('Flux as a function of wavelengh')
fig.set_size_inches(18.5, 10, forward = True)
plt.plot(np.log10(df_r["wavel"]), np.log10(df_r["model"]), 'b+')
plt.plot(np.log10(df_r["wavel"]),np.log10(df_r["flux"]),  'r+')
plt.legend(["flux modelisé", "fux observé"], prop={'size': 22})
plt.title('Flux observé et modelisé de ' +star_name,  size=20)
plt.xlabel('lambda', size=20)
plt.ylabel('Flux', size=20)

plt.savefig('/home/nbadolo/Bureau/Aymard/These/for_papers/plots/log_paper_plots/' + star_name+ '_flux_model_obs.pdf', 
                dpi=100, bbox_inches ='tight')


plt.savefig('/home/nbadolo/Bureau/Aymard/These/for_papers/plots/log_paper_plots/' + star_name+ '_flux_model_obs.png', 
                dpi=100, bbox_inches ='tight')
plt.tight_layout()
plt.show()
#%%

fig = plt.figure('Flux as a function of wavelengh_')
fig.set_size_inches(18.5, 10, forward = True)
plt.plot(df_r["wavel"], (df_r["model"]), 'b+')
plt.plot(df_r["wavel"], (df_r["flux"]),  'r+')
plt.legend(["flux modelisé", "fux observé"], prop={'size': 22})
plt.title('Flux observé et modelisé de 17_Lep')
plt.xlabel('lambda', size=20)
plt.ylabel('Flux', size=20)