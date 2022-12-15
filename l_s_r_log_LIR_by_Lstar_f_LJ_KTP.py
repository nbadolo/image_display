#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 14:03:50 2022

@author: nbadolo
"""

import xlrd
import pandas as pd
from xlrd import*
from openpyxl import load_workbook
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# Nous representons ici, pour le grand tableau log, l'excès à 12 micron en fonction 
# de quelques paramètres clés: la luminosité, le J-K, la temperature effective et 
# la periode  
# =============================================================================

"""
# Definitions bloc:
    file_path_l: chemin gros tableau
    file_path_s: chemin du petit tableau
    file_path_r: chein des étoiles resolues
"""

file_path_l = '/home/nbadolo/Bureau/Aymard/These/for_papers/large_log_par.ods'
file_path_s = '/home/nbadolo/Bureau/Aymard/These/for_papers/short_log_par.ods'
file_path_r = '/home/nbadolo/Bureau/Aymard/These/for_papers/resol_log_par.ods'

df_l = pd.read_excel(file_path_l)
df_s = pd.read_excel(file_path_s)
df_r = pd.read_excel(file_path_r)
# xmin = np.min(df["J-Ks"]) 
# xmax = np.max(df["J-Ks"])
# ymin = np.min(df["E_IR"])
# ymax = np.max(df["E_IR"])
#%%
fig = plt.figure('LIR/L* ratio in relation to key parameters')
fig.set_size_inches(18.5, 10, forward = True)
plt.subplot(2,2,1)
plt.plot(df_l["L"], df_l["E12"], 'b+')
plt.plot(df_s["L"], df_s["E12"], 'r+')
plt.plot(df_r["L"], df_r["E12"], 'ok')
plt.legend(["large table", "the sample", "resolved"])
plt.xlabel('L', size=10)
plt.ylabel('LIR/L*', size=10)

plt.subplot(2,2,2)
plt.plot(df_l["J-K"], df_l["LIR/L*"], 'b+')
plt.plot(df_s["J-K"], df_s["LIR/L*"], 'r+')
plt.plot(df_r["J-K"], df_r["LIR/L*"], 'ok')
plt.legend(["large table", "the sample", "resolved"])
plt.xlabel('J-K', size=10)
plt.ylabel('LIR/L*', size=10)


plt.subplot(2,2,3)
plt.plot(df_l["Teff"], df_l["LIR/L*"], 'b+')
plt.plot(df_s["Teff"], df_s["LIR/L*"], 'r+')
plt.plot(df_r["Teff"], df_r["LIR/L*"], 'ok')
plt.legend(["large table", "the sample", "resolved"])
plt.xlabel('Teff', size=10)
plt.ylabel('LIR/L*', size=10)

plt.subplot(2,2,4)
plt.plot(df_l["P"], df_l["LIR/L*"], 'b+')
plt.plot(df_s["P"], df_s["LIR/L*"], 'r+')
plt.plot(df_r["P"], df_r["LIR/L*"], 'ok')
plt.legend(["large table", "the sample", "resolved"])
plt.xlabel('P', size=10)
plt.ylabel('LIR/L*', size=10)

#plt. title( label = "Excess infrared in relation to key parameters" ,fontsize = 12 ,color = "k" )

plt.savefig('/home/nbadolo/Bureau/Aymard/These/for_papers/plots/log_paper_plots/' + 'l_s_r_log_LIR_Lstr.pdf', 
                dpi=100, bbox_inches ='tight')


plt.savefig('/home/nbadolo/Bureau/Aymard/These/for_papers/plots/log_paper_plots/' + 'l_s_r_log_LIR_Lstr.png', 
                dpi=100, bbox_inches ='tight')
plt.tight_layout()
plt.show()