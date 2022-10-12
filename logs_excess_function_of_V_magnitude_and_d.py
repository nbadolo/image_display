#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 10:02:04 2022

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
# Nous representons ici, pour le grand tableau log, l'échantillon et les objets  
# resolus les execes infra rouge et le ratio L_IR/L_str en fonction de la magnitude en V
# des étoiles  
# =============================================================================

"""
# Definitions bloc:
    file_path_l: chemin gros tableau
    file_path_s: chemin du petit tableau
    file_path_r: chemin des étoiles resolues
"""

file_path_l = '/home/nbadolo/Bureau/Aymard/These/for_papers/large_log.ods'
file_path_s = '/home/nbadolo/Bureau/Aymard/These/for_papers/short_log.ods'
file_path_r = '/home/nbadolo/Bureau/Aymard/These/for_papers/resol_log.ods'

df_l = pd.read_excel(file_path_l)
df_s = pd.read_excel(file_path_s)
df_r = pd.read_excel(file_path_r)
# xmin = np.min(df["J-Ks"]) 
# xmax = np.max(df["J-Ks"])
# ymin = np.min(df["E_IR"])
# ymax = np.max(df["E_IR"])
#%%
fig = plt.figure('Excess  as a function of the V-magnitude and d')
plt.clf()
fig.set_size_inches(18.5, 10, forward = True)
plt.subplot(2,2,1)
plt.plot(df_l["V"], df_l["E_IR"], 'b+')
plt.plot(df_s["V"], df_s["E_IR"], 'r+')
plt.plot(df_r["V"], df_r["E_IR"], 'k+')
plt.legend(["large table", "the sample", "resolved objects"])
plt.title('E_IR = f(V)')
plt.xlabel('V', size=10)
plt.ylabel('E_IR', size=10)

plt.subplot(2,2,2)
plt.plot(df_l["V"], df_l["E12"], 'b+')
plt.plot(df_s["V"], df_s["E12"], 'r+')
plt.plot(df_r["V"], df_r["E12"], 'k+')
plt.legend(["large table", "the sample", "resolved objects"])
plt.title('E12 = f(V)')
plt.xlabel('V', size=10)
plt.ylabel('E12', size=10)


plt.subplot(2,2,3)
plt.plot(df_l["d"], np.log10(df_l["E_IR"]), 'b+')
plt.plot(df_s["d"], np.log10(df_s["E_IR"]), 'r+')
plt.plot(df_r["d"], np.log10(df_r["E_IR"]), 'k+')
plt.legend(["large table", "the sample", "resolved objects"])
plt.title('E_IR = f(d)')
plt.xlabel('d', size=10)
plt.ylabel('E_IR', size=10)

plt.subplot(2,2,4)
plt.plot(df_l["V"], df_l["LIR/L*"], 'b+')
plt.plot(df_s["V"], df_s["LIR/L*"], 'r+')
plt.plot(df_r["V"], df_r["LIR/L*"], 'k+')
plt.legend(["large table", "the sample","resolved objects"])
plt.title('LIR/L* = f(V)')
plt.xlabel('V', size=10)
plt.ylabel('LIR/L*', size=10)

#plt.title( label = "Excess infrared in relation to key parameters" ,fontsize = 12 ,color = "k" )

plt.savefig('/home/nbadolo/Bureau/Aymard/These/for_papers/plots/log_paper_plots/' + 'E_IR_f_V_magn.pdf', 
                dpi=100, bbox_inches ='tight')


plt.savefig('/home/nbadolo/Bureau/Aymard/These/for_papers/plots/log_paper_plots/' + 'E_IR_f_V_magn.png', 
                dpi=100, bbox_inches ='tight')
plt.tight_layout()
plt.show()
