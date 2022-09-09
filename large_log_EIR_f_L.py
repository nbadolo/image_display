#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 10:46:02 2022

@author: nbadolo
"""

import xlrd
import pandas as pd
from xlrd import*
from openpyxl import load_workbook
import matplotlib.pyplot as plt
import numpy as np


file_path_l = '/home/nbadolo/Bureau/Aymard/These/for_papers/large_log_par.ods'
file_path_s = '/home/nbadolo/Bureau/Aymard/These/for_papers/short_log_par.ods'

df_l = pd.read_excel(file_path_l)
df_s = pd.read_excel(file_path_s)
# xmin = np.min(df["J-Ks"]) 
# xmax = np.max(df["J-Ks"])
# ymin = np.min(df["E_IR"])
# ymax = np.max(df["E_IR"])
fig, axs = plt.subplots(2, 2)

axs[0, 0].plot(df_l["L"], df_l["E_IR"], 'b+')
axs[0, 0].plot(df_s["L"], df_s["E_IR"], 'r+')
axs[0,0].text(0.4, 0.8,'E_IR = f(L)', horizontalalignment='center', verticalalignment='center', transform = axs[0, 0].transAxes, color='m')

#   fontsize='large', ha='center')

axs[0, 1].plot(df_l["J-Ks"], df_l["E_IR"], 'k+')
axs[0,1].text(0.8, 0.8,'E_IR = f(J-Ks)', horizontalalignment='center', verticalalignment='center', transform = axs[0, 1].transAxes, color='m')
plt.xlabel( 'J-Ks',size=15)
plt.ylabel('E_IR',size=15)

axs[1, 0].plot(df_l["Teff"], df_l["E_IR"], 'r--')
axs[1,0].text(0.4, 0.8,'E_IR = f(Teff)', horizontalalignment='center', verticalalignment='center', transform = axs[1, 0].transAxes, color='m')
# plt.text(0.2, 0.8,'e[B-V] = f(B-V)', color='m',
#       fontsize='large', ha='center')
plt.xlabel( 'Teff',size=15)
plt.ylabel('E_IR',size=15)

axs[1, 1].plot(df_l["P"], df_l["E_IR"], 'c+')
axs[1,1].text(0.8, 0.8,'E_IR = f(P)', horizontalalignment='center', verticalalignment='center', transform = axs[1, 1].transAxes, color='m')
# plt.text(0.2, 0.8,'e[U-B] = f(U-B)', color='m',
#       fontsize='large', ha='center')
plt.xlabel( 'Period',size=15)
plt.ylabel('E_IR',size=15)

fig.suptitle(' Relations entre certains paramètres des étoiles du tableau excel', fontsize=11)


plt.savefig('/home/nbadolo/Bureau/Aymard/These/for_papers/plots/log_paper_plots/' + 'large_log_E_IR.pdf', 
                dpi=100, bbox_inches ='tight')


plt.savefig('/home/nbadolo/Bureau/Aymard/These/for_papers/plots/log_paper_plots/' + 'large_log_E_IR.png', 
                dpi=100, bbox_inches ='tight')
plt.tight_layout()
plt.show()