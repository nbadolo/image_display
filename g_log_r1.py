#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 13:35:53 2022

@author: nbadolo
"""
# =============================================================================
# Lecture du tableau excel des étoiles massives
# =============================================================================

import xlrd
import pandas as pd
from xlrd import*
from openpyxl import load_workbook
import matplotlib.pyplot as plt



file_path = '/home/nbadolo/Bureau/Aymard/These/for_papers/hip.xls'

df = pd.read_excel(file_path)
# cols = df.columns.tolist() #mettre les éléments des colonnes dans une liste
# cols.sort() #ordonné cette liste

plt.figure('B-V = f(d)')
plt.plot(df["d"], df["B-V"])
plt.show()


fig, axs = plt.subplots(2, 2)

axs[0, 0].plot(df["d"], df["B-V"], 'b')
axs[0,0].text(0.4, 0.8,'B-V = f(d)', horizontalalignment='center', verticalalignment='center', transform = axs[0, 0].transAxes, color='m')
#plt.text(50, 5,'B-V = f(d)', color='m',
   #   fontsize='large', ha='center')
axs[0, 1].plot(df["d"], df["U"], 'k')
axs[0,1].text(0.8, 0.8,'U = f(d)', horizontalalignment='center', verticalalignment='center', transform = axs[0, 1].transAxes, color='m')
#plt.text(0.2, 0.8,'eV = f(d)', color='m',
     # fontsize='large', ha='center')
axs[1, 0].plot(df["d"], df["B"], 'r--')
axs[1,0].text(0.8, 0.8,'B = f(d)', horizontalalignment='center', verticalalignment='center', transform = axs[1, 0].transAxes, color='m')
# plt.text(0.2, 0.8,'e[B-V] = f(B-V)', color='m',
#       fontsize='large', ha='center')
axs[1, 1].plot(df["B-V"], df["U"], )
axs[1,1].text(0.8, 0.8,'U = f(B-V)', horizontalalignment='center', verticalalignment='center', transform = axs[1, 1].transAxes, color='m')
# plt.text(0.2, 0.8,'e[U-B] = f(U-B)', color='m',
#       fontsize='large', ha='center')
fig.suptitle(' Relations entre certains paramètres des étoiles du tableau excel', fontsize=11)


plt.savefig('/home/nbadolo/Bureau/Aymard/These/for_papers/plots/' + 'test1.pdf', 
                dpi=100, bbox_inches ='tight')


plt.savefig('/home/nbadolo/Bureau/Aymard/These/for_papers/plots/' + 'test1.png', 
                dpi=100, bbox_inches ='tight')
plt.tight_layout()
plt.show()