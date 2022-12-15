#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 18:18:04 2022

@author: nbadolo
"""



import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

file_path_l = '/home/nbadolo/Bureau/Aymard/These/for_papers/used_tables/large_log.ods'
file_path_s = '/home/nbadolo/Bureau/Aymard/These/for_papers/used_tables/short_log.ods'
file_path_r = '/home/nbadolo/Bureau/Aymard/These/for_papers/used_tables/resol_log.ods'

df_l = pd.read_excel(file_path_l)
df_s = pd.read_excel(file_path_s)
df_r = pd.read_excel(file_path_r)





# Fixing random state for reproducibility


#colors = np.random.rand(N)
fact = 10000
area_l = 1/df_l["d"]*fact 
area_s = 1/df_s["d"]*fact
area_r = 1/df_r["d"]*fact


fig = plt.figure('Excess  as a function of  d and the V-magnitude')
plt.clf()
fig.set_size_inches(18.5, 10, forward = True)
plt.scatter(df_l["V"], df_l["E_IR"], s=area_l, c='y', marker= "*" )
plt.scatter(df_s["V"], df_s["E_IR"], s=area_s, c='g' , marker= "*")
plt.scatter(df_r["V"], df_r["E_IR"], s=area_r, c='b' )
plt.legend(["McDonald et al. 2017, 2012", "our sample", "resolved objects"], prop={'size': 22})
plt.title('The width of the dots is proportional to the distance from the star which is <= 300pc', fontsize= 20)
plt.xlabel('V', size=20)
plt.ylabel('E_IR', size=20)

plt.savefig('/home/nbadolo/Bureau/Aymard/These/for_papers/plots/log_paper_plots/' + 'E_IR_f_d_V_magn.pdf', 
                dpi=100, bbox_inches ='tight')


plt.savefig('/home/nbadolo/Bureau/Aymard/These/for_papers/plots/log_paper_plots/' + 'E_IR_f_d_V_magn.png', 
                dpi=100, bbox_inches ='tight')
plt.tight_layout()
plt.show()

#%%
N = 45
x, y = np.random.rand(2, N)
c = np.random.randint(1, 5, size=N)
s = np.random.randint(10, 220, size=N)

fig, ax = plt.subplots()

scatter = ax.scatter(x, y, c=c, s=s)

# produce a legend with the unique colors from the scatter
legend1 = ax.legend(*scatter.legend_elements(),
                    loc="lower left", title="Classes")
ax.add_artist(legend1)

# produce a legend with a cross section of sizes from the scatter
handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6)
legend2 = ax.legend(handles, labels, loc="upper right", title="Sizes")
ax.xaxis.set_tick_params(labelsize=24)
ax.yaxis.set_tick_params(labelsize=24)
print(handles)
print(labels)
plt.show()

#%%
# open file
file_path_t = '/home/nbadolo/Bureau/Aymard/These/for_papers/used_tables/test_log.ods'
df_t = pd.read_excel(file_path_t)
N = len(df_t["d"])
print(N)
area_t = 1/df_t["d"]*fact
s = area_t 
c = np.random.randint(1, 5, size = N)
#x, y = np.random.rand(2, 575)
x, y = df_t["V"], df_t["E_IR"]
fig, ax = plt.subplots()

scatter = ax.scatter(x, y, c='y', s=s)
handles, labels = scatter.legend_elements(prop="sizes", num = 4, alpha = 0.6)
legend2 = ax.legend(handles, labels, loc="upper right", title="Sizes")
print(handles)
print(labels)

# là c'est okay
#%%  A travailler !!!!!
plt.clf()
fig.set_size_inches(18.5, 10, forward = True)
plt.scatter(df_t["V"], df_t["E_IR"], s=area_t, c='g') #, marker= "*" )
plt.scatter(1.5*df_t["V"], 2*df_t["E_IR"], s=area_t*2, c='y') # , marker= "*")
# plt.scatter(df_r["V"], df_r["E_IR"], s=area_r, c='b')
plt.legend(["McDonald et al. 2017, 2012", "our sample"], prop={'size': 22})
#ax.add_artist(legend1)
legend2 = ax.legend(handles, labels, loc="upper right", title="Sizes")
#legend2 = plt.legend(handles, labels, loc="upper right", title="Sizes")


#fig, ax = plt.subplots()


# produce a legend with a cross section of sizes from the scatter

