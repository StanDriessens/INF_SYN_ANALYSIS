# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 12:57:01 2023

analye Ih experimetn events 

@author: sdr267
"""

import os 
from pathlib import Path
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#%%load in the data 
path_recov = r'C:\Users\sdr267\Documents\PhD\ProjectSynapticConnections\Ih_experiment\events\new_tool\human\control'
files = Path(path_recov).glob('*.csv')
# load in the files
dfs = list()
for f in files:
    data = pd.read_csv(f)
    data['file']= f.stem
    dfs.append(data)
df_human_control = pd.concat(dfs, ignore_index=True)
df_human_control['condition'] = 'control'

path_recov = r'C:\Users\sdr267\Documents\PhD\ProjectSynapticConnections\Ih_experiment\events\new_tool\human\zd'
files = Path(path_recov).glob('*.csv')
# load in the files
dfs = list()
for f in files:
    data = pd.read_csv(f)
    data['file']= f.stem
    dfs.append(data)
df_human_zd = pd.concat(dfs, ignore_index=True)
df_human_zd['condition'] = 'zd'

path_recov = r'C:\Users\sdr267\Documents\PhD\ProjectSynapticConnections\Ih_experiment\events\new_tool\mouse\control'
files = Path(path_recov).glob('*.csv')
# load in the files
dfs = list()
for f in files:
    data = pd.read_csv(f)
    data['file']= f.stem
    dfs.append(data)
df_mouse_control = pd.concat(dfs, ignore_index=True)
df_mouse_control['condition'] = 'control'

path_recov = r'C:\Users\sdr267\Documents\PhD\ProjectSynapticConnections\Ih_experiment\events\new_tool\mouse\zd'
files = Path(path_recov).glob('*.csv')
# load in the files
dfs = list()
for f in files:
    data = pd.read_csv(f)
    data['file']= f.stem
    dfs.append(data)
df_mouse_zd = pd.concat(dfs, ignore_index=True)
df_mouse_zd['condition'] = 'zd'
#%%concentante to tota 
df_total = pd.concat([df_human_control, df_human_zd, df_mouse_control, df_mouse_zd], axis=0 , ignore_index=True)
df_total['file'] = df_total['file_name'].apply(lambda x: os.path.basename(x))
#only select qualtiy pass

df_total_pass = df_total.loc[df_total['qc'] == 'p']

#load in the metadata 

metadata=pd.read_csv(r'C:/Users/sdr267/Documents/PhD/ProjectSynapticConnections/Ih_experiment/metadata.csv')
#triangulate distance using pythagoras 
metadata['distance'] = metadata['x_dif'] + metadata['y_dif']

#merge metadta with the ephys data

df_merge = pd.merge(df_total, metadata, on='file') 

#calculate uM/ms 

df_merge['speed']= df_merge['distance']/df_merge['0']

#get normalized speeds aswell

#%%
#plot some data distributions
import seaborn as sns
import matplotlib.pyplot as plt
mouse_data = df_merge[df_merge['species_x']=='mouse']
human_data = df_merge[df_merge['species_x']== 'human']

mouse_data_average = mouse_data.groupby(['condition', 'file_name', 'layer']).mean().reset_index()
human_data_average = human_data.groupby(['condition', 'file_name', 'layer']).mean().reset_index()

#plot average lines 
for i in mouse_data_average.file_name:
    temp = mouse_data_average[mouse_data_average['file_name'] == i]
    sns.pointplot(data=temp, x='condition', y = '0')
plt.savefig(r'C:\Users\sdr267\Documents\PhD\ProjectSynapticConnections\Ih_experiment\mouse_cells.eps')


for i in human_data_average.file_name:
    temp = human_data_average[human_data_average['file_name'] == i].reset_index()
    if temp.layer[0] == 'deep L3':
        sns.pointplot(data=temp, x='condition', y = '0', color = 'magenta')
    elif temp.layer[0] == 'supp L2/3':
        sns.pointplot(data=temp, x='condition', y = '0', color = 'b')
    
plt.savefig(r'C:\Users\sdr267\Documents\PhD\ProjectSynapticConnections\Ih_experiment\human_cells_new.eps')

sns.pointplot(data=human_data_average, x='condition', y='0', hue='layer')
    



colors_scheme = {'#C875C4', '#32cd32'}


plt.figure(figsize=(18,18))
sns.lineplot(x='condition', y='speed', data=mouse_data, palette=(colors_scheme))
sns.boxplot(x='condition', y='speed', data=mouse_data, medianprops=dict(color="black", alpha=1, linewidth=4), boxprops=dict(linewidth=1.5, facecolor='w', edgecolor='w', alpha=1), 
                 whiskerprops=dict(linewidth=4, color='k', alpha=1), capprops=dict(linewidth=4, color='k', alpha=1))
plt.ylabel('synaptic delay (ms)', fontsize=30)
plt.yticks(fontsize=30)
plt.xticks(fontsize=30)
plt.xlabel('condition', fontsize=30)
plt.tick_params(length=8)
plt.tick_params(width=4)




sns.boxplot(mouse_data, x='condition', y='0')
sns.boxplot(human_data, x='condition', y='0')

sns.pointplot(data=mouse_data, x='n_pulse', y='0', hue='condition', errorbar=('ci', 95), capsize=1)


plt.figure(figsize=(15,25))

unique_files = mouse_data['file_name'].unique()

for file in unique_files:
    file_data = mouse_data[mouse_data['file_name'] == file]
    sns.lineplot(data=file_data, x='condition', y='normalized_speed', marker='o', linewidth =5,
                 dashes=True, markers=True, linestyle='--', markersize=20)
plt.ylabel('normalized signaling speed', fontsize=20)
#plt.ylim(2.0, 6.0)
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
plt.xlabel('condition', fontsize=20)
plt.tick_params(length=8)
plt.tick_params(width=4)

plt.figure(figsize=(15,25))
    
unique_files = human_data['file_name'].unique()

for file in unique_files:
    file_data = human_data[human_data['file_name'] == file]
    sns.lineplot(data=file_data, x='condition', y='normalized_speed', marker='o', linewidth =5,
                 dashes=True, markers=True, linestyle='--', markersize=20, legend='full')
plt.ylabel('normalized signaling speed', fontsize=20)
#plt.ylim(2.0, 6.0)
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
plt.xlabel('condition', fontsize=20)
plt.tick_params(length=8)
plt.tick_params(width=4)


#%%calculate normalized speeds
avg_speeds_control = human_data[human_data['condition'] == 'control'].groupby('file_name')['speed'].mean()
human_data = human_data.merge(avg_speeds_control.reset_index(), on='file_name', suffixes=('', '_avg_control'))

human_data['normalized_speed'] = human_data['normalized_speed'] = human_data['speed'] / human_data['speed_avg_control']
#mouse
avg_speeds_control = mouse_data[mouse_data['condition'] == 'control'].groupby('file_name')['speed'].mean()
mouse_data = mouse_data.merge(avg_speeds_control.reset_index(), on='file_name', suffixes=('', '_avg_control'))

mouse_data['normalized_speed'] = mouse_data['normalized_speed'] = mouse_data['speed'] / mouse_data['speed_avg_control']

human_zd = human_data.loc[human_data['condition']== 'zd'].reset_index()
mouse_zd = mouse_data.loc[mouse_data['condition']== 'zd'].reset_index()

human_zd['species']= 'human'
mouse_zd['species']= 'mouse'

total_zd = pd.concat([mouse_zd, human_zd], axis=0)


sns.lmplot(data=human_data, x='distance', y = 'speed', hue='condition', ci = None)



sns.lmplot(data=human_zd, x='distance', y = '0', ci = None, markers='x', palette='k')
plt.ylabel('normalized signaling speed', fontsize=15)
plt.xlabel('distance from stimulation electrode to soma (um)', fontsize =15)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)

sns.lmplot(data=mouse_zd, x='distance', y = 'normalized_speed', ci = None)

sns.lmplot(data=total_zd, x='distance', y = 'normalized_speed', hue='species', ci=None)


#differences for pulse numbers\
mouse_data_pulse = df_total_pass.loc[df_total_pass['species_y'] == 'mouse'].groupby(['file_name', 'condition', 'n_pulse']).mean().reset_index()
human_data_pulse = df_total_pass.loc[df_total_pass['species_y'] == 'human'].groupby(['file_name', 'condition', 'n_pulse']).mean().reset_index()

sns.lineplot(data=mouse_data_pulse, x='n_pulse', y = '0', hue='condition', errorbar=('ci', 95), style="condition", markers=True)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)

plt.ylabel('signal latency (ms)', fontsize= 15)
plt.xlabel('external stimulation pulse number', fontsize= 15)

sns.lineplot(data=human_data_pulse, x='n_pulse', y = 'speed', hue='condition', errorbar=('ci', 0), style="condition", markers=True)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)

unique_files = human_data_pulse['file_name'].unique()

for file in unique_files:
    file_data = human_data_pulse[human_data_pulse['file_name'] == file]
    sns.lineplot(data=file_data, x='n_pulse', y='speed', marker='o', linewidth =1,
                 dashes=True, markers=True, linestyle='--', markersize=1, legend='full', hue = 'condition')


plt.ylabel('signal latency (ms)', fontsize= 15)
plt.xlabel('external stimulation pulse number', fontsize= 15)

#%%comapre speeds 


df_control_total = df_total_pass[df_total_pass['condition'] == 'control'].groupby(['file_name', 'species_x']).mean().reset_index()
 

colors_scheme = {'#C875C4', '#32cd32'}


plt.figure(figsize=(10,10))
sns.stripplot(x='species_x', y='speed', data=df_control_total, size=12, hue='species_x', palette=colors_scheme)
sns.boxplot(x='species_x', y='speed', data=df_control_total, medianprops=dict(color="black", alpha=1, linewidth=4), boxprops=dict(linewidth=1.5, facecolor='w', edgecolor='w', alpha=1), 
                 whiskerprops=dict(linewidth=4, color='k', alpha=1), capprops=dict(linewidth=4, color='k', alpha=1))
plt.ylabel('signaling speed (um/ms)', fontsize=30)
plt.yticks(fontsize=30)
plt.xticks(fontsize=30)
plt.xlabel('species', fontsize=30)
plt.tick_params(length=8)
plt.tick_params(width=4)


import pingouin as pg 

pg.mwu(df_control_total[df_control_total['species_x']=='human'].speed, df_control_total[df_control_total['species_x']=='mouse'].speed)



#%%plot the decay and rise times 
df_total.tau = (1/df_total.tau)/500000

df_total.tau = df_total.tau*1e3
df_total.rise_time_0_100 = df_total.rise_time_0_100*1e3
df_total.rise_time_10_90 = df_total.rise_time_10_90*1e3

#plot
df_human = df_total[df_total['species'] == 'human']
df_mouse = df_total[df_total['species'] == 'mouse']

human_data = df_human.groupby(['file', 'condition', 'layer']).mean().reset_index()
mouse_data = df_mouse.groupby(['file', 'condition', 'layer']).mean().reset_index()


human_data.rise_time_10_90 = human_data.rise_time_10_90 *1e3
human_data.rise_time_0_100 = human_data.rise_time_0_100 *1e3



mouse_data.rise_time_10_90 = mouse_data.rise_time_10_90 *1e3
mouse_data.rise_time_0_100 = mouse_data.rise_time_0_100 *1e3

human_data = mouse_data

human_data = human_data[human_data['tau'] < 150]

unique_files = np.unique(human_data.file)
for file in unique_files:
    file_data = human_data[human_data['file'] == file]
    sns.lineplot(data=file_data, x='condition', y='tau', marker='o', linewidth =5,
                 dashes=True, markers=True, linestyle='--', markersize=20)
plt.ylabel('Tau', fontsize=20)
#plt.ylim(2.0, 6.0)
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
plt.xlabel('condition', fontsize=20)
plt.tick_params(length=8)
plt.tick_params(width=4)

plt.figure(figsize=(15,25))
    
unique_files = human_data['file_name'].unique()

#%% all data

sns.lineplot(data=human_data, x='condition', y='tau', marker='o', linewidth =5,
             dashes=True, markers=True, linestyle='--', markersize=20, errorbar=(None))
plt.ylabel('Tau', fontsize=20)
plt.ylim(20, 55)





