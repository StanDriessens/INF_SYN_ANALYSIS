# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 08:57:04 2023

@author: sdr267
"""

from pathlib import Path
import seaborn as sns
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from pynwb import NWBFile, TimeSeries, NWBHDF5IO
from pynwb.epoch import TimeIntervals
from pynwb.file import Subject
from pynwb.behavior import SpatialSeries, Position
from datetime import datetime
from dateutil import tz
from ipfx.dataset.create import create_ephys_data_set
from ipfx.feature_extractor import SpikeFeatureExtractor
import h5py
import matplotlib.pyplot as plt 
from scipy.signal import order_filter
from scipy.signal import butter,filtfilt
import scipy as sp
import pandas as pd
import seaborn as sns
import os
from scipy.signal import find_peaks, peak_prominences
import h5py
import tkinter 
from tkinter import *
from tkinter import filedialog
from tkinter import messagebox as mb
from tkinter import Radiobutton

#%% SINGLE EVENT DETECION
print('select file from the list')
file = filedialog.askopenfilename()
condition = []
while condition != 'ZD' and condition != 'control':
    condition = input('which protocol are you analyzing?, ZD or control?')
    if condition == 'ZD':
        cond = 'TTL_test_ZD'
    elif condition == 'control':
        cond = 'TTL_test_gap'
    else:
        print('wrong input, try again')
species = input('please provide species (human/mouse)')

print('loading dataset')
dataset= create_ephys_data_set(file)
all_sweeps = dataset.filtered_sweep_table()
print('fetching TTL data')
#where is TTL??
f = h5py.File(file, 'r')
stim = pd.DataFrame(list(f['stimulus']['presentation'].keys()), columns=['ttl'])
#create list with only TTL file names
data_TTL = stim[stim['ttl'].str.contains('TTL')]
#select subset of sweeps 
subset = list(all_sweeps.loc[all_sweeps['stimulus_code']==cond].sweep_number)
#convert to string for search
strings = [str(x) for x in subset]
#only select TTLs for subselected sweep numbers 
selected_TTLs = data_TTL[data_TTL['ttl'].str.contains('|'.join(strings))].reset_index(drop=True)
"""
work in progress to make sure TTL matches the correct sweep incase wrong gap free was used 

works now 
"""


if condition == 'control':
    ttl_accept = []
    sweeps_to_analyze = []
    for i in selected_TTLs.ttl:
        temp = str(f['stimulus']['presentation'][i].attrs['stimulus_description'])
        if '10hz' in temp:
            sweep_number = f['stimulus']['presentation'][i].attrs['sweep_number']
            DA_sweep = all_sweeps.iloc[sweep_number]
            if 'TTL' in DA_sweep.stimulus_code and 'ZD' not in DA_sweep.stimulus_code and  'inbetwe' not in DA_sweep.stimulus_code:
                ttl_accept.append(i)
                sweeps_to_analyze.append(int(DA_sweep.sweep_number))
                
if condition == 'ZD':
    ttl_accept = []
    sweeps_to_analyze = []
    for i in selected_TTLs.ttl:
        temp = str(f['stimulus']['presentation'][i].attrs['stimulus_description'])
        if '10hz' in temp:
            sweep_number = f['stimulus']['presentation'][i].attrs['sweep_number']
            DA_sweep = all_sweeps.iloc[sweep_number]
            if 'ZD' in DA_sweep.stimulus_code  and  'inbetwe' not in DA_sweep.stimulus_code:
                ttl_accept.append(i)
                sweeps_to_analyze.append(int(DA_sweep.sweep_number))           




#select sweeps to analyze
#sweeps_to_analyze = all_sweeps.loc[all_sweeps['stimulus_code']==cond].sweep_number.reset_index(drop=True)

response_onset = pd.DataFrame()
for i in range(0,len(sweeps_to_analyze)):
    sweeps = []
    ttl = selected_TTLs.ttl[i]
    ttl = np.array(f['stimulus']['presentation'][ttl]['data'])
    sweep = dataset.sweep(sweeps_to_analyze[i])
    #filter the signal please
    nyq = 250*1e3
    cutoff = 2000
    order = 3
    fs = 500*1e3
    def butter_lowpass_filter(data, cutoff, fs, order):
        normal_cutoff = cutoff / nyq
         # Get the filter coefficients 
        b, a = butter(order, normal_cutoff, btype='low', analog=False)
        y = filtfilt(b, a, data)
        return y
    #find peaks in the ttl signal 
    peaks_ttl, nonesne = find_peaks(ttl, 1) 
    v_filt = butter_lowpass_filter(sweep.v, cutoff, fs, order)
    dv = np.gradient(v_filt, edge_order=1)
    dt = np.gradient(sweep.t, edge_order=1)
    dvdt = dv/dt
    fig,ax = plt.subplots(3,1, sharex=True)
    ax[0].plot(sweep.t*1e3, ttl)
    ax[0].plot(sweep.t[peaks_ttl]*1e3, ttl[peaks_ttl],  'x')
    ax[1].plot(sweep.t*1e3, v_filt)
    ax[2].plot(sweep.t*1e3, dvdt)
    t=sweep.t*1e3
    #find differences between TTL peak and second peak in the derivative after peak
    onsets = []
    for j in range(0,len(peaks_ttl)):
        left_border = peaks_ttl[j]+650
        right_border =  peaks_ttl[j]+2000
        aoi = dvdt[left_border:right_border]
        threshold = int(input("Please provide derravitve threshold for event detection: "))
        new_onset = []
        while threshold > max(aoi) or new_onset == 'f' :
            print('dvdt does not exceed threshold, please change thershold')
            threshold = int(input("Please provide derravitve threshold for event detection: "))
            onset_index = np.where(aoi>threshold)[0][0] 
            ax[1].plot(t[left_border:right_border][onset_index], v_filt[left_border:right_border][onset_index], 'x')
            print('printing new onset, press a button to continue')
            while True:
                if plt.waitforbuttonpress():
                    break 
            new_onset = input('New onset good? (p/f)')
        print('analyzing event, press a button to continue')
        onset = 'f'
        while onset == 'f':
            onset_index = np.where(aoi>threshold)[0][0] 
            ax[1].plot(t[left_border:right_border][onset_index], v_filt[left_border:right_border][onset_index], 'x')
            while True:
                if plt.waitforbuttonpress():
                    break 
            onset = input('Onset good? (p/f)')
            if onset == 'f':
                threshold = int(input("Please provide derravitve threshold for event detection: "))
        onsets.append(t[left_border:right_border][onset_index]-t[peaks_ttl[j]])
        sweeps.append(j)
    onsets_df = pd.DataFrame(onsets, sweeps)
onsets_df['condition'] = condition
onsets_df['file_name'] = f.filename
onsets_df['species']   = species

save_name = f.filename[82:101]+'.csv'
onsets_df.to_csv(r'C:\Users\sdr267\Documents\PhD\ProjectSynapticConnections\Ih_experiment\events\\' + save_name)

#%%AVERAGE OVER SWEEPS
#initialize analysis 

#open file
file = filedialog.askopenfilename()
#select species
species_select = ['human', 'mouse']
def save_species():
    global species
    species = species_select[x.get()]
    window.destroy()
window = Tk()
x = IntVar()
for index in range(len(species_select)):
    radiobutton = Radiobutton(window, text=species_select[index], variable=x, value=index, padx = 25)
    radiobutton.pack(anchor=W)
okay_button = Button(window, text="Select Species", command=save_species())
okay_button.pack(pady=10)    
window.mainloop()    



print('loading dataset')
dataset= create_ephys_data_set(file)
all_sweeps = dataset.filtered_sweep_table()

#get the correct sweep name from tkinter button
sets_unique = list(np.unique(all_sweeps.stimulus_code))

def save_selected_value():
    global cond
    cond = sets_unique[x.get()]
    window.destroy()
    

window = Tk()
x = IntVar()
for index in range(len(sets_unique)):
    radiobutton = Radiobutton(window, text=sets_unique[index], variable=x, value=index, padx = 25)
    radiobutton.pack(anchor=W)
okay_button = Button(window, text="Select sweep protocol", command=save_selected_value)
okay_button.pack(pady=10)    

window.mainloop()
    

print('fetching TTL data')
#where is TTL??
f = h5py.File(file, 'r')
stim = pd.DataFrame(list(f['stimulus']['presentation'].keys()), columns=['ttl'])
#create list with only TTL file names
data_TTL = stim[stim['ttl'].str.contains('TTL')]

#select subset of sweeps 
subset = list(all_sweeps.loc[all_sweeps['stimulus_code']==cond].sweep_number)
#convert to string for search
strings = [str(x) for x in subset]
#only select TTLs for subselected sweep numbers 
selected_TTLs = data_TTL[data_TTL['ttl'].str.contains('|'.join(strings))].reset_index(drop=True)
#select sweeps to analyze

if condition == 'control':
    ttl_accept = []
    sweeps_to_analyze = []
    for i in selected_TTLs.ttl:
        temp = str(f['stimulus']['presentation'][i].attrs['stimulus_description'])
        if '10hz' in temp:
            sweep_number = f['stimulus']['presentation'][i].attrs['sweep_number']
            DA_sweep = all_sweeps.iloc[sweep_number]
            if 'TTL' in DA_sweep.stimulus_code and 'ZD' not in DA_sweep.stimulus_code and  'inbetwe' not in DA_sweep.stimulus_code:
                ttl_accept.append(i)
                sweeps_to_analyze.append(int(DA_sweep.sweep_number))
                
if condition == 'ZD':
    ttl_accept = []
    sweeps_to_analyze = []
    for i in selected_TTLs.ttl:
        temp = str(f['stimulus']['presentation'][i].attrs['stimulus_description'])
        if '10hz' in temp:
            sweep_number = f['stimulus']['presentation'][i].attrs['sweep_number']
            DA_sweep = all_sweeps.iloc[sweep_number]
            if 'ZD' in DA_sweep.stimulus_code  and  'inbetwe' not in DA_sweep.stimulus_code:
                ttl_accept.append(i)
                sweeps_to_analyze.append(int(DA_sweep.sweep_number))              

fig,ax = plt.subplots(3,1, sharex=True)     
average_list = []
average_list_dvdt = []
for i in range(0,len(sweeps_to_analyze)):
    sweeps = []
    ttl = ttl_accept[i]
    ttl = np.array(f['stimulus']['presentation'][ttl]['data'])
    sweep = dataset.sweep(sweeps_to_analyze[i])
    #filter the signal please
    nyq = 250*1e3
    cutoff = 5000
    order = 2
    fs = 500*1e3
    def butter_lowpass_filter(data, cutoff, fs, order):
        normal_cutoff = cutoff / nyq
         # Get the filter coefficients 
        b, a = butter(order, normal_cutoff, btype='low', analog=False)
        y = filtfilt(b, a, data)
        return y
    #find peaks in the ttl signal 
    peaks_ttl, nonesne = find_peaks(ttl, 1) 
    v_filt = butter_lowpass_filter(sweep.v, cutoff, fs, order)
    dv = np.gradient(v_filt, edge_order=1)
    dt = np.gradient(sweep.t, edge_order=1)
    dvdt = dv/dt
    ax[0].plot(sweep.t*1e3, ttl, color='k')
    ax[0].plot(sweep.t[peaks_ttl]*1e3, ttl[peaks_ttl],  'x', color='red')
    ax[1].plot(sweep.t*1e3, v_filt, color='lightgrey')
    ax[2].plot(sweep.t*1e3, dvdt, color='k')
    #calculate and make an average
    average_list.append(v_filt)
    average_list_dvdt.append(dvdt)
average = np.mean(average_list, axis=0)
average_dvdt = np.mean(average_list_dvdt, axis=0)
ax[1].plot(sweep.t*1e3, average, color='red')
ax[2].plot(sweep.t*1e3, average_dvdt, color='red')

#%%detect average events per set 
fig,ax = plt.subplots(3,1, sharex=True)     
ax[0].plot(sweep.t*1e3, ttl, color='k')
ax[1].plot(sweep.t*1e3, average, color='k')
ax[2].plot(sweep.t*1e3, average_dvdt, color='k')
ax[2].set_ylim(-1000, 1000)
peaks_ttl, nonesne = find_peaks(ttl, 1) 
t=sweep.t*1e3
pulse_number = []
onsets = []
rise   = []
qc_list = []
rise_0_100 = []
rise_10_90 = []
amplitude  = []
for j in range(0,len(peaks_ttl)):
    left_border = peaks_ttl[j]+1000
    right_border =  peaks_ttl[j]+5000
    right_border2 = right_border+20000
    aoi = average_dvdt[left_border:right_border]
    aoi2_v = average[left_border:(right_border+20000)]
    aoi2_t = t[left_border:right_border2]
    threshold = int(input("Please provide derravitve threshold for event detection: "))
    new_onset = 'f'
    if threshold > max(aoi):
        print('dvdt does not exceed threshold, please change thershold')
        print('max aoi is do not exceed this value again', max(aoi))
        while new_onset == 'f':
            threshold = int(input("Please provide derravitve threshold for event detection: "))
            onset_index = np.where(aoi>threshold)[0][0] 
            ax[1].plot(t[left_border:right_border][onset_index], average[left_border:right_border][onset_index], 'x')
            peak_temp, nonsense = find_peaks(aoi2_v, prominence=.3)
            ax[1].plot(t[left_border:right_border2][peak_temp], average[left_border:right_border2][peak_temp], 'x')
            #plot rise times (0-100%)
            #adjust rise time so we take 10-90% of maxium in V not t. 
            amp = average[left_border:right_border2][peak_temp] - average[left_border:right_border][onset_index]
            amp_10 = average[left_border:right_border][onset_index] + 0.1*amp
            amp_90 = average[left_border:right_border2][onset_index] + 0.9*amp 
            #cooridanted of amp_10 - ampg_90
            amp_10_x = np.where(average[left_border:right_border] > amp_10)[0][0]
            amp_90_x = np.where(average[left_border:right_border2] > amp_90)[0][0]
            
            #plot total rise time
            ax[1].plot(t[left_border:right_border2][onset_index:peak_temp[0]], average[left_border:right_border2][onset_index:peak_temp[0]], color = 'r')
            rise0_100  =  t[left_border:right_border2][peak_temp[0]]- t[left_border:right_border2][onset_index]
            """
            code below assumed 10 - 90 % of time not the time 10%90% of amplitude 
            
            #plot rise time (10 - 90%)
            #total_rise = t[left_border:right_border2][onset_index:peak_temp[0]], average[left_border:right_border2][onset_index:peak_temp[0]]
            #index_start = int(.1*len(t[left_border:right_border2][onset_index:peak_temp[0]]))
            #index_end   = int(.9*len(t[left_border:right_border2][onset_index:peak_temp[0]]))
            #ax[1].plot(t[left_border:right_border2][onset_index:peak_temp[0]][index_start:index_end], average[left_border:right_border2][onset_index:peak_temp[0]][index_start:index_end], color = 'b')
            #rise10_90 = t[left_border:right_border2][onset_index:peak_temp[0]][index_end] - t[left_border:right_border2][onset_index:peak_temp[0]][index_start]
            #amp       = average[left_border:right_border2][onset_index:peak_temp[0]][index_end] -  average[left_border:right_border2][onset_index:peak_temp[0]][index_start]
            """
            
            #plot 
            ax[1].plot(t[left_border:right_border2][amp_10_x:amp_90_x], average[left_border:right_border2][amp_10_x:amp_90_x], color = 'b')
            rise10_90 = t[left_border:right_border2][amp_90_x] - t[left_border:right_border2][amp_10_x]
            
            print('printing new onset, press a button to continue')
            while True:
                if plt.waitforbuttonpress():
                    break 
            new_onset = input('New onset good? (p/f)')
            qc = input('Does the event qc pass? ')
            qc_list.append(qc)
            onsets.append(t[left_border:right_border][onset_index]-t[peaks_ttl[j]])
            pulse_number.append(j)
            rise_0_100.append(rise0_100)
            rise_10_90.append(rise10_90)
            amplitude.append(amp)
            print('analyzing event, press a button to continue')
    elif threshold < max(aoi):
        onset = 'f'
        while onset == 'f':
            onset_index = np.where(aoi>threshold)[0][0] 
            ax[1].plot(t[left_border:right_border][onset_index], average[left_border:right_border][onset_index], 'x')
            onset_index = np.where(aoi>threshold)[0][0] 
            ax[1].plot(t[left_border:right_border][onset_index], average[left_border:right_border][onset_index], 'x')
            peak_temp, nonsense = find_peaks(aoi2_v, prominence=.3)
            ax[1].plot(t[left_border:right_border2][peak_temp], average[left_border:right_border2][peak_temp], 'x')
            amp = average[left_border:right_border2][peak_temp] - average[left_border:right_border][onset_index]
            amp_10 = average[left_border:right_border][onset_index] + 0.1*amp
            amp_90 = average[left_border:right_border2][onset_index] + 0.9*amp 
            #cooridanted of amp_10 - ampg_90
            amp_10_x = np.where(average[left_border:right_border] > amp_10)[0][0]
            amp_90_x = np.where(average[left_border:right_border2] > amp_90)[0][0]
            
            #plot total rise time
            ax[1].plot(t[left_border:right_border2][onset_index:peak_temp[0]], average[left_border:right_border2][onset_index:peak_temp[0]], color = 'r')
            rise0_100  =  t[left_border:right_border2][peak_temp[0]]- t[left_border:right_border2][onset_index]
            """
            code below assumed 10 - 90 % of time not the time 10%90% of amplitude 
            
            #plot rise time (10 - 90%)
            #total_rise = t[left_border:right_border2][onset_index:peak_temp[0]], average[left_border:right_border2][onset_index:peak_temp[0]]
            #index_start = int(.1*len(t[left_border:right_border2][onset_index:peak_temp[0]]))
            #index_end   = int(.9*len(t[left_border:right_border2][onset_index:peak_temp[0]]))
            #ax[1].plot(t[left_border:right_border2][onset_index:peak_temp[0]][index_start:index_end], average[left_border:right_border2][onset_index:peak_temp[0]][index_start:index_end], color = 'b')
            #rise10_90 = t[left_border:right_border2][onset_index:peak_temp[0]][index_end] - t[left_border:right_border2][onset_index:peak_temp[0]][index_start]
            #amp       = average[left_border:right_border2][onset_index:peak_temp[0]][index_end] -  average[left_border:right_border2][onset_index:peak_temp[0]][index_start]
            """
            
            #plot 
            ax[1].plot(t[left_border:right_border2][amp_10_x:amp_90_x], average[left_border:right_border2][amp_10_x:amp_90_x], color = 'b')
            rise10_90 = t[left_border:right_border2][amp_90_x] - t[left_border:right_border2][amp_10_x]
            while True:
                if plt.waitforbuttonpress():
                    break 
            onset = input('Onset good? (p/f)')
            if onset == 'f':
                threshold = int(input("Please provide derravitve threshold for event detection: "))
        qc = input('Does the event qc pass? ')
        onsets.append(t[left_border:right_border][onset_index]-t[peaks_ttl[j]])
        pulse_number.append(j)
        rise_0_100.append(rise0_100)
        rise_10_90.append(rise10_90)
        qc_list.append(qc)
        amplitude.append(amp)
onsets_df = pd.DataFrame(onsets)
onsets_df['condition'] = condition
onsets_df['file_name'] = f.filename
onsets_df['species']   = species
onsets_df['n_pulse']   = pulse_number
onsets_df['qc']        = qc_list
onsets_df['rise']      = rise_0_100
onsets_df['rise%10_90']= rise_10_90
onsets_df['amplitude'] = amplitude

if species == 'human':
    save_name = f.filename[82:101]+'.csv'
    onsets_df.to_csv(r'C:\Users\sdr267\Documents\PhD\ProjectSynapticConnections\Ih_experiment\events\\' + save_name)
elif species == 'mouse':
    save_name = f.filename[77:102]+'.csv'
    onsets_df.to_csv(r'C:\Users\sdr267\Documents\PhD\ProjectSynapticConnections\Ih_experiment\events\\' + save_name)

#%%

from  get_PSP_param import get_synaptic_params