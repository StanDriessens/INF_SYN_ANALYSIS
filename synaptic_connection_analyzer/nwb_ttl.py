# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 16:40:44 2024

This function opens and analyzes TTL to post synaptic responses 

Ususally this is called from intialize.py 

For ABF spike analysis see abf_spike.py 

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
from tkinter import simpledialog
from get_epsp_parameters_abf_new import get_epsp_parameters_ttl_exc
from get_epsp_parameters_abf_new import get_epsp_parameters_spikes_inh
from get_epsp_parameters_abf_new import get_epsp_parameters_spikes_exc 
from scipy.signal import savgol_filter


def run_nwb_ttl_analyzer(file, species):
    #load in the dataset 
    print('loading dataset')
    dataset= create_ephys_data_set(file)
    all_sweeps = dataset.filtered_sweep_table()
    
    #get the correct sweep name from tkinter button
    sets_unique = list(np.unique(all_sweeps.stimulus_code))
    def save_selected_value():
        global cond
        cond = sets_unique[x.get()]
        window2.destroy()
    window2 = Tk()
    x = IntVar()
    for index in range(len(sets_unique)):
        radiobutton = Radiobutton(window2, text=sets_unique[index], variable=x, value=index, padx = 25)
        radiobutton.pack(anchor=W)
    okay_button = Button(window2, text="Select sweep protocol", command=save_selected_value)
    okay_button.pack(pady=10)    

    window2.mainloop()

    print(cond)
    
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
    
    if 'ZD' not in cond:
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
                    
    if 'ZD' in cond:
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
    # else:
    #     print('analysis has stopped no correct protocol selected')
    #     print('do you want to analyze sweeps manually? ')
    #     sweeps_to_analyze = list(simpledialog.askstring('Sweeps', 'Please enter sweep numbers to analyze in list format [21, 23 , etc]'))
        
    
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
        ax[0].plot(sweep.t*1e3, ttl, color='w')
        ax[0].plot(sweep.t[peaks_ttl]*1e3, ttl[peaks_ttl],  'x', color='red')
        ax[1].plot(sweep.t*1e3, v_filt, color='lightgrey')
        ax[2].plot(sweep.t*1e3, dvdt, color='w')

        #calculate and make an average
        average_list.append(v_filt)
        average_list_dvdt.append(dvdt)
    average = np.mean(average_list, axis=0)
    average_dvdt = np.mean(average_list_dvdt, axis=0)
    ax[1].plot(sweep.t*1e3, average, color='red')
    ax[2].plot(sweep.t*1e3, average_dvdt, color='red')
    tkinter.messagebox.showinfo('average tracing', message='press OK to continue event analysis')
    #detect average events per set 
    fig,ax = plt.subplots(3,1, sharex=True)
    fig.set_size_inches(12,6)     
    ax[0].plot(sweep.t, ttl, color='w')
    ax[1].plot(sweep.t, average, color='w')
    ax[2].plot(sweep.t, average_dvdt, color='w')
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
    
    post_syn_x = sweep.t
    post_syn_y = average
    
    time_onset = get_epsp_parameters_ttl_exc(post_syn_x=post_syn_x, post_syn_y=post_syn_y, peaks_ttl=peaks_ttl,
                                         average_dvdt=average_dvdt, ax=ax, fig=fig)

    if species == 'human':
        save_name = f.filename[82:101]+'.csv'
        time_onset.to_csv(r'C:\Users\sdr267\Documents\PhD\ProjectSynapticConnections\Ih_experiment\events\\' + save_name)
    elif species == 'mouse':
         save_name = f.filename[77:102]+'.csv'
         time_onset.to_csv(r'C:\Users\sdr267\Documents\PhD\ProjectSynapticConnections\Ih_experiment\events\\' + save_name)
         

    
    
    
    
    
    
    
        
            
        
    
 
    
        
   
     
     