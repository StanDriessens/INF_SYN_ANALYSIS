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
#%%          
def multipatch_nwb_analyzer(file, connection): 
    #load the file and make the sweep table 
    f = h5py.File(file, 'r')
    
    f_sweeps =  pd.DataFrame(list(f['acquisition'].keys()), columns=['sweeps'])
    #get the channels from the sweeps 
    # Initialize an empty list to collect channel names
    channels = []
    
    # Extract channel names from sweep names and add them to the list
    for sweep in f_sweeps['sweeps']:
        channel = sweep.split('_')[-1]
        if channel not in channels:
            channels.append(channel)
    
    # Sort the list of channel names
    channels.sort()
    channel_1 = list(channels)[0]
    channel_2 = list(channels)[1]
    channel_3 = list(channels)[2]
    channel_4 = list(channels)[3]
    f_sweeps_channel1 = f_sweeps[f_sweeps['sweeps'].str.contains(channel_1)].reset_index(drop=True)
    f_sweeps_channel2 = f_sweeps[f_sweeps['sweeps'].str.contains(channel_2)].reset_index(drop=True)
    f_sweeps_channel3 = f_sweeps[f_sweeps['sweeps'].str.contains(channel_3)].reset_index(drop=True)
    f_sweeps_channel4 = f_sweeps[f_sweeps['sweeps'].str.contains(channel_4)].reset_index(drop=True)
    
    #f_sweeps_channel1= f_sweeps[f_sweeps['sweeps'].str.contains('AD2')].reset_index(drop=True)
    #f_sweeps_channel2 = f_sweeps[f_sweeps['sweeps'].str.contains(channel_2)].reset_index(drop=True)
    f_stimulus = pd.DataFrame(list(f['stimulus']['presentation'].keys()), columns=['sweeps'])
    f_stimulus_channel1 = f_stimulus[f_stimulus['sweeps'].str.contains('DA0')].reset_index(drop=True)
    f_stimulus_channel2 = f_stimulus[f_stimulus['sweeps'].str.contains('DA1')].reset_index(drop=True)
    f_stimulus_channel3 = f_stimulus[f_stimulus['sweeps'].str.contains('DA2')].reset_index(drop=True)
    f_stimulus_channel4 = f_stimulus[f_stimulus['sweeps'].str.contains('DA3')].reset_index(drop=True)
    
    # Initialize lists to store stimulus codes
    stimulus_code_1 = []
    stimulus_code_2 = []
    stimulus_code_3 = []
    stimulus_code_4 = []
    
    # Retrieve stimulus codes for each channel
    for i in f_stimulus_channel1.sweeps:
        stimulus_code_1.append(str(f['stimulus']['presentation'][i].attrs['stimulus_description']))
        
    for i in f_stimulus_channel2.sweeps:
        stimulus_code_2.append(str(f['stimulus']['presentation'][i].attrs['stimulus_description']))
        
    for i in f_stimulus_channel3.sweeps:
        stimulus_code_3.append(str(f['stimulus']['presentation'][i].attrs['stimulus_description']))
        
    for i in f_stimulus_channel4.sweeps:
        stimulus_code_4.append(str(f['stimulus']['presentation'][i].attrs['stimulus_description']))
    
    
       
        # Stimulus codes per channel 1
    stimulus_code_1 = []
    for sweep_name in f_stimulus_channel1.sweeps:
        stimulus_description = f['stimulus']['presentation'][sweep_name].attrs['stimulus_description']
        stimulus_code_1.append(str(stimulus_description))
    
    # Stimulus codes per channel 2
    stimulus_code_2 = []
    for sweep_name in f_stimulus_channel2.sweeps:
        stimulus_description = f['stimulus']['presentation'][sweep_name].attrs['stimulus_description']
        stimulus_code_2.append(str(stimulus_description))
    
    # Stimulus codes per channel 3
    stimulus_code_3 = []
    for sweep_name in f_stimulus_channel3.sweeps:
        stimulus_description = f['stimulus']['presentation'][sweep_name].attrs['stimulus_description']
        stimulus_code_3.append(str(stimulus_description))
    
    # Stimulus codes per channel 4
    stimulus_code_4 = []
    for sweep_name in f_stimulus_channel4.sweeps:
        stimulus_description = f['stimulus']['presentation'][sweep_name].attrs['stimulus_description']
        stimulus_code_4.append(str(stimulus_description))
        
    # create the sweep table
    # Channel 1
    sweep_table_channel1 = pd.DataFrame({
       'sweep_file_name': f_sweeps_channel1.sweeps,
       'stim_file_name': f_stimulus_channel1.sweeps,
       'stimulus_code': stimulus_code_1,
       'channel': 'channel_1'
    })
    
    # Channel 2
    sweep_table_channel2 = pd.DataFrame({
       'sweep_file_name': f_sweeps_channel2.sweeps,
       'stim_file_name': f_stimulus_channel2.sweeps,
       'stimulus_code': stimulus_code_2,
       'channel': 'channel_2'
    })
    
    # Channel 3
    sweep_table_channel3 = pd.DataFrame({
       'sweep_file_name': f_sweeps_channel3.sweeps,
       'stim_file_name': f_stimulus_channel3.sweeps,
       'stimulus_code': stimulus_code_3,
       'channel': 'channel_3'
    })
    
    # Channel 4
    sweep_table_channel4 = pd.DataFrame({
       'sweep_file_name': f_sweeps_channel4.sweeps,
       'stim_file_name': f_stimulus_channel4.sweeps,
       'stimulus_code': stimulus_code_4,
       'channel': 'channel_4'
    })
    
    sweep_table = pd.concat([sweep_table_channel1, sweep_table_channel2, sweep_table_channel3, sweep_table_channel4], ignore_index=True)
    
    # Optionally, if you want to convert sweep_file_name and stim_file_name to string type
    sweep_table['stim_file_name'] = sweep_table['stim_file_name'].astype(str)
    
    #get sweep nubmers 
    pattern = r'data_(\d{5})_AD*'
    # Extract the sweep numbers
    sweep_table['sweep_number'] = sweep_table['sweep_file_name'].str.extract(pattern)
    
    """
    prompt the user for channel and stimulus 
    """
    #select pre channel
    channels = list(np.unique(sweep_table.channel))
    def select_pre_channel():
        global pre_channel
        pre_channel = channels[x.get()]
        window6.destroy()
    window6 = Tk()
    x = IntVar()
    for index in range(len(channels)):
        radiobutton = Radiobutton(window6, text=channels[index], variable=x, value=index, padx=25)
        radiobutton.pack(anchor=W)
    okay_button = Button(window6, text="Select pre channel", command=select_pre_channel)
    okay_button.pack(pady=10)    
    window6.mainloop()  
    #select post channel 
    def select_post_channel():
        global post_channel
        post_channel = channels[x.get()]
        window7.destroy()
    window7 = Tk()
    x = IntVar()
    for index in range(len(channels)):
        radiobutton = Radiobutton(window7, text=channels[index], variable=x, value=index, padx=25)
        radiobutton.pack(anchor=W)
    okay_button = Button(window7, text="Select post channel", command=select_post_channel)
    okay_button.pack(pady=10)    
    window7.mainloop()  
    
    #get the list of stimulus codes
    stimuli = list(np.unique(sweep_table.stimulus_code))
    def select_pre_stimulus():
        global pre_stimulus
        pre_stimulus = stimuli[x.get()]
        window8.destroy()
    window8 = Tk()
    x = IntVar()
    for index in range(len(channels)):
        radiobutton = Radiobutton(window8, text=stimuli[index], variable=x, value=index, padx=25)
        radiobutton.pack(anchor=W)
    okay_button = Button(window8, text="Select pre stimulus", command=select_pre_stimulus)
    okay_button.pack(pady=10)    
    window8.mainloop()  
    
    def select_post_stimulus():
        global post_stimulus
        post_stimulus = stimuli[x.get()]
        window9.destroy()
    window9 = Tk()
    x = IntVar()
    for index in range(len(channels)):
        radiobutton = Radiobutton(window9, text=stimuli[index], variable=x, value=index, padx=25)
        radiobutton.pack(anchor=W)
    okay_button = Button(window9, text="Select post stimulus", command=select_post_stimulus)
    okay_button.pack(pady=10)    
    window8.mainloop()  
    
    """
    make a selction of sweeps based on user input 
    """
    #make a sweep table for the selected channels and stimulus 
    sweep_tab_pre  = sweep_table[(sweep_table['channel'] == pre_channel) & (sweep_table['stimulus_code'] == pre_stimulus)]
    sweep_tab_post = sweep_table[(sweep_table['channel'] == post_channel) & (sweep_table['stimulus_code'] == post_stimulus)]
    #get the correct sweep_number from the pre_sweep 
    pre_sweeps = list(sweep_tab_pre.sweep_number)
    #only select those sweeps from the post_sweeps
    sweep_tab_post = sweep_tab_post[sweep_tab_post['sweep_number'].isin(pre_sweeps)]
    #now the other way around
    post_sweeps = list(sweep_tab_post.sweep_number)
    sweep_tab_pre =sweep_tab_pre[sweep_tab_pre['sweep_number'].isin(post_sweeps)]
    #error check if sweep numbers match 
    check = set( sweep_tab_pre.sweep_number) ==  set(sweep_tab_post.sweep_number)
    if check ==True:
        print('sweeps foud analyzing')
        
        
    v_pre_tot = []
    v_post_tot = []    
    fig,ax = plt.subplots(4,1, sharex=True)
    for i in sweep_tab_pre.sweep_number: 
        #get pre sweep
        key_pre = sweep_tab_pre[sweep_tab_pre['sweep_number'] == i].sweep_file_name.reset_index()
        v_pre = np.array(f['acquisition'][key_pre.sweep_file_name[0]]['data'])
        #get_poswet sweep 
        key_post = sweep_tab_post[sweep_tab_post['sweep_number'] == i].sweep_file_name.reset_index()
        v_post = np.array(f['acquisition'][key_post.sweep_file_name[0]]['data'])
        
        """
        allow space for filering of the signal ?? 
        v_pre_filt = .... 
        
        """
        # to consturct the time array we need
        # sample rate 
        # length of the array 
        
        t = np.arange(0, len(f['acquisition'][key_pre.sweep_file_name[0]]['data'])) / f['acquisition'][key_pre.sweep_file_name[0]]['starting_time'].attrs['rate']
        #use ipffx to detect spikes in the sweeps 
        if max(v_pre) < 100:
            ax[0].plot(t, v_pre, color='w')
            ax[1].plot(t, v_post, color='w')
            #append arrays for averaging 
            v_pre_tot.append(v_pre)
            v_post_tot.append(v_post)
            #create average 
    average_pre = np.mean(v_pre_tot, axis=0)
    average_post = np.mean(v_post_tot, axis=0)
    #plot average 
    ax[2].plot(t, average_pre)
    ax[3].plot(t, average_post)
        
            
        
    
 
    
        
   
     
     