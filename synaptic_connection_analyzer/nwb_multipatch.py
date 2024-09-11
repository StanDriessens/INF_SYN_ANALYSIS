# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 09:42:23 2024

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
from inspect_select_traces import inspect_select_traces

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
    print('appended stimulus to channel 1', np.unique(stimulus_code_1))
    
    # Stimulus codes per channel 2
    stimulus_code_2 = []
    for sweep_name in f_stimulus_channel2.sweeps:
        stimulus_description = f['stimulus']['presentation'][sweep_name].attrs['stimulus_description']
        stimulus_code_2.append(str(stimulus_description))
    print('appended stimulus to channel 2', np.unique(stimulus_code_2))

    
    # Stimulus codes per channel 3
    stimulus_code_3 = []
    for sweep_name in f_stimulus_channel3.sweeps:
        stimulus_description = f['stimulus']['presentation'][sweep_name].attrs['stimulus_description']
        stimulus_code_3.append(str(stimulus_description))
    print('appended stimulus to channel 3', np.unique(stimulus_code_3))

    
    # Stimulus codes per channel 4
    stimulus_code_4 = []
    for sweep_name in f_stimulus_channel4.sweeps:
        stimulus_description = f['stimulus']['presentation'][sweep_name].attrs['stimulus_description']
        stimulus_code_4.append(str(stimulus_description))
    print('appended stimulus to channel 4', np.unique(stimulus_code_4))

        
   
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
    for index in range(len(stimuli)):
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
    for index in range(len(stimuli)):
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
    
    #include exclude sweeps 
    
    sweeps_inc = inspect_select_traces(f, sweep_tab_pre, sweep_tab_post)
    
    
        
        
    v_pre_tot = []
    v_post_tot = []    
    fig,ax = plt.subplots(3,1, sharex=True)
    for i in sweeps_inc: 
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
            # ax[0].plot(t, v_pre, color='w')
            # ax[1].plot(t, v_post, color='w')
            #append arrays for averaging 
            v_pre_tot.append(v_pre)
            v_post_tot.append(v_post)
            #create average 
    average_pre = np.mean(v_pre_tot, axis=0)
    average_post = np.mean(v_post_tot, axis=0)
    
    #filter the average 
    
    average_post = savgol_filter(average_post, 251 ,2)
    #plot average 
    ax[0].plot(t, average_pre)
    ax[1].plot(t, average_post)
    dv = np.gradient(average_post, edge_order = 1)
    dt = np.gradient(t, edge_order=1)
    dvdt = dv/dt
    ax[2].plot(t, dvdt)
    #here get the requiered data for PSP feature extraction 
    #derravitive for thersholding
 
    #spikes for onset
    ext = SpikeFeatureExtractor()
    I = t 
    results = ext.process(t, average_pre, I)
    #run the PSP analyzers 
    if connection == 'inh':
        time_onset = get_epsp_parameters_spikes_inh(t, average_post, results, dvdt, ax, fig)
    elif connection =='exc':
        time_onset = get_epsp_parameters_spikes_exc(t, average_post, results, dvdt, ax, fig)
        
        
    #save the data 
    #file     
    window10 = Tk()
    save_folder = filedialog.askdirectory()
    path = file
    filename = os.path.basename(path)
    time_onset.to_csv(save_folder + '/' + filename + '_' + 'analyzed' + '.csv')

        
    return time_onset