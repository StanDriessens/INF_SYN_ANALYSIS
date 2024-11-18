# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 13:32:48 2024

@author: sdr267
"""
#import packages
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
from run_cc_steps import SpikeFeats, InputResistance, RunSag

def run_nwb_SpikeFeats(file, species):
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
    print('Channels found in file: ', channels)
    # Sort the list of channel names
    
    
    channels.sort()

    # Dictionary to store sweeps and stimulus information
    f_sweeps_channels = {}
    stimulus_codes = {}
    sweep_tables = []
    
    # Get stimulus information
    f_stimulus = pd.DataFrame(list(f['stimulus']['presentation'].keys()), columns=['sweeps'])
    stimulus_channels = ['DA0', 'DA1', 'DA2', 'DA3']
    f_stimulus_channels = {ch: f_stimulus[f_stimulus['sweeps'].str.contains(ch)].reset_index(drop=True) for ch in stimulus_channels}
    
    # Loop over available channels and populate the sweep and stimulus information
    for idx, channel in enumerate(channels):
        # Extract sweeps for this channel
        f_sweeps_channel = f_sweeps[f_sweeps['sweeps'].str.contains(channel)].reset_index(drop=True)
        f_sweeps_channels[channel] = f_sweeps_channel
        
        # Find the corresponding stimulus channel (e.g., DA0, DA1, ...)
        stim_channel_name = f'DA{idx}'
        f_stimulus_channel = f_stimulus_channels.get(stim_channel_name, pd.DataFrame(columns=['sweeps']))
    
        # Retrieve stimulus codes for each sweep if the stimulus channel exists
        stimulus_code = []
        if not f_stimulus_channel.empty:
            for sweep_name in f_stimulus_channel.sweeps:
                stimulus_description = f['stimulus']['presentation'][sweep_name].attrs['stimulus_description']
                stimulus_code.append(str(stimulus_description))
            stimulus_codes[stim_channel_name] = stimulus_code
    
        # Only create the DataFrame if both sweep and stimulus data are available
        if not f_sweeps_channel.empty and stimulus_code:
            sweep_table_channel = pd.DataFrame({
                'sweep_file_name': f_sweeps_channel.sweeps,
                'stim_file_name': f_stimulus_channel.sweeps,
                'stimulus_code': stimulus_code,
                'channel': f'channel_{idx + 1}'
            })
            sweep_tables.append(sweep_table_channel)
            print(f'Appended stimulus to {stim_channel_name}:', np.unique(stimulus_code))
    
    # Combine all sweep tables into a single DataFrame
    sweep_table = pd.concat(sweep_tables, ignore_index=True) if sweep_tables else pd.DataFrame()
    
    #make radio buttn with protocols 
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
    okay_button = Button(window6, text="Select channel", command=select_pre_channel)
    okay_button.pack(pady=10)    
    window6.mainloop()  
    
    sweep_table_channel = sweep_table[sweep_table['channel'] == pre_channel]
    stimuli = list(np.unique(sweep_table.stimulus_code))
    def select_stim_set():
        global StimSet
        StimSet = stimuli[x.get()]
        window8.destroy()
    window8 = Tk()
    x = IntVar()
    for index in range(len(stimuli)):
        radiobutton = Radiobutton(window8, text=stimuli[index], variable=x, value=index, padx=25)
        radiobutton.pack(anchor=W)
    okay_button = Button(window8, text="Select StimSet", command=select_stim_set)
    okay_button.pack(pady=10)    
    window8.mainloop()  
    
    
    
    if 'tep' in StimSet :
        print('Stimset Selected:', StimSet)
        sweep_tab_pre  = sweep_table[(sweep_table['channel'] == pre_channel) & (sweep_table['stimulus_code'] == StimSet)]
        print(sweep_tab_pre)
        results   = SpikeFeats(f, sweep_tab_pre)
        Rin       = InputResistance(f, sweep_tab_pre)
        #get sag df 
        SagDf     = RunSag(f, file, sweep_tab_pre, 1000)
        #calcualte one sag feature at the highest current injection 
        SagFeat = SagDf[SagDf['StimAmp'] == min(SagDf.StimAmp)].reset_index()
        #make sure one is selected 
        Sag=  (SagFeat.iloc[0].MinVoltage - SagFeat.iloc[0].SteadState)/ SagFeat.iloc[0].VoltageDeflect
        #make and save a dataframe  
        results['Rin']  =  Rin
        results['Sag']  =  Sag
        results['cahnnel'] = pre_channel
        SagDf['channel']   = pre_channel 
        path = file
        filename = os.path.basename(path)
        results['cell'] =   filename 
        save_folder = filedialog.askdirectory()
        #save to csv
        results.to_csv(save_folder + '/' + filename + '_' + 'analyzed' + pre_channel + '.csv')
        SagDf.to_csv(save_folder + '/' + filename + '_' + 'Sag' + pre_channel + '.csv')
    else:
        print('No CC steps found in this channel')
    return  results, SagDf
        

        
    
    