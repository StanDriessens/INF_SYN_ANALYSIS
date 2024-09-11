# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 10:19:52 2024

@author: sdr267
"""
import pandas as pd
import numpy as np
import h5py

def load_from_h5py(file, channels_to_analyze):
    f =  h5py.File(file, 'r')
    f_sweeps =  pd.DataFrame(list(f['acquisition'].keys()), columns=['sweeps'])
    #get the channels from the sweeps 
    channels = set()
    for sweep in f_sweeps['sweeps']:
        channel = sweep.split('_')[-1]
        channels.add(channel)
    #sweeps for different channels 
    #select one channel for analyses 
    
    channel_1 = list(channels)[0]
    channel_2 = list(channels)[1]
    f_sweeps_channel1= f_sweeps[f_sweeps['sweeps'].str.contains('AD2')].reset_index(drop=True)
    f_sweeps_channel2 = f_sweeps[f_sweeps['sweeps'].str.contains(channel_2)].reset_index(drop=True)
    f_stimulus = pd.DataFrame(list(f['stimulus']['presentation'].keys()), columns = ['sweeps']) 
    f_stimulus_channel1 = f_stimulus[f_stimulus['sweeps'].str.contains('DA2')].reset_index(drop=True)
    stimulus_code = []
    for i in f_stimulus_channel1.sweeps:
        temp = i 
        stimulus_code.append(str(f['stimulus']['presentation'][i].attrs['stimulus_description']))
    #get sweep nubmers
    sweep_n = np.unique(list(f['general']['intracellular_ephys']['sweep_table']['sweep_number']))  
    #create the sweep table 
    sweep_table = pd.DataFrame(columns=['sweep_file_name','stim_file_name', 'stimulus_code', 'sweep_number'])
    sweep_table.sweep_file_name      = f_sweeps_channel1.sweeps
    sweep_table.stim_file_name       = f_stimulus_channel1.sweeps
    sweep_table.stimulus_code        = stimulus_code
    sweep_table.sweep_number         = sweep_n
    #sweep_table = sweep_table.reset_index()
    sweep_table['stimulus_code'] = sweep_table['stimulus_code'].str[2:-1]     
    return  f, sweep_table