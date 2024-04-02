# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 16:12:51 2024

analyzing connections from abf file from spike to synaptic response

Usually this file is called from intialze.py 


Stan Driessens 2024/01/08

@author: sdr267
"""


from ipfx.dataset.create import create_ephys_data_set
from ipfx.feature_extractor import SpikeFeatureExtractor
import h5py
import matplotlib.pyplot as plt 
from scipy.signal import order_filter
from scipy.signal import butter,filtfilt
import pyabf
from scipy.signal import savgol_filter
import pandas as pd 
import tkinter
from tkinter import *
from tkinter import filedialog
from tkinter import messagebox as mb
import numpy as np
import time
import os 
from get_epsp_parameters_abf_new import get_epsp_parameters_spikes_exc
from get_epsp_parameters_abf_new import get_epsp_parameters_spikes_inh
from sweep_select import sweep_select
from get_channels_abf import make_radio_channels_abf


def run_abf_spike_analyzer(file, species, connection):
    root =tkinter.Tk()
    pre_in = tkinter.simpledialog.askstring('Pre_channel', 'Please provide channel for pre-synaptic neuron', parent=root)
    post_in = tkinter.simpledialog.askstring('Post_channel', 'Please provide channel for pre-synaptic neuron', parent=root)
 
    train  =  tkinter.simpledialog.askinteger("Amount of spikes in train ", 'amount of spikes in train (6)')
    abf = pyabf.ABF(file)
    channel_names = abf.adcNames
    pre_channel = channel_names.index(pre_in) 
    post_channel = channel_names.index(post_in) 
    onsets = pd.DataFrame(columns=['latency'])
    plt.ion()
    for i in range(0,abf.sweepCount):
        abf.setSweep(i,pre_channel)
        print('analyzing sweep {} out of {}'.format(i, abf.sweepCount-1))
        #fitler  the pre_syanptic signal 
        pre_syn_y = savgol_filter(abf.sweepY, 99 ,5)
        pre_syn_x = abf.sweepX
        pre_syn_i = abf.sweepC
        abf.setSweep(i,post_channel)
        #filter the post synaptic signal
        post_syn_y = savgol_filter(abf.sweepY, 251 ,2)
        post_syn_x = abf.sweepX
        #get spike data using the IPFX spike extractor (AIBS)
        ext = SpikeFeatureExtractor()
        results = ext.process(pre_syn_x, pre_syn_y, pre_syn_i)
        #check if the amount of spikes matches the amount of spikes expected in the sweep 
        if len(results) == train:
            time_max_dvdt = results.upstroke_index
            #get index of post synaptic response onset
            dv = np.gradient(post_syn_y)
            dt = np.gradient(post_syn_x)
            dvdt = dv/dt
            fig, ax = plt.subplots(3, 1, sharex=True )
            fig.set_size_inches(12,6) 
            fig.suptitle('Sweep {}'.format(i))

            #plot the pre synaptic spikes
            #and the max dvdt of the spike
            ax[0].plot(pre_syn_x, pre_syn_y, color = 'w')
            ax[0].plot(pre_syn_x[time_max_dvdt], pre_syn_y[time_max_dvdt], 'x')
            ax[0].set_ylabel('pre_synaptic membrane potential (mV)')
            #plot the post synaptic response 
            ax[1].plot(post_syn_x, post_syn_y, color = 'w')
            ax[1].set_ylim(np.mean(post_syn_y) - 1, np.mean(post_syn_y) + 1)
            ax[1].set_ylabel('post_synaptic membrane potential (mV)')
            # plot the post synaptic membrane derravative 
            ax[2].plot(post_syn_x, dvdt, color = 'lime')
            ax[2].set_ylabel('post_synaptic membrane derivative (mV/ms)')
            #get 2*std of the membrane derative for event detection
            if connection == 'exc':
               time_onset = get_epsp_parameters_spikes_exc(post_syn_x, post_syn_y, results,  dvdt, ax, fig)
            elif connection == 'inh':
               time_onset =  get_epsp_parameters_spikes_inh(post_syn_x, post_syn_y, results,  dvdt, ax, fig)
        else:
            print('wrong amount of spikes detected, skipping sweep')
            time_onset = pd.DataFrame()
            time_onset['latency'] = [None] 
            time_onset['amplitude'] = [None]
            time_onset['sweep'] = [i]
            time_onset['qc'] = [0]    
        onsets = pd.concat([onsets, time_onset], axis = 0)
        onsets['cell'] = file
        
        #save the data  
        
        #select save directory 
    path = file    
    save_folder = filedialog.askdirectory()
    save_name = file[99:109] + '_' +  file[110:113] + '_' + pre_in + '_' + 'to' +'_'+ post_in + '.csv'
    
    save_name = os.path.basename(path) + '_' + pre_in + '_' + post_in
    
    onsets.to_csv(save_folder + '//' + save_name + '.csv')
    
    
"""
to do: 

    make average and run functions in average script. 

"""

def run_abf_spike_analyzer_average(file, species, connection):
    root =tkinter.Tk()
    #get the channels from the file 
    pre_in, post_in = make_radio_channels_abf(file)
    # train  =  tkinter.simpledialog.askinteger("Amount of spikes in train ", 'amount of spikes in train (6)')
    abf = pyabf.ABF(file)
    channel_names = abf.adcNames
    pre_channel = channel_names.index(pre_in) 
    post_channel = channel_names.index(post_in) 
    onsets = pd.DataFrame(columns=['latency'])
    post_sweeps_y = []
    post_sweeps_x = []
    pre_sweeps_y = []
    plt.ion()    
    
    #sweeps selecting
    included_sweeps = sweep_select(abf, pre_channel, post_channel)

    for i in included_sweeps:
        abf.setSweep(i,pre_channel)
        print('analyzing sweep {} out of {}'.format(i, abf.sweepCount-1))
        #fitler  the pre_syanptic signal 
        pre_syn_y = savgol_filter(abf.sweepY, 99 ,5)
        pre_syn_x = abf.sweepX
        pre_syn_i = abf.sweepC
        abf.setSweep(i,post_channel)
        #filter the post synaptic signal
        post_syn_y = savgol_filter(abf.sweepY, 251 ,2)
        post_syn_x = abf.sweepX
        #get spike data using the IPFX spike extractor (AIBS)
        
        #collect and average the sweeps 
        post_sweeps_x.append(post_syn_x)
        post_sweeps_y.append(post_syn_y)
        pre_sweeps_y.append(pre_syn_y)
    post_average_sweep_y = np.mean(post_sweeps_y, axis = 0)
    post_average_sweep_x = np.mean(post_sweeps_x, axis = 0)
    pre_average_sweep_y = np.mean(pre_sweeps_y, axis =0)
    ext = SpikeFeatureExtractor()
    results = ext.process(post_average_sweep_x, pre_average_sweep_y, pre_syn_i)
    #check if the amount of spikes matches the amount of spikes expected in the sweep 
    
    time_max_dvdt = results.upstroke_index
    #get index of post synaptic response onset
    dv = np.gradient(post_average_sweep_y)
    dt = np.gradient(post_average_sweep_x)
    dvdt = dv/dt
    fig, ax = plt.subplots(3, 1, sharex=True )
    fig.set_size_inches(12,6) 
    fig.suptitle('Sweep {}'.format(i))

    #plot the pre synaptic spikes
    #and the max dvdt of the spike
    ax[0].plot(post_average_sweep_x, pre_average_sweep_y, color = 'w')
    ax[0].plot(pre_syn_x[time_max_dvdt], pre_average_sweep_y[time_max_dvdt], 'x')
    ax[0].set_ylabel('pre_synaptic membrane potential (mV)')
    #plot the post synaptic response 
    ax[1].plot(post_average_sweep_x, post_average_sweep_y, color = 'w')
    ax[1].set_ylim(np.mean(post_average_sweep_y) - 1, np.mean(post_average_sweep_y) + 1)
    ax[1].set_ylabel('post_synaptic membrane potential (mV)')
    # plot the post synaptic membrane derravative 
    ax[2].plot(post_average_sweep_x, dvdt, color = 'lime')
    ax[2].set_ylabel('post_synaptic membrane derivative (mV/ms)')
    #get 2*std of the membrane derative for event detection
    if connection == 'exc':
       time_onset = get_epsp_parameters_spikes_exc(post_average_sweep_x, post_average_sweep_y, results,  dvdt, ax, fig)
    elif connection == 'inh':
       time_onset =  get_epsp_parameters_spikes_inh(post_average_sweep_x, post_average_sweep_y, results,  dvdt, ax, fig)
    onsets = pd.concat([onsets, time_onset], axis = 0)
    onsets['cell'] = file
        
    #save the data  
    
    #select save directory 
    path = file    
    save_folder = filedialog.askdirectory()
    save_name = file[99:109] + '_' +  file[110:113] + '_' + pre_in + '_' + 'to' +'_'+ post_in + '_'+ 'average' + '.csv'
    
    save_name = os.path.basename(path) + '_' + pre_in + '_' + post_in
    
    onsets.to_csv(save_folder + '//' + save_name + '.csv')
    
    