# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 13:55:22 2024

@author: sdr267
"""

import matplotlib.pyplot as plt
import numpy as np
from ipfx.dataset.create import create_ephys_data_set
from ipfx.feature_extractor import SpikeFeatureExtractor
import pandas as pd 
import scipy as sp
from scipy.signal import find_peaks
import os 

def SpikeFeats(f, sweep_table):
    fig, ax = plt.subplots(2,1, sharex=True)
    df = pd.DataFrame(columns=['spikes', 'sweep', 'stim'])
    for i,j in zip(sweep_table.sweep_file_name, sweep_table.stim_file_name):
        v = np.array(f['acquisition'][i]['data'])
        I = np.array(f['stimulus']['presentation'][j]['data'])
        t = np.arange(0, len(f['acquisition'][i]['data'])) / f['acquisition'][i]['starting_time'].attrs['rate']
        #plot the data 
        ax[0].plot(t,v)
        ax[1].plot(t,I)
        #run results to find first sweep after rheobase
        if str(v[0])!='nan':
            ext = SpikeFeatureExtractor()
            results = ext.process(t,v,I)
            if results.empty:
                # Create a consistent row with None for other parameters
                print('no spikes found in sweep')
                results = pd.DataFrame({
                    'spikes': [0],
                    'sweep': [i],
                    'stim': [j]
                })
            else:   
                results['spikes'] = len(results)
                results['sweep']  = i
                results['stim']   = j   
            df = pd.concat([df,results],axis = 0)
        else:
            print('empty array found skipping')
    return df 

def InputResistance(f, sweep_table, sharex=True):
    fig, ax = plt.subplots(2,1)  
    V = []
    I2 = []
    plt.figure()
    for i,j in zip(sweep_table.sweep_file_name, sweep_table.stim_file_name):
        v = np.array(f['acquisition'][i]['data'])
        I = np.array(f['stimulus']['presentation'][j]['data'])
        t = np.arange(0, len(f['acquisition'][i]['data'])) / f['acquisition'][i]['starting_time'].attrs['rate']
        #plot the data 
        if min(I[10000:]) < 0:
            print('Stim = hyperpol')
            ax[0].plot(t,v)
            ax[1].plot(t,I)
    # get the adta for all
            V.append(min(v)*1e-3)
            I2.append(min(I)*1e-12)
            plt.plot(I2,V, 'x')
            plt.xlabel('current injection (A)')
            plt.ylabel('membrane voltage (V)')
    slope, intercept = np.polyfit(I2, V, 1)
    #plot  the fitted line 
    x_range = np.linspace(min(I2), max(I2), len(I2))
    y_range = slope * x_range + intercept
    plt.plot(x_range, y_range)
    Rin = slope/1e6
    #convert ohm to MegaOhm 
    return Rin
    
    
def RunSag(f, file, sweep_table, cutoff):
    df = pd.DataFrame()
    Sweep              = []
    StimAmp            = []
    MinVoltage         = []
    VoltageDeflect     = []
    SteadState         = []
    BaselineVolt       = []
    df  = pd.DataFrame(columns=['Sweep', 'StimAmp', 'MinVoltage', 'VoltageDeflect', 'SteadState','BaselineVolt'])
    fig, ax = plt.subplots(2,1, sharex=True)
    for i,j in zip(sweep_table.sweep_file_name, sweep_table.stim_file_name):
        v = np.array(f['acquisition'][i]['data'])
        I = np.array(f['stimulus']['presentation'][j]['data'])
        t = np.arange(0, len(f['acquisition'][i]['data'])) / f['acquisition'][i]['starting_time'].attrs['rate']
        #filter the data 
        #cutoff = 30000 #Hz
        order = 1
        fs    =  f['acquisition'][i]['starting_time'].attrs['rate']
        b, a = sp.signal.butter(order, cutoff / (fs/2), btype='low')
        v_filt = sp.signal.filtfilt(b, a, v)
        #plot the data 
        if min(I[10000:]) < 0:
            print('Stim = hyperpol')
            ax[0].plot(t,v_filt)
            ax[1].plot(t,I)
            #flip signal to get nice peak 
            Iflip = I*-1 
            #if there is a test pulse index peak 1 not zero 
            peaks, properties = find_peaks(Iflip, prominence = 10)
            ax[1].plot(t[properties['left_bases'][1]], I[properties['left_bases'][1]], 'x')
            #find and plot maximum voltage deflection 
            MinV = min(v_filt[properties['left_bases'][1]:(properties['left_bases'][1] + 8000)]) 
            ax[0].plot(t[np.where(v_filt == MinV)[0][0]], v_filt[np.where(v_filt == MinV)[0][0]], 'x', color='r')
            #find the steady state 
            SteadyV = v_filt[properties['right_bases'][1]-250]
            ax[0].plot(t[properties['right_bases'][1]], v[properties['right_bases'][1]], 'x', color='b')
            #get sag parameters 
            Sweep.append(i)
            StimAmp.append(min(I))
            BaselineVolt.append(v_filt[properties['left_bases'][1]-100])
            MinVoltage.append(MinV)
            VoltageDeflect.append(v_filt[properties['left_bases'][1]-100] - MinV)
            SteadState.append(SteadyV)
    df['Sweep']          = Sweep
    df['StimAmp']        = StimAmp 
    df['MinVoltage']     = MinVoltage  
    df['VoltageDeflect'] = VoltageDeflect
    df['SteadState']     = SteadState
    df['BaselineVolt']   = BaselineVolt
    df['cell']           = os.path.basename(file)
    return df  



  
       