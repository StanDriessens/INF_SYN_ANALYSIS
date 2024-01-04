# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 15:14:48 2023

@author: sdr267
"""

def get_synaptic_params(average,sweep,average_dvdt, ttl, condition,f, species):
    import matplotlib.pyplot as plt
    from scipy.signal import find_peaks, peak_prominences
    import pandas as pd
    import numpy as np
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
                ax[1].plot(t[left_border:right_border2][onset_index:peak_temp[0]], average[left_border:right_border2][onset_index:peak_temp[0]], color = 'r')
                rise0_100  =  t[left_border:right_border2][peak_temp[0]]- t[left_border:right_border2][onset_index]
                #plot rise time (10 - 90%)
                total_rise = t[left_border:right_border2][onset_index:peak_temp[0]], average[left_border:right_border2][onset_index:peak_temp[0]]
                index_start = int(.1*len(t[left_border:right_border2][onset_index:peak_temp[0]]))
                index_end   = int(.9*len(t[left_border:right_border2][onset_index:peak_temp[0]]))
                ax[1].plot(t[left_border:right_border2][onset_index:peak_temp[0]][index_start:index_end], average[left_border:right_border2][onset_index:peak_temp[0]][index_start:index_end], color = 'b')
                rise10_90 = t[left_border:right_border2][onset_index:peak_temp[0]][index_end] - t[left_border:right_border2][onset_index:peak_temp[0]][index_start]
                amp       = average[left_border:right_border2][onset_index:peak_temp[0]][index_end] -  average[left_border:right_border2][onset_index:peak_temp[0]][index_start]
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
                #plot rise times (0-100%)
                ax[1].plot(t[left_border:right_border2][onset_index:peak_temp[0]], average[left_border:right_border2][onset_index:peak_temp[0]], color = 'r')
                rise0_100  =  t[left_border:right_border2][peak_temp[0]]- t[left_border:right_border2][onset_index]
                #plot rise time (10 - 90%)
                total_rise = t[left_border:right_border2][onset_index:peak_temp[0]], average[left_border:right_border2][onset_index:peak_temp[0]]
                index_start = int(.1*len(t[left_border:right_border2][onset_index:peak_temp[0]]))
                index_end   = int(.9*len(t[left_border:right_border2][onset_index:peak_temp[0]]))
                ax[1].plot(t[left_border:right_border2][onset_index:peak_temp[0]][index_start:index_end], average[left_border:right_border2][onset_index:peak_temp[0]][index_start:index_end], color = 'b')
                rise10_90 = t[left_border:right_border2][onset_index:peak_temp[0]][index_end] - t[left_border:right_border2][onset_index:peak_temp[0]][index_start]
                amp       = average[left_border:right_border2][onset_index:peak_temp[0]][index_end] -  average[left_border:right_border2][onset_index:peak_temp[0]][index_start]
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
    return onsets_df
    
