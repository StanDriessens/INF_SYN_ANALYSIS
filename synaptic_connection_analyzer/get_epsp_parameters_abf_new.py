# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 10:36:53 2024

function to get synaptic parameters from abf trace

@author: sdr267
"""

import numpy as np
import scipy as sp 
import matplotlib.pyplot as plt
from tkinter import * 
from tkinter import filedialog
from tkinter import messagebox as mb
from tkinter import simpledialog
import pandas as pd 
from fit_decay import fit_decay

def get_epsp_parameters_spikes_exc(post_syn_x, post_syn_y, results,  dvdt, ax, fig):
    stdev = np.std(dvdt)
    ax[2].axhline(y = (np.mean(dvdt) + 3*stdev), color = 'r', linestyle = '--')
    threshold = 3*stdev 
    time_onset = []
    amplitude  = []
    qc         = []
    rise_time_100 = []
    rise_time_90  = []
    root = Tk()
    decay = []
    for j in range(0,len(results)):
        left_border = results.iloc[j].upstroke_index-10
        #onset
        right_border =  results.iloc[j].upstroke_index+400
        right_border2 = right_border + 2000
        right_border3 = int(right_border + 8e3)
        aoi = dvdt[left_border:right_border]
        #find the point where dvdt exceeds 30 mv/ms
        if max(aoi) >= threshold:
            epsp_onset = np.where(aoi > threshold)[0][0]
            amplitude_time = np.where(post_syn_y[left_border:right_border] == max(post_syn_y[left_border:right_border]))[0][0]
            ax[1].plot(post_syn_x[left_border:right_border][epsp_onset], post_syn_y[left_border:right_border][epsp_onset], 'x')
            ax[1].plot(post_syn_x[left_border:right_border][amplitude_time],max(post_syn_y[left_border:right_border]), 'x')
            #get rise times
            amp = post_syn_y[left_border:right_border2][amplitude_time] - post_syn_y[left_border:right_border][epsp_onset]
            amp_10 = post_syn_y[left_border:right_border][epsp_onset] + 0.1*amp
            amp_90 = post_syn_y[left_border:right_border2][epsp_onset] + 0.9*amp 
            #cooridanted of amp_10 - ampg_90
            amp_10_x = np.where(post_syn_y[left_border:right_border] > amp_10)[0][0]
            amp_90_x = np.where(post_syn_y[left_border:right_border2] > amp_90)[0][0]
            ax[1].plot(post_syn_x[left_border:right_border2][epsp_onset:amplitude_time], post_syn_y[left_border:right_border2][epsp_onset:amplitude_time], color = 'r')
            #0-100 rise time
            rise0_100  =  post_syn_x[left_border:right_border2][amplitude_time]- post_syn_x[left_border:right_border2][epsp_onset]
            #plot and get 10 - 90 % rise time 
            ax[1].plot(post_syn_x[left_border:right_border2][amp_10_x:amp_90_x], post_syn_y[left_border:right_border2][amp_10_x:amp_90_x], color = 'cyan')
            rise10_90 = post_syn_x[left_border:right_border2][amp_90_x] - post_syn_y[left_border:right_border2][amp_10_x]
            #plot and get decat decay area
            ax[1].plot(post_syn_x[left_border:right_border3][amplitude_time:right_border3], post_syn_y[left_border:right_border3][amplitude_time:right_border3], color='magenta')
            tau = fit_decay(amplitude_time, left_border, right_border3, post_syn_y, post_syn_x)
            fig.canvas.draw()    
            root.deiconify()
            pass_fail = simpledialog.askstring('threshold detection', 'does the onset pass ? (p/f)', parent=root)
            root.withdraw()
            if pass_fail == 'p':
                time_onset.append((post_syn_x[left_border:right_border][epsp_onset] - results.iloc[j].upstroke_t)*1e3)
                amplitude.append(max(post_syn_y[left_border:right_border])-post_syn_y[left_border:right_border][epsp_onset])
                qc.append('pass')
                rise_time_100.append(rise0_100)
                rise_time_90.append(rise10_90)
                decay.append(tau)

            elif pass_fail == 'f':
               time_onset.append((post_syn_x[left_border:right_border][epsp_onset] - results.iloc[j].upstroke_t)*1e3)
               amplitude.append(max(post_syn_y[left_border:right_border])-post_syn_y[left_border:right_border][epsp_onset])
               qc.append('fail')
               rise_time_100.append(rise0_100)
               rise_time_90.append(rise10_90)
               decay.append(tau)

        else:
            messagebox.showinfo(title='no event', message = 'no event exceeding threshold found')
            time_onset.append(None)
            amplitude.append(None)
            qc.append(None)
            rise_time_100.append(None)
            rise_time_90.append(None)
            

    plt.close()
    Time_Onset = {'latency':time_onset, 'amplitude':amplitude, 'qc':qc, 'rise_time_0_100': rise_time_100, 'rise_time_10_90': rise_time_90, 'tau': decay}   
    time_onset = pd.DataFrame(Time_Onset)
    return time_onset
    
         

def get_epsp_parameters_spikes_inh(post_syn_x, post_syn_y, results,  dvdt, ax, fig):
    stdev = np.std(dvdt)
    ax[2].axhline(y = (np.mean(dvdt) - 2*stdev), color = 'r', linestyle = '--')
    threshold = - 2*stdev
    global time_onset
    time_onset = []
    amplitude  = []
    qc         = []

    rise_time_100 = []
    rise_time_90  = []
    root = Tk()
    for j in range(0,len(results)):
        onset_change = 'f'
        while True:
            left_border = results.iloc[j].upstroke_index-10
            right_border =  results.iloc[j].upstroke_index+800
            right_border2 = right_border + 2000
            aoi = dvdt[left_border:right_border]
            if min(aoi) <= threshold:
                epsp_onset = np.where(aoi < threshold)[0][0] 
                amplitude_time = np.where(post_syn_y[left_border:right_border] == min(post_syn_y[left_border:right_border]))[0][0]
                ax[1].plot(post_syn_x[left_border:right_border][epsp_onset], post_syn_y[left_border:right_border][epsp_onset], 'x')
                ax[1].plot(post_syn_x[left_border:right_border][amplitude_time],min(post_syn_y[left_border:right_border]), 'x')
                #now get the times
                #get rise times
                amp = post_syn_y[left_border:right_border2][amplitude_time] - post_syn_y[left_border:right_border][epsp_onset]
                amp_10 = post_syn_y[left_border:right_border][epsp_onset] + 0.1*amp
                amp_90 = post_syn_y[left_border:right_border2][epsp_onset] + 0.9*amp 
                #cooridanted of amp_10 - ampg_90
                amp_10_x = np.where(post_syn_y[left_border:right_border] < amp_10)[0][0]
                amp_90_x = np.where(post_syn_y[left_border:right_border2] < amp_90)[0][0]
                ax[1].plot(post_syn_x[left_border:right_border2][epsp_onset:amplitude_time], post_syn_y[left_border:right_border2][epsp_onset:amplitude_time], color = 'r')
                #0-100 rise time
                rise0_100  =  post_syn_x[left_border:right_border2][amplitude_time]- post_syn_x[left_border:right_border2][epsp_onset]
                #plot and get 10 - 90 % rise time 
                ax[1].plot(post_syn_x[left_border:right_border2][amp_10_x:amp_90_x], post_syn_y[left_border:right_border2][amp_10_x:amp_90_x], color = 'cyan')
                rise10_90 = post_syn_x[left_border:right_border2][amp_90_x] - post_syn_y[left_border:right_border2][amp_10_x]
                #ask for user input regarding the events 
                fig.canvas.draw()    
                root.deiconify()
                res=mb.askquestion('sweep approval', 'Do you like the onset?')
                if res == 'no':
                    onset_change = 'f'
                    factor = simpledialog.askfloat('stdv', 'change std threshold level')
                    threshold =np.mean(dvdt)- factor *stdev
                    ax[2].axhline(y = (np.mean(dvdt) - factor*stdev), color = 'g', linestyle = '--')
                elif res == 'yes': 
                    onset_change = 'yes'
                    break
                pass_fail = simpledialog.askstring('threshold detection', 'does the onset pass ? (p/f)')
                root.withdraw()
                if pass_fail == 'p':
                    time_onset.append((post_syn_x[left_border:right_border][epsp_onset] - results.iloc[j].upstroke_t)*1e3)
                    amplitude.append(max(post_syn_y[left_border:right_border])-post_syn_y[left_border:right_border][epsp_onset])
                    qc.append('pass')
                    rise_time_100.append(rise0_100)
                    rise_time_90.append(rise10_90)
                    term = 'done'

                elif pass_fail == 'f':
                    time_onset.append((post_syn_x[left_border:right_border][epsp_onset] - results.iloc[j].upstroke_t)*1e3)
                    amplitude.append(max(post_syn_y[left_border:right_border])-post_syn_y[left_border:right_border][epsp_onset])
                    qc.append('fail')
                    rise_time_100.append(rise0_100)
                    rise_time_90.append(rise10_90)
                    term = 'done'

            else:
                messagebox.showinfo(title='no event', message = 'no event exceeding threshold found')
                time_onset.append(None)
                amplitude.append(None)
                qc.append(None)
                rise_time_100.append(None)
                rise_time_90.append(None)
                break
    plt.close()
    Time_Onset = {'latency':time_onset, 'amplitude':amplitude, 'qc':qc, 'rise_time_0_100': rise_time_100, 'rise_time_10_90': rise_time_90}   
    time_onset = pd.DataFrame(Time_Onset)
    return time_onset
  

def get_epsp_parameters_ttl_exc(post_syn_x, post_syn_y, peaks_ttl, average_dvdt, ax, fig):
    stdev = np.std(average_dvdt[50000:55000])
    ax[2].axhline(y = (np.mean(average_dvdt) + 4*stdev), color = 'r', linestyle = '--')
    threshold = 4*stdev
    time_onset = []
    amplitude  = []
    qc         = []
    rise_time_100 = []
    rise_time_90  = []
    root = Tk()
    decay = []
    for j in range(0,len(peaks_ttl)):
        left_border = peaks_ttl[j]+1000
        #onset
        right_border =  peaks_ttl[j]+10000
        right_border2 = right_border + 20000
        right_border3 = int(right_border + 39000)
        aoi = average_dvdt[left_border:right_border]
        #find the point where dvdt exceeds 30 mv/ms
        if max(aoi) >= threshold:
            epsp_onset = np.where(aoi > threshold)[0][0]
            amplitude_time = np.where(post_syn_y[left_border:right_border] == max(post_syn_y[left_border:right_border]))[0][0]
            ax[1].plot(post_syn_x[left_border:right_border][epsp_onset], post_syn_y[left_border:right_border][epsp_onset], 'x')
            ax[1].plot(post_syn_x[left_border:right_border][amplitude_time],max(post_syn_y[left_border:right_border]), 'x')
            #get rise times
            amp = post_syn_y[left_border:right_border2][amplitude_time] - post_syn_y[left_border:right_border][epsp_onset]
            amp_10 = post_syn_y[left_border:right_border][epsp_onset] + 0.1*amp
            amp_90 = post_syn_y[left_border:right_border2][epsp_onset] + 0.9*amp 
            #cooridanted of amp_10 - ampg_90
            amp_10_x = np.where(post_syn_y[left_border:right_border] > amp_10)[0][0]
            amp_90_x = np.where(post_syn_y[left_border:right_border2] > amp_90)[0][0]
            ax[1].plot(post_syn_x[left_border:right_border2][epsp_onset:amplitude_time], post_syn_y[left_border:right_border2][epsp_onset:amplitude_time], color = 'r')
            #0-100 rise time
            rise0_100  =  post_syn_x[left_border:right_border2][amplitude_time]- post_syn_x[left_border:right_border2][epsp_onset]
            #plot and get 10 - 90 % rise time 
            ax[1].plot(post_syn_x[left_border:right_border2][amp_10_x:amp_90_x], post_syn_y[left_border:right_border2][amp_10_x:amp_90_x], color = 'cyan')
            rise10_90 = post_syn_x[left_border:right_border2][amp_90_x] - post_syn_x[left_border:right_border2][amp_10_x]
            #plot and get decat decay area
            ax[1].plot(post_syn_x[left_border:right_border3][amplitude_time:right_border3], post_syn_y[left_border:right_border3][amplitude_time:right_border3], color='magenta')
            tau = fit_decay(amplitude_time, left_border, right_border3, post_syn_y, post_syn_x)
            fig.canvas.draw()    
            root.deiconify()
            pass_fail = simpledialog.askstring('threshold detection', 'does the onset pass ? (p/f)', parent=root)
            root.withdraw()
            if pass_fail == 'p':
                time_onset.append((post_syn_x[left_border:right_border][epsp_onset] - peaks_ttl[j])*1e3)
                amplitude.append(max(post_syn_y[left_border:right_border])-post_syn_y[left_border:right_border][epsp_onset])
                qc.append('pass')
                rise_time_100.append(rise0_100)
                rise_time_90.append(rise10_90)
                decay.append(tau)

            elif pass_fail == 'f':
               time_onset.append((post_syn_x[left_border:right_border][epsp_onset] - peaks_ttl[j])*1e3)
               amplitude.append(max(post_syn_y[left_border:right_border])-post_syn_y[left_border:right_border][epsp_onset])
               qc.append('fail')
               rise_time_100.append(rise0_100)
               rise_time_90.append(rise10_90)
               decay.append(tau)

        elif  max(aoi) <= threshold:
            messagebox.showinfo(title='no event', message = 'no event exceeding threshold found')
            time_onset.append(None)
            amplitude.append(None)
            qc.append(None)
            rise_time_100.append(None)
            rise_time_90.append(None)
            

    plt.close()
    Time_Onset = {'latency':time_onset, 'amplitude':amplitude, 'qc':qc, 'rise_time_0_100': rise_time_100, 'rise_time_10_90': rise_time_90, 'tau': decay}   
    time_onset = pd.DataFrame(Time_Onset)
    return time_onset  

