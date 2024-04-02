# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 15:19:20 2024

@author: sdr267
"""

import matplotlib.pyplot as plt 
import numpy as np 
import scipy as sp
from scipy.signal import savgol_filter
from tkinter import messagebox as mb
import pyabf

def sweep_select(abf, pre_channel, post_channel):
    sweep_set = []
    for i in range(0,abf.sweepCount):
        fig,ax = plt.subplots(2,1, sharex=True)
        abf.setSweep(i,pre_channel)
        #fitler  the pre_syanptic signal 
        pre_syn_y = savgol_filter(abf.sweepY, 99 ,5)
        pre_syn_x = abf.sweepX
        pre_syn_i = abf.sweepC
        abf.setSweep(i,post_channel)
        #filter the post synaptic signal
        post_syn_y = savgol_filter(abf.sweepY, 251 ,2)
        post_syn_x = abf.sweepX
        ax[0].plot(pre_syn_x, pre_syn_y)
        ax[1].plot(post_syn_x, post_syn_y)
        res=mb.askquestion('sweep approval', 'Do you want to include your sweep?')
        if res == 'yes':
            sweep_set.append(i)
        elif res == 'no':
            print('sweep not included: ', i)
        plt.close()
    return sweep_set
        
        