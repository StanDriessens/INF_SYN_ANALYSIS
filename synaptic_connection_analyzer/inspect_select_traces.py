# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 11:34:09 2024

select single traces 



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



def inspect_select_traces(f, sweep_tab_pre, sweep_tab_post):
    sweeps = [] 
    #plot all the sweeps indiviually 
    for i in sweep_tab_pre.sweep_number:
        fig,ax = plt.subplots(2,1)
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
        #plot the sweeps
        ax[0].plot(t, v_pre)
        ax[1].plot(t, v_post)
        #prompt for sweep acceptance 
        root = tkinter.Tk()
        include = tkinter.messagebox.askquestion(
       "Accept Sweep?",
       "Do you want to include sweep in your analysis?",
           icon="warning")
        
        if include == 'yes':
            sweeps.append(i)
        root.destroy()
        plt.close(fig)

    return(sweeps)
     
        
    
    
     
