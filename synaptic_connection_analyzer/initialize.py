# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 15:27:14 2024


Initialize synaptic data analyzer

prior to nwb ttl - synaptic or abf ttl synapic connection data. 


Stan Driessens 2024/01/08 

Note make sure you are in the correct directory when running the script. 

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
from nwb_ttl import run_nwb_ttl_analyzer
from nwb_ttl import multipatch_nwb_analyzer
from abf_spike import run_abf_spike_analyzer
from abf_spike import run_abf_spike_analyzer_average

 
plt.style.use('dark_background')

#tkinter.Tk().withdraw() # prevents an empty tkinter window from appearing
#open file
#open file
window = Tk()
file = filedialog.askopenfilename()
#window.destroy()
species_select = ['human', 'mouse']
def save_species():
    global species
    species = species_select[x.get()]
    window1.destroy()
window1 = Tk()
x = IntVar()
for index in range(len(species_select)):
    radiobutton = Radiobutton(window1, text=species_select[index], variable=x, value=index, padx=25)
    radiobutton.pack(anchor=W)
okay_button = Button(window1, text="Select Species", command=save_species)
okay_button.pack(pady=10)    
window1.mainloop()    

file_select = ['nwb', 'abf']
def file_selecting():
    global file_type
    file_type = file_select[x.get()]
    window2.destroy()
window2 = Tk()
x = IntVar()
for index in range(len(file_select)):
    radiobutton = Radiobutton(window2, text=file_select[index], variable=x, value=index, padx=25)
    radiobutton.pack(anchor=W)
okay_button = Button(window2, text="Select file type", command=file_selecting)
okay_button.pack(pady=10)    
window2.mainloop()

connection_types = ['exc', 'inh']
def select_connection():
    global connection
    connection = connection_types[x.get()]
    window3.destroy()
window3 = Tk()
x = IntVar()
for index in range(len(connection_types)):
    radiobutton = Radiobutton(window3, text=connection_types[index], variable=x, value=index, padx=25)
    radiobutton.pack(anchor=W)
okay_button = Button(window3, text="Select connection", command=select_connection)
okay_button.pack(pady=10)    
window3.mainloop()    


window4= Tk() 
average = mb.askyesno('average', 'Do you want to analyze an average trace?')
window4.destroy()

connection_types = ['multipatch', 'no multipatch']
def select_connection():
    global multipatch
    multipatch = connection_types[x.get()]
    window5.destroy()
window5 = Tk()
x = IntVar()
for index in range(len(connection_types)):
    radiobutton = Radiobutton(window5, text=connection_types[index], variable=x, value=index, padx=25)
    radiobutton.pack(anchor=W)
okay_button = Button(window5, text="Select multipatch", command=select_connection)
okay_button.pack(pady=10)    
window5.mainloop()  




if (average == True) and (file_type == 'abf'):
    run_abf_spike_analyzer_average(file, species, connection)
elif (file_type == 'nwb') and (multipatch == 'no multipatch'):
    run_nwb_ttl_analyzer(file, species)
    
elif (file_type == 'nwb') and (multipatch == 'multipatch'):
    multipatch_nwb_analyzer(file, connection)
    
elif file_type == 'abf':
    run_abf_spike_analyzer(file, species, connection)



