# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 11:17:41 2024

make radiobutton channels

@author: sdr267
"""


"""

below create a radio buttont to select pre and post channels from found channels 
in abf

to do: similar for nwb files. 
"""

import pyabf
from tkinter import*

def make_radio_channels_abf(file):
    abf = pyabf.ABF(file)
    channel_names = abf.adcNames
    def get_pre_channel():
        global pre_channel
        pre_channel = channel_names[x.get()]
        window1.destroy()
    window1 = Tk()
    x = IntVar()
    for index in range(0,len(channel_names)):
        radiobutton = Radiobutton(window1, text=channel_names[index], variable=x, value=index, padx =25)
        radiobutton.pack(anchor=W)
    okay_button = Button(window1, text="Select pre Channel", command=get_pre_channel)
    okay_button.pack(pady=10)    
    window1.mainloop()    
    
    def get_post_channel():
        global post_channel
        post_channel = channel_names[x.get()]
        window2.destroy()
    window2 = Tk()
    x = IntVar()
    for index in range(0,len(channel_names)):
        radiobutton = Radiobutton(window2, text=channel_names[index], variable=x, value=index, padx =25)
        radiobutton.pack(anchor=W)
    okay_button = Button(window2, text="Select post Channel", command=get_post_channel)
    okay_button.pack(pady=10)    
    window1.mainloop()    
    
    return pre_channel, post_channel 
            
    
            
    
