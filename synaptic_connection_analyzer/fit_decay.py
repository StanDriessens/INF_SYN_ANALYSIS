# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 13:33:39 2024

fit decay time

@author: sdr267
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def fit_decay(amplitude_time, left_border, right_border3, post_syn_y, post_syn_x):
    def exponential_decay(x, m, t, b):
        return m*np.exp(-t*x) + b
        
    #FIND THE MAXIMUM IN THE CURRENT PLOT 
    response_data = post_syn_y[left_border:right_border3][amplitude_time:right_border3]
    t_data = np.arange(0 , len(response_data),1)
    initial_guess = (2000, .1, 50)
    params, covariance = curve_fit(exponential_decay, t_data, response_data, p0=initial_guess, maxfev=4000)
    decay_time_constant = params[1]
    plt.figure()
    plt.plot(t_data, response_data, label='Original Data', color = 'y')
    plt.plot(t_data, exponential_decay(t_data, *params), label='Fitted Curve', color='red')
    return decay_time_constant

    
