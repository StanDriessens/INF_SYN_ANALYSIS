# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 13:33:39 2024

fit decay time

@author: sdr267
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def fit_decay(left_border, right_border3, post_syn_y, post_syn_x):
    def exponential_decay(t, A, tau, C):
        return A*np.exp(-t/tau) + C
        
    #FIND THE MAXIMUM IN THE CURRENT PLOT 
    max_index = np.where(post_syn_y[int(left_border):int(right_border3)] == max(post_syn_y[int(left_border):int(right_border3)]))[0][0] + left_border
    t_data = post_syn_x[max_index:int(right_border3)]
    response_data = post_syn_y[max_index:int(right_border3)]
    initial_guess = (0.4, 0.150, -69.1)
    params, covariance = curve_fit(exponential_decay, t_data, response_data, p0=initial_guess, maxfev=4000)
    decay_time_constant = params[1]
    plt.figure()
    plt.plot(t_data, response_data, label='Original Data', color = 'y')
    plt.plot(t_data, exponential_decay(t_data, *params), label='Fitted Curve', color='red')
    return decay_time_constant

    
