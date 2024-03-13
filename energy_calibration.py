'''this code calibrates spectrum binning from channel to energy, can be used with .Chn and .Spe file types
    requires use of read_chn.py (written by Oscar Tegmyr) for reading in .Chn files'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.stats import norm
from scipy import optimize
import read_chn as rc

# read_data function to get x/y data
# fit gaussian to get mean
# fit energy eq. and plot results
def read_data(filepath,head,foot):
    '''read in .Chn and .Spe data, gives (x,y)=(channel,data output) arrays for fitting'''
    
    if filepath[-4:]=='.Spe':
        data = pd.read_csv(filepath, dtype=np.int32, engine='python', skiprows=head, skipfooter=foot)
        data_ar = np.asarray(data)
        y = data_ar.flatten()
        x = np.arange(0,len(y),1)
        
    if filepath[-4:]=='.Chn':
        data = rc.gamma_data(filepath)
        y = data.hist_array
        x = np.arange(0,len(y),1)
        
    return x,y

# functions to be used in the curve fit and energy fit
def gauss(x, mean, amp, stdev, A, center):
    '''used to fit spectrum peak - added step function'''
    '''A=step amplitude, center=step mean'''
    fit = []
    for i in range(len(x)):
        if x[i]<=center:
            step = A
        if x[i]>center:
            step = 0
        g = amp * np.exp(-((x[i]-mean)**2) / (2 * stdev**2)) + step
        fit.append(g)
    return fit

def efit(x,m,b): # could be improved by including higher-order polynomial but need more data points
    '''energy = m*channel + b'''
    return m*x + b
    
# optimization functions
def gauss_opt(data,upper_bound,lower_bound,p0_guess):
    '''select data, apply curve_fit and print results'''
    # select data to curve-fit
    data_y = data[lower_bound:upper_bound]        
    # creating array of x-values for fit
    data_x = np.arange(lower_bound,upper_bound,1)
    
    # curve-fit to find peak
    var, cov = optimize.curve_fit(gauss, data_x, data_y, p0=p0_guess) #p0_guess = (mean, amp, stedev)
    # calculate full-width half-mass
    fwhm = 2.35*var[2]
    # calculate fitted curve
    data_fit = gauss(data_x, *var)
    # print results
    print('curve-fit results')
    print(var, '[mean, amp, stdev, step amp, step mean]')
    print('errors')
    print(cov)
    print('FWHM [channel]')
    print(fwhm)
    
    return data_x,data_y,data_fit,var[0],fwhm #x,y,fit,mean,fwhm

def energy_opt(peak1, peak2, mean1, mean2, x1, x2):
    '''apply curve_fit to energy function and return energy calibrated data'''
    xvals = np.array([mean1, mean2])
    yvals = np.array([peak1, peak2])
    values, cov = optimize.curve_fit(efit, xvals, yvals, p0=[1,5])
    
    print('curve-fit results:')
    print(values,'[m,b]')
    
    '''convert x-arrays from channel to energy'''
    x1_energy = efit(x1,*values)
    x2_energy = efit(x2,*values)
    
    return x1_energy,x2_energy,values
    
# plotting function
def plot_gauss(x,y,fit,mean,fwhm,source_name):
    
    fig, ax1 = plt.subplots(1,1, figsize=(8,4))

    ax1.plot(x, y, c='black', linewidth=1, label='data')
    ax1.plot(x, fit, c='red', lw=1, label='fit')
    ax1.axvline(mean, c='blue', linestyle='dashed', lw=0.8, label='mean')
    ax1.set_ylabel('count')
    ax1.set_xlabel('channel')
    ax1.minorticks_on()
    ax1.legend()

    fig.tight_layout()

    plt.show()
    
def plot_efit(x1,y1,x2,y2,peak1,peak2,source1,source2):
    
    fig, (ax1, ax2) = plt.subplots(2,1, figsize=(8,8))

    ax1.plot(x1, y1, c='black', linewidth=1, label='data')
    ax1.axvline(peak1, c='red', linestyle='dashed', lw=0.8, label=source1)
    ax1.axvline(peak2, c='orange', linestyle='dashed', lw=0.8, label=source2)
    ax1.set_ylabel('count')
    ax1.set_xlim(0,max(x1))
    ax1.minorticks_on()
    ax1.legend()

    ax2.plot(x2, y2, c='black', linewidth=1, label='data')
    ax2.axvline(peak2, c='orange', linestyle='dashed', lw=0.8, label=source2)
    ax2.axvline(peak1, c='red', linestyle='dashed', lw=0.8, label=source1)
    ax2.set_xlabel('energy [keV]')
    ax2.set_ylabel('counts')
    ax2.set_xlim(0,max(x2))
    ax2.minorticks_on()
    ax2.legend()

    fig.tight_layout()

    plt.show()
    
# fit functions 
def spect(source_name, data, peak, upper_bound, lower_bound, p0_guess, plotting=False):
    '''read in data and plot before inputting into function'''
    '''data inputting here refers to y-data from read_data function'''
    data_x,data_y,data_fit,mean,fwhm = gauss_opt(data,upper_bound,lower_bound,p0_guess)
    
    # plot fitted spectrum
    if plotting==True:
        plot_gauss(data_x,data_y,data_fit,mean,fwhm,source_name)
    
    return mean,fwhm
    
def energy_fit(peak1, peak2, mean1, mean2, x1, x2, y1, y2, fwhm1, fwhm2, source1, source2, plotting=False):
    '''peak1-2 known, mean1-2 are corresponding peak channels'''
    
    x1_eval,x2_eval,values = energy_opt(peak1,peak2,mean1,mean2,x1,x2)
    fwhm1_e = efit(fwhm1,*values)
    fwhm2_e = efit(fwhm2,*values)
    print('FWHM [keV]')
    print(source1,fwhm1_e)
    print(source2,fwhm2_e)
    
    '''use source1-2 to label plots of energy calibrated spectrum'''
    if plotting==True:
        plot_efit(x1_eval,y1,x2_eval,y2,peak1,peak2,source1,source2)
    
    return values