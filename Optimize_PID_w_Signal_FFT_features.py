'''
Calculating the optimal PID settings, via spectral analysis of signal 
to find the resonant frequencies and computing their power spectrum
'''

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from scipy import optimize, signal, stats
from statistics import mean
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


import math
import os
import sys


sys.path.append(os.path.expanduser('/Data')) 

date = '20210208' # date of data collection
date_today = '20210301' # for saving

colors = ['k+','b.','gx','c*','ms','yD', 'r>', 'ko','b^','gv','cd', 
'mP','yo', 'r+', 'kp', 'bs', 'gX', 'c<', 'mo', 'y*', 'r^']

line_style = ['b-','r-','g-','b--','r--', 'g--', 'b-.', 'r-.','g-.','c:']

font_size = 18

gen_data_folder = '/Data' + '/' + date
data_suffix = '0031' 
spec_data_folder = '/All' + data_suffix
file_name_prefix = '/F' + data_suffix + 'CH'

## Experimental settings
laser_swp_frq = 10 # in Hz
laser_swp_Vpp = 4 # in Vpp; used in plot label
laser_wavelen = 1535.33 # in nm; used in plot label
laser_freq = 195310.0 # in GHz; used in plot label
eom_freq = 40.0 # in MHz; used in plot label
eom_Vpp = 4.5 # in Vpp; used in plot label
laser_power = 8 # in mW; used in plot label

err_sig_slope = 7.31e-03 # V/GHz, slope of error signal for 2/8/2021 1509nm
err_sig_slope_2 = 1.72e-4 # V/GHz, slope of error signal for 12/16/2020 1535nm


print("Looking at ", gen_data_folder + spec_data_folder) 

########################################################################

## Reading data and putting it into arrays
fileShortName = '1s'

fileName = gen_data_folder + spec_data_folder + '/' + fileShortName + '.txt'

fullFileName = os.getcwd() + fileName
print("Looking at ", fileName) 

# contains the entire csv
pd_df = np.loadtxt(fullFileName)

nparr_time_vals = pd_df[:,0]
nparr_voltage_vals = pd_df[:,1]

time_int = nparr_time_vals[1] - nparr_time_vals[0]
print("time_int: ", time_int)

maxval = max(nparr_voltage_vals)
mean_1 = mean(nparr_voltage_vals) 
nparr_voltage_vals = (nparr_voltage_vals - mean_1)/maxval

total_time = nparr_time_vals[-1] - nparr_time_vals[0]
print("total_time: ", total_time)

num_vals = len(nparr_time_vals)
print("numvals: ", num_vals)

# Plotting
fig = plt.figure(figsize=(12,8))
ax1 = plt.subplot()

plotTitle = "Channel Data"
xAxisLabel = "Time (s)"
yAxis1Label = "Voltage (V)"

ax1.plot(nparr_time_vals, nparr_voltage_vals, line_style[0], label="Signal")

ax1.set_xlabel(xAxisLabel)
ax1.tick_params(axis='both', which='both', direction='in', 
    right=True,top=True)
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.set_ylabel(yAxis1Label)
ax1.yaxis.set_minor_locator(AutoMinorLocator())

ax1.set_title(plotTitle)

lines, labels = ax1.get_legend_handles_labels()
ax1.legend(lines, labels, loc=0)

plt.show()
fig.savefig('figs/err_sig_lock_laser2gas_cell_' 
 + date_today + '.pdf', bbox_inches = 'tight')

########################################################################

## Saving data from second experiment

## Reading data and putting it into arrays
fileShortName_2 = '100s'

fileName = gen_data_folder + spec_data_folder + '/' + fileShortName_2 + '.txt'

fullFileName_2 = os.getcwd() + fileName_2
print("Looking at ", fileName_2) 

# contains the entire csv
pd_df = pd.read_csv(fullFileName, sep='\t', lineterminator='\r', skiprows=None, header = None)
pd_df_2 = np.loadtxt(fullFileName_2)

nparr_time_vals_2 = pd_df_2[:,0]
nparr_voltage_vals_2 = pd_df_2[:,1]

time_int_2 = nparr_time_vals_2[1] - nparr_time_vals_2[0]
print("time_in_2: ", time_int_2)
print("time_int - time_int_2: ", time_int - time_int_2)


for i in range(1,num_vals-1):
    if np.isclose(nparr_time_vals[i-1] + time_int, nparr_time_vals[i], rtol = 1e-08, atol = 1e-08) is False:
        print("i = ", i, "; nparr_time_vals[i-1]: ", nparr_time_vals[i-1], "; time_int: ", time_int, "; nparr_time_vals[i]: ", nparr_time_vals[i])
        raise Exception("Time axis is not uniformly spaced\n")


total_time_2 = nparr_time_vals_2[-1] - nparr_time_vals_2[0]
print("total_time_2: ", total_time_2)

num_vals_2 = len(nparr_time_vals_2)
print("numvals_2: ", num_vals_2)

num_vals = len(nparr_time_vals_2)

# Plotting
fig_2 = plt.figure(figsize=(12,8))
ax1_2 = plt.subplot()

plotTitle_2 = "Channel Data"
xAxisLabel_2 = "Time (s)"
yAxis1Label_2 = "Voltage (V)"

ax1_2.plot(nparr_time_vals_2, nparr_voltage_vals_2, line_style[0], label="Signal")

ax1_2.set_xlabel(xAxisLabel_2)
ax1_2.tick_params(axis='both', which='both', direction='in', 
    right=True,top=True)
ax1_2.xaxis.set_minor_locator(AutoMinorLocator())
ax1_2.set_ylabel(yAxis1Label_2)
ax1_2.yaxis.set_minor_locator(AutoMinorLocator())

ax1_2.set_title(plotTitle_2)

lines_2, labels_2 = ax1_2.get_legend_handles_labels()
ax1_2.legend(lines_2, labels_2, loc=0)

plt.show()
fig.savefig('figs/err_sig_lock_laser2gas_cell_' 
 + date_today + '.pdf', bbox_inches = 'tight')


########################################################################

# Plotting
fig_3 = plt.figure(figsize=(12,8))
ax1_3 = plt.subplot()

plotTitle_3 = "Unlocked and Locked noise time trace"
xAxisLabel_3 = "Time (s)"
yAxis1Label_3 = "Voltage (V)"

#ax1_3.plot(nparr_time_vals, nparr_voltage_vals, line_style[0], label="Unlocked")
ax1_3.plot(nparr_time_vals_2, nparr_voltage_vals_2, line_style[1], label="Locked")

ax1_3.set_xlabel(xAxisLabel_3)
ax1_3.tick_params(axis='both', which='both', direction='in', 
    right=True,top=True)
ax1_3.xaxis.set_minor_locator(AutoMinorLocator())
ax1_3.set_ylabel(yAxis1Label_3)
ax1_3.yaxis.set_minor_locator(AutoMinorLocator())

ax1_3.set_title(plotTitle_3)

lines_3, labels_3 = ax1_3.get_legend_handles_labels()
ax1_3.legend(lines_3, labels_3, loc=0)
fig_3.savefig('figs/overlayed_unlocked_noise_time_trace_' + date + '_col_' + date_today + '_saved_' + '.pdf', bbox_inches = 'tight')



########################################################################

def transfer_fn(x, kp, ki):
    return 20*np.log10(np.sqrt(np.power(kp,2.0) + np.power(ki/x,2.0)))


########################################################################

## Taking Fourier Transform

## testing Fourier Transform
start = 0
stop = 100
num_samples = 5000
x_even = np.linspace(start,stop, num_samples+1)
x = x_even[:-1]
dt = x[1]-x[0]
print(dt)
x_even_len = len(x_even)
freq1 = 2
freq2 = 14
freq3 = 17

y = 4*np.sin(2 * np.pi * freq1 * x) + 2 * np.sin(2 * np.pi * freq2 * x + np.pi/3) + 5 * np.sin(2*np.pi*freq3*x)

print(y.size)

fft_y = np.fft.rfft(y)*np.power(dt/num_samples,0.5) # this is sqrt of {power spectrum = fft^2 * (dt^2/ total time)} exact
fft_y_len = len(fft_y)
xf = np.fft.rfftfreq(num_samples, d=(stop-start)/num_samples)

print(x)
print(x_even_len)
print(fft_y_len)

sumP1 = np.sum(np.power(y,2))
sumP2 = (np.abs(fft_y[0])**2 + 2.0*np.sum(np.power(np.abs(fft_y[1:]),2)))/num_samples # this fft_y must be np.fft.rfft(y)

ifft_y = np.fft.irfft(fft_y,num_samples)

print(xf)


fig = plt.figure(figsize=(12,8))
ax = plt.subplot()

ydata = np.abs(fft_y[1:-1])
ax.plot(xf[1:-1], ydata)
plt.show()



## Fourier Tranform of data
time_lbnd_fft = 0.0254
time_ubnd_fft = 0.02595

lbnd_fft_regime = np.where(nparr_time_vals > time_lbnd_fft)[0][0]
print("lbnd_fft_regime: ", lbnd_fft_regime)
ubnd_fft_regime = np.where(nparr_time_vals > time_ubnd_fft)[0][0]
print("ubnd_fft_regime: ", ubnd_fft_regime)


lbnd_fft_regime = 0
ubnd_fft_regime = num_vals

fft_regime = nparr_voltage_vals[lbnd_fft_regime:ubnd_fft_regime] 
print("shape of fft_regime ", nparr_voltage_vals.shape)
fft_regime_len = fft_regime.size
print("This is fft_regime_len: ", fft_regime_len)

fft_data = np.fft.rfft(fft_regime)* np.power(time_int/fft_regime_len,0.5)*2.0 # normalized (see comment above in testing area)
fft_data_len = len(fft_data)
print("this is len of fft_data: ", fft_data_len)

xf = np.fft.rfftfreq(int(fft_regime_len), d=total_time/int(fft_regime_len)) # assuming taking fft of entire data
#xf = np.linspace(0,int(fft_data_len), int(fft_data_len))/total_time
print("This is xf: ", xf)


figfft = plt.figure(figsize=(12,8))
axfft = plt.subplot() 

ydata = np.abs(fft_data[1:-1])/err_sig_slope
print("len ydata: ", len(ydata))
print("this is ydata: ", ydata)


X_lin_fit = np.linspace(xf[lbnd_trans_regime], xf[ubnd_trans_regime], 600)
### Finding linear part of dispersion curve
popt, pcov = optimize.curve_fit(transfer_fn, xf[lbnd_trans_regime:ubnd_trans_regime], y_data[lbnd_trans_regime:ubnd_trans_regime])
print("fit parameters transferfn ")
print(popt)
print(pcov)


axfft.plot(xf[1:-1], ydata, 
    label="Signal \nLaser wvlngth: " + str(laser_wavelen) 
    + " nm \nLaser Freq: " + str(laser_freq) + " GHz\n")
axfft.set_xscale("log")
axfft.set_yscale("log")

axfft.tick_params(axis='both', which='both', direction='in', 
    right=True,top=True)
xAxisLabel = "Frequency (Hz)"
axfft.set_xlabel(xAxisLabel)
yAxis2Label = "Freq/sqrt(freq) (GHz/" + r'$\sqrt{Hz}$' + ")"
axfft.set_ylabel(yAxis2Label)
plotTitle = "FFT of " + fileName
axfft.set_title(plotTitle)

axfft.legend()

plt.show()
figfft.savefig('figs/fft_noisefreqspec_' + date + '_' + fileShortName + '_col_' + date_today + '_saved_' + '.pdf', bbox_inches = 'tight')


########################################################################
