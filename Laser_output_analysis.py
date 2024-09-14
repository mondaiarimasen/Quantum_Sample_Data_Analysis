'''
Calculates parameters from laser locking data, e.g. slope of linear fit 
Parameters are used for calibration of laser

Used hypothesis test to see if slope is significant or not
'''

import matplotlib.pyplot as plt
import matplotlib.ticker
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

date = '20210226' # date of data collection
date_today = '20210516' # for saving

colors = ['k+','b.','gx','c*','ms','yD', 'r>', 'ko','b^','gv','cd', 
'mP','yo', 'r+', 'kp', 'bs', 'gX', 'c<', 'mo', 'y*', 'r^']

line_style = ['b-','r-','g-','b--','r--', 'g--', 'b-.', 'r-.','g-.','c:']

font_size=18

gen_data_folder = '/Data' + '/' + date
data_suffix = '0073' 
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


avg_wavelen_lbnd = 195120.310 # in GHz
avg_wavelen_ubnd = 195121.629 # in GHz

print()
print("Looking at ", gen_data_folder + spec_data_folder)


########################################################################

# Below is code to see laser drift / scanning range


fileName = '/20mHz5min1509nm.csv'

pd_df_laser_dri_data = pd.read_csv(os.getcwd() + gen_data_folder + fileName)

nparr_wavelen = pd_df_laser_dri_data.iloc[:,1].values
nparr_timestmp = pd_df_laser_dri_data.iloc[:,0].values

num_nparr_wavelen = len(nparr_wavelen)
print("This is num_nparr_wavelen: ", num_nparr_wavelen)

wavelen_lbnd_est =  198618.5 #1535.44 #
wavelen_ubnd_est =  198648.08 #1535.65 #
arr_lbnd_wavelens = np.empty(0)
arr_ubnd_wavelens = np.empty(0)

for i in range(0, num_nparr_wavelen):
    test = nparr_wavelen[i]
    if test <= wavelen_lbnd_est:
        arr_lbnd_wavelens = np.append(arr_lbnd_wavelens, test)
    if test >= wavelen_ubnd_est:
        arr_ubnd_wavelens = np.append(arr_ubnd_wavelens, test)

print("This is arr_lbnd_wavelens:", arr_lbnd_wavelens)
print("This is arr_ubnd_wavelens:", arr_ubnd_wavelens)

avg_wavelen_lbnd = np.mean(arr_lbnd_wavelens)
avg_wavelen_ubnd = np.mean(arr_ubnd_wavelens)
swp_rng = avg_wavelen_ubnd - avg_wavelen_lbnd

print("This is avg_wavelen_lbnd (GHz): ", avg_wavelen_lbnd)
print("This is avg_wavelen_ubnd (GHz): ", avg_wavelen_ubnd)

#print(pd_df_laser_dri_data)

# Plotting
fig_laser_drift = plt.figure(figsize=(12,8))
ax_laser_drift = plt.subplot()

plotTitle = "Laser Drift / Scanning range"
xAxisLabel = "Time (s)"
yAxisLabel = "Wavelength (GHz)"

ax_laser_drift.plot(nparr_timestmp, nparr_wavelen, line_style[0])
ax_laser_drift.plot(nparr_timestmp, [avg_wavelen_lbnd for i in range(0, 
    num_nparr_wavelen)], line_style[1], label = "Est. lower bnd.")
ax_laser_drift.plot(nparr_timestmp, [avg_wavelen_ubnd for i in range(0, 
    num_nparr_wavelen)], line_style[2], label = "Est. upper bnd.")

ax_laser_drift.set_title(plotTitle)
ax_laser_drift.set_ylabel(yAxisLabel)
ax_laser_drift.set_xlabel(xAxisLabel)

ax_laser_drift.legend()
plt.show()


########################################################################

time_lbnd_lin = -0.002
time_ubnd_lin = -0.00

# array of data folder names 
arr_fileName = np.array([])
arr_fileName = np.append(arr_fileName, os.getcwd() + gen_data_folder + spec_data_folder + file_name_prefix + '1' + '.CSV')
arr_fileName = np.append(arr_fileName, os.getcwd() + gen_data_folder + spec_data_folder + file_name_prefix + '3' + '.CSV')
arr_fileName = np.append(arr_fileName, os.getcwd() + gen_data_folder + spec_data_folder + file_name_prefix + '4' + '.CSV')

# contains the entire csv
pd_df_allCh1 = pd.read_csv(arr_fileName[0])
pd_df_allCh3 = pd.read_csv(arr_fileName[1])
pd_df_allCh4 = pd.read_csv(arr_fileName[1])

# stores the specified column (0 indexing) as a NP array
nparr_Time_CH1 = pd_df_allCh1.iloc[:,3].values
nparr_Time_CH3 = pd_df_allCh3.iloc[:,3].values
nparr_Time_CH4 = pd_df_allCh4.iloc[:,3].values

#if np.array_equal(nparr_Time_CH1, nparr_Time_CH3) is False:
#    raise Exception("Time axis for Channel 1 and 3 do not match")
if np.array_equal(nparr_Time_CH1, nparr_Time_CH4) is False:
    raise Exception("Time axis for Channel 1 and 4 do not match")


time_len = nparr_Time_CH4.size
print("first time element: ", nparr_Time_CH4[0])
time_int = nparr_Time_CH4[1] - nparr_Time_CH4[0]
print("time_int: ", time_int)
for i in range(1,time_len):
    if np.isclose(nparr_Time_CH4[i-1] + time_int, nparr_Time_CH4[i], rtol = 1e-08, atol = 1e-08) is False:
        print("i = ", i, "; nparr_Time_CH4[i-1]: ", nparr_Time_CH4[i-1], "; time_int: ", time_int, "; nparr_Time_CH4[i]: ", nparr_Time_CH4[i])
        raise Exception("Time axis is not uniformly spaced\n")


lbnd_lin_regime = np.where(nparr_Time_CH4 > time_lbnd_lin)[0][0]
ubnd_lin_regime = np.where(nparr_Time_CH4 > time_ubnd_lin)[0][0]


nparr_Voltage_CH1 = pd_df_allCh1.iloc[:,4].values#/y_scale_CH1 
nparr_Voltage_CH3 = pd_df_allCh3.iloc[:,4].values#/y_scale_CH3
nparr_Voltage_CH4 = pd_df_allCh4.iloc[:,4].values#/y_scale_CH4
print("nparr_Voltage_CH1.size = ", nparr_Voltage_CH1.size)
print("nparr_Voltage_CH4.size = ", nparr_Voltage_CH4.size)

########################################################################

## Reading from second folder to overlay two graphs

print("")

time_lbnd_lin_2 = -0.004
time_ubnd_lin_2 = -0.001

# array of data folder names 
arr_fileName_2 = np.array([])
arr_fileName_2 = np.append(arr_fileName_2, os.getcwd() + '/Data/20210226/All0071/F0071CH1.CSV')
#arr_fileName = np.append(arr_fileName, os.getcwd() + gen_data_folder + spec_data_folder + file_name_prefix + '3' + '.CSV')
arr_fileName_2 = np.append(arr_fileName_2, os.getcwd() + '/Data/20210226/All0071/F0071CH4.CSV')

# contains the entire csv
pd_df_allCh1_2 = pd.read_csv(arr_fileName_2[0])
#pd_df_allCh3 = pd.read_csv(arr_fileName[1])
pd_df_allCh4_2 = pd.read_csv(arr_fileName_2[1])

# stores the specified column (0 indexing) as a NP array
nparr_Time_CH1_2 = pd_df_allCh1_2.iloc[:,3].values
#nparr_Time_CH3 = pd_df_allCh3.iloc[:,3].values
nparr_Time_CH4_2 = pd_df_allCh4_2.iloc[:,3].values

#if np.array_equal(nparr_Time_CH1, nparr_Time_CH3) is False:
#    raise Exception("Time axis for Channel 1 and 3 do not match")
if np.array_equal(nparr_Time_CH1_2, nparr_Time_CH4_2) is False:
    raise Exception("Second data set: Time axis for Channel 1 and 4 do not match")


time_len_2 = nparr_Time_CH4_2.size
print("second dataset: first time element: ", nparr_Time_CH4_2[0])
time_int_2 = nparr_Time_CH4_2[1] - nparr_Time_CH4_2[0]
print("time_int_2: ", time_int_2)
for i in range(1,time_len_2):
    if np.isclose(nparr_Time_CH4_2[i-1] + time_int_2, nparr_Time_CH4_2[i], rtol = 1e-08, atol = 1e-08) is False:
        print("i = ", i, "; nparr_Time_CH4_2[i-1]: ", nparr_Time_CH4_2[i-1], "; time_int_2: ", time_int_2, "; nparr_Time_CH4_2[i]: ", nparr_Time_CH4_2[i])
        raise Exception("Time axis is not uniformly spaced\n")


lbnd_lin_regime_2 = np.where(nparr_Time_CH4_2 > time_lbnd_lin_2)[0][0]
ubnd_lin_regime_2 = np.where(nparr_Time_CH4_2 > time_ubnd_lin_2)[0][0]


nparr_Voltage_CH1_2 = pd_df_allCh1_2.iloc[:,4].values#/y_scale_CH1 
#nparr_Voltage_CH3 = pd_df_allCh3.iloc[:,4].values#/y_scale_CH3
nparr_Voltage_CH4_2 = pd_df_allCh4_2.iloc[:,4].values#/y_scale_CH4
#print("nparr_Voltage_CH1.size = ", nparr_Voltage_CH1.size)
print("nparr_Voltage_CH4_2.size = ", nparr_Voltage_CH4_2.size)





########################################################################

# defines functions
def linear(x, a,b):
  return a*x + b

def lorentizan_sin(x, a, w, l, b, f, p, c, d): # a is amplitude,  w is full width at half max, l is resonance absorption wavelength/freq of dip
    return a * np.power(w/2.0,2.0) / (np.power(w/2.0,2.0) + np.power(x-l,2.0) ) + b * np.sin(2.0*np.pi * f * x + p) + c * x + d

# for aligning y-axes
def align_yaxis(ax1, v1, ax2, v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny+dy, maxy+dy)



########################################################################

# dispersion slope time
slp_time = nparr_Time_CH4[ubnd_lin_regime] - nparr_Time_CH4[lbnd_lin_regime]


# Plotting
fig = plt.figure(figsize=(12,8))
ax1 = plt.subplot()

plotTitle = "Channel Data for " + gen_data_folder + spec_data_folder + " and ALL0071"
xAxisLabel = "Time (s)"
ax1.set_xlabel(xAxisLabel, fontsize=font_size)
ax1.tick_params(axis='both', which='both', direction='in', 
    right=True, top=True, labelsize=font_size)

swp_rng = (avg_wavelen_ubnd - avg_wavelen_lbnd)
print("sweeping range (GHz): ", swp_rng)

X_lin_fit = np.linspace(nparr_Time_CH4[lbnd_lin_regime] * 1.0, nparr_Time_CH4[ubnd_lin_regime] * 1.0, 600)

### Finding linear part of the dispersion curve
popt, pcov = optimize.curve_fit(linear, nparr_Time_CH4[lbnd_lin_regime:ubnd_lin_regime], 
                                nparr_Voltage_CH1[lbnd_lin_regime:ubnd_lin_regime])
print("fit parameters linear:")
print(popt)
print(pcov)

# Extracting the slope and its standard error from the fit
slope = popt[0]
slope_std_error = np.sqrt(pcov[0, 0])

# Hypothesis test (t-test) for slope = 0
t_statistic = slope / slope_std_error
df = len(nparr_Time_CH4[lbnd_lin_regime:ubnd_lin_regime]) - 2  # degrees of freedom
p_value = 2 * (1 - stats.t.cdf(abs(t_statistic), df))

print(f"T-statistic for the slope: {t_statistic}")
print(f"P-value for the slope significance test: {p_value}")

if p_value < 0.05:
    print("The slope is significantly different from 0 (reject H0)")
else:
    print("The slope is not significantly different from 0 (fail to reject H0)")

# Plot the linear fit
ax1.plot(nparr_Time_CH4[lbnd_lin_regime:ubnd_lin_regime], nparr_Voltage_CH1[lbnd_lin_regime:ubnd_lin_regime], 'bo')
ax1.plot(X_lin_fit, linear(X_lin_fit, *popt), 'r-', label='Linear Fit')

# Repeat the process for the second dataset
slp_time_2 = nparr_Time_CH4_2[ubnd_lin_regime_2] - nparr_Time_CH4_2[lbnd_lin_regime_2]
X_lin_fit_2 = np.linspace(nparr_Time_CH4_2[lbnd_lin_regime_2]*1.1, nparr_Time_CH4_2[ubnd_lin_regime_2]*1.1, 600)

popt_2, pcov_2 = optimize.curve_fit(linear, nparr_Time_CH4_2[lbnd_lin_regime_2:ubnd_lin_regime_2], 
                                    nparr_Voltage_CH1_2[lbnd_lin_regime_2:ubnd_lin_regime_2])
print("Second data set fit parameters linear:")
print(popt_2)
print(pcov_2)

# Extract the slope and its standard error for the second fit
slope_2 = popt_2[0]
slope_std_error_2 = np.sqrt(pcov_2[0, 0])

# Hypothesis test for the second dataset slope = 0
t_statistic_2 = slope_2 / slope_std_error_2
df_2 = len(nparr_Time_CH4_2[lbnd_lin_regime_2:ubnd_lin_regime_2]) - 2
p_value_2 = 2 * (1 - stats.t.cdf(abs(t_statistic_2), df_2))

print(f"T-statistic for the second slope: {t_statistic_2}")
print(f"P-value for the second slope significance test: {p_value_2}")

if p_value_2 < 0.05:
    print("The second slope is significantly different from 0 (reject H0)")
else:
    print("The second slope is not significantly different from 0 (fail to reject H0)")

# Plot the linear fit for the second dataset
ax1.plot(nparr_Time_CH4_2[lbnd_lin_regime_2:ubnd_lin_regime_2], nparr_Voltage_CH1_2[lbnd_lin_regime_2:ubnd_lin_regime_2], 'go')
ax1.plot(X_lin_fit_2, linear(X_lin_fit_2, *popt_2), 'r-', label='Second Linear Fit')

plt.show()






########################################################################


## Fitting to Lorentizan

# Plotting
fig_lor = plt.figure(figsize=(12,8))
ax1_lor = plt.subplot()

plotTitle = "Channel Data"
xAxisLabel = "Time (s)"
ax1_lor.set_xlabel(xAxisLabel)
ax1_lor.tick_params(axis='both', which='both', direction='in')
ax1_lor.xaxis.set_minor_locator(AutoMinorLocator())


### Convert time of sweep to frequency range
time_ind_min = np.argmin(nparr_Voltage_CH4)
time_ind_max = np.argmax(nparr_Voltage_CH4)
print("nparr_Time_CH4[np.argmin(nparr_Time_CH4)]: ", nparr_Time_CH4[time_ind_min])
print("nparr_Time_CH4[np.argmax(nparr_Time_CH4)]: ", nparr_Time_CH4[time_ind_max])

time_min = nparr_Time_CH4[time_ind_min]
time_max = nparr_Time_CH4[time_ind_max]

sweep_lowbnd = 195186.216 # GHz, from wavemeter
sweep_upbnd = 195210.0 # GHz, from wavemeter

def time_to_freq(t,tl, tu): # converts between time of sweep and frequency
    return 1.0 * (sweep_upbnd - sweep_lowbnd) / ((tu-tl) * 1.0) * t + sweep_lowbnd - 1.0 * (sweep_upbnd - sweep_lowbnd) / ((tu-tl)*1.0) * tl

middle_freq = time_to_freq((time_min + time_max)/2.0, time_min, time_max)
print("middle freq: ", middle_freq)
if np.isclose(middle_freq, (sweep_lowbnd + sweep_upbnd)/2.0, rtol = 1e-08, atol = 1e-08) is False:
        raise Exception("middle_freq is not close\n")

a = 100
b = -100

freq_array = time_to_freq(nparr_Time_CH4[time_ind_min + a :time_ind_max+b], time_min, time_max)
print("freq_array[0]: ", freq_array[0])
print("freq_array[-1]: ", freq_array[-1])


print(type(freq_array))
print(type(nparr_Voltage_CH4))

X_lor_fit = np.linspace(freq_array[0], freq_array[-1], 2000)


### Fitting Lorentizan to Absorption curve
popt, pcov = optimize.curve_fit(lorentizan_sin, freq_array, nparr_Voltage_CH1[time_ind_min + a :time_ind_max+b], maxfev=1000, p0=[-0.18, 6, 195198, 0.07, 800, 0, 0.035, 0.2])
print("These are parameters of Lorentzian fit")
print(popt)
print(pcov)


#ax1.plot(nparr_Time_CH3, nparr_Voltage_CH3, line_style[-1], label="Sweep of Laser (Ch3)")
yAxis1Label = "Ch1 Voltage (V)"
ax1_lor.set_ylabel(yAxis1Label)
ax1_lor.yaxis.set_minor_locator(AutoMinorLocator())

ax2_lor = ax1_lor.twinx() # instantiate a second axes that shares the same x-axis
ax2_lor.tick_params(axis='both', which='both', direction='in')



ax1_lor.plot(freq_array, nparr_Voltage_CH1[time_ind_min + a :time_ind_max+b], line_style[0], label="Ch1")
ax2_lor.plot(freq_array, nparr_Voltage_CH4[time_ind_min + a :time_ind_max+b], line_style[1], label="Error signal (Ch4, uncropped)")

ax1_lor.plot(X_lor_fit, lorentizan_sin(X_lor_fit, *popt), line_style[2], label='amp = %5.3f, FWHM = %5.3f GHz, l = %5.3f GHz, sin amp = %5.3f V, sin frq = %5.3f GHz, phase = %5.3f rad, offset slope = %5.3f V, offset inter = %5.3f' % tuple(popt))
yAxis2Label = "Ch4 Voltage (V)"
ax2_lor.set_ylabel(yAxis2Label)
ax2_lor.yaxis.set_minor_locator(AutoMinorLocator())

ax1_lor.set_title(plotTitle)

lines, labels = ax1_lor.get_legend_handles_labels()
lines2, labels2 = ax2_lor.get_legend_handles_labels()
ax2_lor.legend(lines + lines2, labels + labels2, loc=0)

plt.show()
fig.savefig('figs/err_sig_lock_laser2gas_cell_' 
 + date_today + '.pdf', bbox_inches = 'tight')









########################################################################
