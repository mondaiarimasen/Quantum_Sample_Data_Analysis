'''
Extract sweep data from the wavemeter, and transform from scope time 
to laser wavelength.
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import signal

### Modifiable Parameters ###
# These parameters prevent capturing unwanted local maxima.
# They are in nm, corresponding to the wavelength read by the wavemeter.
peak_threshold_upper = 1535.53  # Upper limit for peak detection (nm)
dip_threshold_lower = 1535.34  # Lower limit for dip detection (nm)
c = 299792458.
### Modifiable Parameters ###
input_file = '...'




df = pd.read_csv(input_file, skiprows=1)
df.columns = ['time_col', 'wavelength_col']
time_values = list(df.time_col)
wavelength_values = list(df.wavelength_col)
print(wavelength_values)

peak_positions = []

# Convert timestamp to seconds (HH:MM:SS to seconds)
def time_in_seconds(time_str):
    time_parts = time_str.split(":")
    hr = int(time_parts[0])
    min = int(time_parts[1])
    sec = float(time_parts[2])
    return sec + (min * 60) + (hr * 3600)

converted_times = []
initial_time = time_values[0].split("T", 1)[1]
zero_time_ref = time_in_seconds(initial_time)

# Convert the times relative to the initial timestamp
for t in time_values:
    curr_time = t.split("T", 1)[1]
    time_in_sec = time_in_seconds(curr_time)
    adjusted_time = time_in_sec - zero_time_ref
    converted_times.append(adjusted_time)

# Identify peaks (maxima) in the wavelength data
peaks_detected = signal.find_peaks(wavelength_values, height=peak_threshold_upper, distance=10)
peak_indices = peaks_detected[0]
real_peak_times = [converted_times[i] for i in peak_indices]

# Invert the wavelength values to find minima (dips)
wavelength_inverted = [-1.0 * w for w in wavelength_values]

# Identify dips (minima) in the inverted wavelength data
dips_detected = signal.find_peaks(wavelength_inverted, height=-1.0 * dip_threshold_lower, distance=10)
dip_indices = dips_detected[0]
real_dip_times = [converted_times[i] for i in dip_indices]

# Gather peak and dip heights
peak_wavelengths = peaks_detected[1].get('peak_heights')
dip_wavelengths = -1.0 * dips_detected[1].get('peak_heights')

# Calculate average and standard deviation for peaks and dips
avg_peak_wavelength = np.mean(peak_wavelengths)
std_peak_wavelength = np.std(peak_wavelengths)
avg_dip_wavelength = np.mean(dip_wavelengths)
std_dip_wavelength = np.std(dip_wavelengths)

num_peaks = len(peak_wavelengths)
num_dips = len(dip_wavelengths)

if __name__ == "__main__":
    plt.plot(converted_times, wavelength_values)
    plt.plot(real_peak_times, peak_wavelengths, 'go')
    plt.plot(real_dip_times, dip_wavelengths, 'yo')
    print('Number of Peaks: ' + str(num_peaks))
    print('Number of Dips: ' + str(num_dips))
    plt.show()
