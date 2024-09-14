'''
Process the relevant data files, apply transformations to center the absorption profile,
normalize values, and reduce noise artifacts.
'''

import numpy as np
import matplotlib.pyplot as plt

### Configuration Section ###
spectrum_filename = 'channel3.csv'  # Absorption spectrum file (csv or txt)
sweep_filename = 'channel2.csv'  # Sweep data file (csv or txt)
plot_title = '1535.2nm'  # Title for the final plot
start_range = 540  # Lower limit of the absorption range
end_range = 2070  # Upper limit of the absorption range
noise_correction = 0.0225  # Noise correction factor
### Configuration Section ###

# Initialize the resolution variable as None
vertical_scale = None

# Function to process the spectrum data and extract relevant parameters
def process_spectrum(file_name, start, end, correction):
    # Load the file and select relevant columns
    spectrum_data = np.genfromtxt(file_name, usecols=(0, 1, 3, 4), 
                                  dtype=('U20', 'U20', float, float), delimiter=',')

    # Extract time and voltage values based on the specified range
    time_segment = spectrum_data['f2'][start:end]
    voltage_segment = spectrum_data['f3'][start:end]
    
    # Apply noise correction to the voltage values
    corrected_voltages = [voltage - correction for voltage in voltage_segment]
    
    return time_segment, corrected_voltages, spectrum_data

# Function to find and set the vertical scale resolution
def find_resolution(data):
    for entry in data['f0']:
        if entry == 'Vertical Scale':
            return float(entry['f1'])
    return None

# Normalize voltage values by the maximum voltage
def normalize_voltages(voltages):
    max_value = max(voltages)
    if max_value != 0:
        return [v / max_value for v in voltages], max_value
    return voltages, max_value

# Direct logic execution starts here
time_data, voltage_data, full_data = process_spectrum(spectrum_filename, start_range, end_range, noise_correction)

# Check if the resolution was found, else print warning
vertical_scale = find_resolution(full_data)
if vertical_scale is None:
    print('Warning: Vertical scale not detected. Please set manually in vertical_scale.')

# Normalize the voltage data
voltage_data, max_voltage = normalize_voltages(voltage_data)

print(f'Max Voltage: {max_voltage}')
print(f'Normalized Voltages: {voltage_data}')

# If you'd like to plot the data, uncomment the following block
'''
plt.plot(time_data, voltage_data)
plt.title(plot_title)
plt.xlabel('Time')
plt.ylabel('Normalized Voltage')
plt.show()
'''
