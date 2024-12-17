# Code to collapse the trace in 1-D 
#%%
import os
import numpy as np
from astropy.io import fits
import glob
from pathlib import Path
from matplotlib import pyplot as plt
from scipy.signal import find_peaks

def collapse_trace_to_spec(file): 
    with fits.open(file) as hdul: 
        trace_data = hdul[0].data 
        # trace_header = hdul[0].header

        spectrum = np.sum(trace_data, axis = 0) 
        pxls = np.linspace(460, 1000, len(spectrum))

        return spectrum, pxls
    
def plot_spec(spectrum, pxls, xlim=None, ylim=None): 
    fig, ax = plt.subplots()
    ax.plot(pxls, spectrum)
    ax.set_xlim(xlim, ylim)
    ax.set_xlabel('Pixel Number')
    ax.set_ylabel('Flux [arbitrary units]')


    
#%% 
file = 'HW7/wavelength_cal/arc0245.a.fits' 
spectrum, pxls = collapse_trace_to_spec(file)
# %%
plot_spec(spectrum, pxls, 500,1000)
#%%
peaks, _  = find_peaks(spectrum, height = 3*np.mean(spectrum))
plt.plot(pxls, spectrum, )
plt.plot(pxls[peaks], spectrum[peaks], "x")
plt.xlim(700, 900)
#plt.yscale('log')
plt.show()
# %%
print(pxls[peaks], spectrum[peaks])
# %%

#%%
# Over plot wl solution 
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
#%%
wl_solution = np.array([[731.31964809, 0.763772070],
                        [757.18475073, 0.96604350]])

# Create the plot
fig, ax1 = plt.subplots(figsize=(8, 6))

# Plot flux vs pixel numbers
ax1.plot(pxls, spectrum, label='Spectrum', color='tab:blue')
ax1.plot(pxls[peaks], spectrum[peaks], "x", color='red', label='Peaks')
ax1.set_xlabel("Pixel Number")
ax1.set_ylabel("Flux")
ax1.set_xlim(700, 900)
ax1.set_title("1D Spectrum with Pixel and Wavelength Axes")

# Add a secondary x-axis for wavelengths 
ax2 = ax1.secondary_xaxis('top')
for index in range(len(wl_solution[:][0])):
    ax1.axvline(wl_solution[index][0], 0, 3e6, color = 'k', ls = '--')
    ax2.set_xticks(wl_solution[index][1])
    ax2.set_xticklabels(str(wl_solution[index][1]))
    ax2.set_xlabel('Wavelength')

# Add legend and show the plot
fig.tight_layout()
plt.legend()
plt.show()

# Print pixel and flux at the peaks
print("Pixel numbers of peaks:", pxls[peaks])
print("Flux at peaks:", spectrum[peaks])
# %%
# Wavelength solution with two pixel-wavelength pairs
wl_solution = np.array([[731.31964809, 0.763772070],
                         [757.18475073, 0.96604350]])

# Create the plot
fig, ax1 = plt.subplots(figsize=(8, 6))

# Plot flux vs pixel numbers
ax1.plot(pxls, spectrum, label='Spectrum', color='tab:blue')
ax1.plot(pxls[peaks], spectrum[peaks], "x", color='red', label='Peaks')
ax1.set_xlabel("Pixel Number")
ax1.set_ylabel("Flux")
ax1.set_xlim(700, 900)
ax1.set_title("1D Spectrum with Pixel and Wavelength Axes")

# Add a secondary x-axis for wavelengths
# Mapping the pixel values to wavelengths from the wl_solution
def pixel_to_wavelength(pixels):
    """Directly map pixel values to wavelengths using the wl_solution array."""
    pixel_values = wl_solution[:, 0]  # Get pixel values from the solution
    wavelength_values = wl_solution[:, 1]  # Get corresponding wavelength values
    
    # Interpolate or directly match the pixel values with their corresponding wavelengths
    wavelengths = np.interp(pixels, pixel_values, wavelength_values)
    return wavelengths

# Set up the secondary x-axis for wavelength
ax2 = ax1.secondary_xaxis('top')
ax2.set_xlabel("Wavelength [nm]")

# Use the pixel-to-wavelength function to set the ticks on the secondary x-axis
wavelengths_at_px = pixel_to_wavelength(pxls)  # Get wavelengths for each pixel
ax2.set_xticks(pxls)  # Set the pixel positions as ticks on the secondary axis
ax2.set_xticklabels([f"{wavelength:.2f}" for wavelength in wavelengths_at_px])  # Set the wavelength values as labels

# Add legend and show the plot
fig.tight_layout()
plt.legend()
plt.show()

# Print pixel and flux at the peaks
print("Pixel numbers of peaks:", pxls[peaks])
print("Flux at peaks:", spectrum[peaks])

# %%
