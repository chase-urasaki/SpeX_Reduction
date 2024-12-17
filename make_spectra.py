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
peaks, _  = find_peaks(spectrum, height = 3*np.median(spectrum))
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
                        [757.18475073, 0.91254710],
                        [765.1026393, 0.96604350],
                        [856.95014663, 1.694521],
                        [871.20234604, 1.7919521]])

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
for index in range(len(wl_solution)):
    ax1.axvline(wl_solution[index][0], 0, 3e6, color = 'k', ls = '--')
    ax2.set_xticks([wl_solution[index][1]])
    ax2.set_xticklabels([wl_solution[index][1]])
    ax2.set_xlabel('Wavelength')

# Add legend and show the plot
fig.tight_layout()
plt.legend()
plt.show()

# Print pixel and flux at the peaks
print("Pixel numbers of peaks:", pxls[peaks])
print("Flux at peaks:", spectrum[peaks])
# %%
