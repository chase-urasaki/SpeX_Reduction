# Code wavelength calibration 
#%% 
import numpy as np 
from astropy.io import fits
import sys
import glob
from pathlib import Path
from matplotlib import pyplot as plt
from scipy.signal import find_peaks
from make_spectra import collapse_trace_to_spec
from make_spectra import plot_spec
#from order_traces import order1
from scipy import interpolate
from trace_extraction import TraceExtraction_old_Spex
#%%
#%% 
x1=390; x2=482; y1=550; y2=950

file = 'HW7/wavelength_cal/arc0245.a_trace.fits' 
#%%
spectrum, pxls = collapse_trace_to_spec(file, y1, y2)
# %%
plot_spec(spectrum, pxls)
#%%
peaks, _  = find_peaks(spectrum, height = 3*np.median(spectrum))
plt.plot(pxls, spectrum, )
plt.plot(pxls[peaks], spectrum[peaks], "x")
plt.show()
# %%
print(pxls[peaks], spectrum[peaks])
#%%
wl_matches = np.array([[514.05228758, 0.763772070],
                        [564.16122004, 0.91254710],
                        [578.19172113, 0.96604350],
                        [754.5751634 , 1.694521],
                        [778.62745098, 1.7919521]])

# wl_matches = np.array([[731.31964809, 0.763772070],
#                         [757.18475073, 0.91254710],
#                         [765.1026393, 0.96604350],
#                         [856.95014663, 1.694521],
#                         [871.20234604, 1.7919521]])
#%%
# Create the plot
fig, ax1 = plt.subplots(figsize=(8, 6))

# Plot flux vs pixel numbers
ax1.plot(pxls, spectrum, label='Spectrum', color='tab:blue')
ax1.plot(pxls[peaks], spectrum[peaks], "x", color='red', label='Peaks')
ax1.set_xlabel("Pixel Number")
ax1.set_ylabel("Flux")
ax1.set_title("1D Spectrum with Pixel and Wavelength Axes")

# Add a secondary x-axis for wavelengths 
ax2 = ax1.secondary_xaxis('top')
for index in range(len(wl_matches)):
    ax1.axvline(wl_matches[index][0], 0, 3e6, color = 'k', ls = '--')
    ax2.set_xticks([wl_matches[index][1]])
    ax2.set_xticklabels([wl_matches[index][1]])
    ax2.set_xlabel('Wavelength')

# Add legend and show the plot
fig.tight_layout()
plt.legend()
plt.show()

# Print pixel and flux at the peaks
print("Pixel numbers of peaks:", pxls[peaks])
print("Flux at peaks:", spectrum[peaks])
# %%
# Now plot the pixel points and their wavelength solutions



# %%
#Interpolate the matches (and extrapolate) to find the wavelength solution
pixel_range = np.arange(y1,y2+1,1)
#%%
# %%
a, b = np.polyfit(wl_matches[:,0], wl_matches[:,1],1)
coeffs = a, b
fit = np.polyval(coeffs, wl_matches[:,0])
wl_solution = np.polyval(coeffs, pixel_range)
#%%
plt.scatter(wl_matches[:,0], wl_matches[:,1])
plt.plot(wl_matches[:,0], fit)
plt.plot(pixel_range, wl_solution)
plt.xlabel('Pixel Number')
plt.ylabel('Wavelength [microns]')
# %%
# Apply calibration to the traces and make the spectrum - test  
file = 'HW7/reduced_science/spc0224.a_trace.fits'
spectrum, pixels = collapse_trace_to_spec(file, y1, y2)
wl = np.polyval(coeffs, pixels)
# %%
plt.plot(wl, spectrum)
#%% 
np.save('wavelength_solution.npy', wl)
# %%