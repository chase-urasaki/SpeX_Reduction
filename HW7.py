#%%
# Import modules 
import numpy as np 
from astropy.io import fits
import sys
import glob
import ccdproc 
from pathlib import Path
from matplotlib import pyplot as plt
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
def ds9(a):
    """ ds9 shrotcut from Mike Bottom"""
    fits.writeto('temp.fits', np.array(a), overwrite=True)
    os.system('ds9 temp.fits -zoom to fit &')
    return
#%%
# Combine flats by median
def flats_combine(flats_directory, counts_threshold = 20000): 
    # Get all FITS files in the directory
    file_list = glob.glob(os.path.join(flats_directory, "*.fits"))
    if len(file_list) == 0:
        return None

    # Read the data from all files into a list
    flat_data = []
    for file in file_list: 
        with fits.open(file) as hdul:
            data = hdul[0].data  # Assuming data is in the primary HDU
            header = hdul[0].header
            flat_data.append(data)

    # Stack the arrays and compute the median along the third axis
    flat_stack = np.stack(flat_data, axis=-1)  # Create a 3D array

    val_distribution = flat_stack.flatten()
    filtered_vals = val_distribution[val_distribution > counts_threshold]
    median_normed_flat = np.median(flat_stack, axis=-1) / np.mean(filtered_vals) # Median combine
    
    return median_normed_flat

def write_fits(data, name, target_dir, header = None): 
    full_path = os.path.join(target_dir, name)
    hdu = fits.PrimaryHDU(data)
    hdu.writeto(full_path, overwrite=True)
#%%
def flats_combine(flats_directory, counts_threshold=20000):
    """
    Combine flat field images by median after normalizing them.
    
    Parameters:
        flats_directory (str): Directory containing flat field FITS files.
        counts_threshold (int): Threshold for valid pixel values.
    
    Returns:
        np.ndarray: Median-combined and normalized flat field.
    """
    import os  # Ensure os is imported
    # Get all FITS files in the directory
    file_list = glob.glob(os.path.join(flats_directory, "*.fits"))
    if len(file_list) == 0:
        print("No FITS files found in the directory.")
        return None

    # Initialize a list to hold flat data
    flat_data = []

    # Read the data from all files into a list
    for file in file_list:
        with fits.open(file) as hdul:
            data = hdul[0].data  # Assuming data is in the primary HDU
            if data is None:
                print(f"File {file} has no data; skipping.")
                continue

            #scaled_data = data / 12

            flat_data.append(data)

    # Ensure all arrays have the same shape
    if not all(flat.shape == flat_data[0].shape for flat in flat_data):
        raise ValueError("Not all flats have the same dimensions. Check your input files.")

    # Stack the arrays and compute the median along the stack axis
    flat_stack = np.stack(flat_data, axis=-1)  # Stack along the third axis

    # Flatten and filter values by threshold
    val_distribution = flat_stack.flatten()
    filtered_vals = val_distribution[val_distribution > counts_threshold]
    if filtered_vals.size == 0:
        raise ValueError("No valid pixels found above the counts_threshold.")

    # Compute the median combined flat normalized exposure time and by the mean of filtered values
    median_normed_flat = np.median(flat_stack/12, axis=-1) / np.mean(filtered_vals)

    return median_normed_flat   
#%%
# Testing
flats_dir = 'HW7/flats'
master_cals_directory = 'HW7/master_cals'
median_normed_flat = flats_combine(flats_directory=flats_dir)
print("Median flat created with shape:", median_normed_flat.shape)
ds9(median_normed_flat)
#%%
write_fits(median_normed_flat, 'master_flat.fits', target_dir=master_cals_directory)
#%% 
def reduce(raw_directory, reduced_directory, master_cals_dir):

    # Load in master flat:
    master_flat_file = glob.glob(os.path.join(master_cals_dir, "master_flat*.fits"))
    if len(master_flat_file) == 0:
        print("Master flat not found!")
        return None

    with fits.open(master_flat_file[0]) as flat_hdu:
        master_flat_data = flat_hdu[0].data

    # Reduce raws
    raw_file_list = glob.glob(os.path.join(raw_directory, "*.fits"))
    if len(raw_file_list) == 0:
        print("No raw files found!")
        return None

    # Ensure reduced directory exists
    if not os.path.exists(reduced_directory):
        os.makedirs(reduced_directory)

    for file in raw_file_list:
        with fits.open(file) as hdul:
            raw_data = hdul[0].data
            header = hdul[0].header

        # Perform reduction
        reduced_data = raw_data / master_flat_data

        # Add history to header
        header.add_history('Reduced with master flat')

        # Define the output file path
        reduced_file_path = os.path.join(reduced_directory, os.path.basename(file))

        # Save to new file in reduced directory
        hdu = fits.PrimaryHDU(reduced_data, header=header)
        hdu.writeto(reduced_file_path, overwrite=True, output_verify='ignore')
#%%
science_directory = 'HW7/science'
reduced_science_directory = 'HW7/reduced_science'

master_cals_directory = 'HW7/master_cals'

reduce(science_directory, reduced_science_directory, master_cals_directory)
#%%
standards_directory = 'HW7/standards'
reduced_standards_directory = 'HW7/reduced_standards'
reduce(standards_directory, reduced_standards_directory, master_cals_directory)
# %%

# Now do the order tracing
class TraceExtraction_old_Spex:
    def __init__(self, directory):
        """
        Initialize the trace extraction class with a directory containing FITS files.
        """
        self.directory = directory

        trace_files = glob.glob(os.path.join(self.directory, "*trace*.fits"))
        for trace_file in trace_files:
            try:
                os.remove(trace_file)
                print(f"Removed existing trace file: {trace_file}")
            except Exception as e:
                print(f"Error removing file {trace_file}: {e}")
                
        self.file_list = glob.glob(os.path.join(directory, "*.fits"))

        if len(self.file_list) == 0:
            raise ValueError("No FITS files found in the specified directory.")

    def manual_extraction(self, x1, x2, y1, y2):
        """
        Manually extract a trace from a region defined by coordinates (x1, x2, y1, y2).
        Saves the extracted trace as a new FITS file with updated headers.
        Skips files that have already been processed.
        """


        for file in self.file_list:
            try:
                base_name = Path(file).stem
                output_name = f"{base_name}_trace.fits"
                output_path = os.path.join(self.directory, output_name)

                # Skip if the trace file already exists
                if os.path.exists(output_path):
                    print(f"Trace file already exists, skipping: {output_name}")
                    continue

                with fits.open(file) as hdul:
                    data = hdul[0].data
                    header = hdul[0].header

                    # Extract the trace from the specified region
                    cut_data = data[x1:x2, :]
                    trace = cut_data[:, y1:y2]

                    # Update header history
                    header.add_history('Order traced manually')

                    # Write the trace to a new file, overwriting if necessary
                    write_fits(trace, output_name, self.directory, header=header)
                    print(f"Extracted trace saved to {output_name}")
            except Exception as e:
                print(f"Error processing file {file}: {e}")

#%%

x1=390; x2=482; y1=490; y2=1000

science_directory = "HW7/reduced_science"  # Update with the correct path\
standard_directory = 'HW7/reduced_standards'
sci_extractor = TraceExtraction_old_Spex(science_directory)
sci_extractor.manual_extraction(x1, x2, y1, y2)

stand_extractor = TraceExtraction_old_Spex(standard_directory)
stand_extractor.manual_extraction(x1, x2, y1, y2)

# %%
# Wavelength Calibration
directory = 'HW7/wavelength_cal'
file = 'HW7/wavelength_cal/arc0245.a_trace.fits' 
x1=390; x2=482; y1=490; y2=1000

extractor = TraceExtraction_old_Spex(directory)
extractor.manual_extraction(x1, x2, y1, y2)
#%%
spectrum, pxls = collapse_trace_to_spec(file)
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
wl_matches = np.array([[514.10018553, 0.763772070],
                        [564.19294991, 0.91254710],
                        [578.21892393, 0.96604350],
                        [754.54545455, 1.694521],
                        [778.58998145, 1.7919521]])

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
#Interpolate the matches (and extrapolate) to find the wavelength solution
pixel_range = np.arange(460,1001,1)
#%%
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
spectrum, pixels = collapse_trace_to_spec(file)
# %%
plt.plot(np.polyval(coeffs, pixels), spectrum)

# %%
