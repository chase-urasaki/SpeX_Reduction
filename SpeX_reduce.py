#%%
# Import modules 
import numpy as np 
from astropy.io import fits
import sys
import glob
import ccdproc 
from pathlib import Path
from matplotlib import pyplot as plt


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
            exposure_time = np.float64(header.get('ITIME') * header.get('CO_ADDS'))
            scaled_data = data / exposure_time
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
reduced_science_direcotry = 'HW7/reduced_science'

master_cals_directory = 'HW7/master_cals'

reduce(science_directory, reduced_science_direcotry, master_cals_directory)
#%%
standards_directory = 'HW7/standards'
reduced_standards_directory = 'HW7/reduced_standards'
reduce(standards_directory, reduced_standards_directory, master_cals_directory)
# %%
