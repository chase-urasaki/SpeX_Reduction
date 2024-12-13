#%%
# Import modules 
import numpy as np 
from astropy.io import fits
import sys
import glob


#%%
def ds9(a):
    """ ds9 shrotcut from Mike Bottom"""
    fits.writeto('temp.fits', np.array(a), overwrite=True)
    os.system('ds9 temp.fits -zoom to fit &')
    return
#%%
# Combine flats by median
def flats_combine(flats_directory): 
    # Get all FITS files in the directory
    file_list = glob.glob(os.path.join(flats_directory, "*.fits"))
    if len(file_list) == 0:
        return None

    # Read the data from all files into a list
    flat_data = []
    for file in file_list: 
        with fits.open(file) as hdul:
            data = hdul[0].data  # Assuming data is in the primary HDU
            flat_data.append(data)

    # Stack the arrays and compute the median along the third axis
    flat_stack = np.stack(flat_data, axis=-1)  # Create a 3D array
    median_flat = np.median(flat_stack, axis=-1)  # Median combine

    return median_flat

# Testing
flats_dir = 'HW7/flats'
median_flat = flats_combine(flats_directory=flats_dir)
print("Median flat created with shape:", median_flat.shape)
ds9(median_flat)

#%%

#%%
# Call routine to get biases and median combine 
def get_biases(): 

# Call routine to get darks and median combine (also determine if you're doing the "fast or slow method")
def get_darks():
    
# bad pixel masking/Cosmic ray correction? 

def reduce(science_directory, standards_directory, wavelength_cal):
    # Reduce science frames 
    science_file_list = glob.glob(os.path.join(science_directory, "*.fits"))
    if len(science_file_list) == 0: 
        return None
    
    os.makedir('reduced_science')
    os.makedir('reduced_standards')
    os.makedir('reduced_wavelength_cal')

    


    # Reduce standards
     
    # Reduce wavelength cal
# Wavelength Calibration/solution
def wavelength_sol(): 
    
    # plot of wavelength solution 

    # Return accuracy 


# Nodding subtraction 
def nod_subtract():

# call routine to perform regular sky subtraction after nodding subtraction
def telluric_subtraction():
    #use molecfit

def flux_calibration(): 
    # Use the standard stars 

# Create trace profile 
    # think astropy can do this 
def get_traces():

if __name__ == "__main__":
    if len(sys.argv) <7: 
        print("Usage: flats_dir, bias_dir, darks_dir, wavelength_sol_dir, standards_dir, science_dir")
    
    flats_dir = str(sys.argv[1])
    bias_dir = str(sys.argv[2])
    darks_dir = str(sys.argv[3])
    wavelength_sol_dir = str(sys.argv[4])
    standards_dir = str(sys.argv[5])
    science_dir = str(sys.argv[6])

    flats_combine(flats_directory=flats_dir)

    flat_combine
    dark_combine
    bias_combine
    bad_pixel_mask
    cr_subtract 
    reduce
    apply_wavelength_sol
    nod_subtract
    sky_subtract
    get_traces
    