# Import modules 
import numpy as np 
from astropy.io import fits

#%% 
# Pre-processing
    # CAll routine to sort all files in a directory \
def sort_files(directory):
    


# Call routine to get flats and median combine
def get_flats():
    for flat_frame in reduction.ini 
    return 
    #for file in directory open object header and search for 'Flat' 

# Call routine to get biases and median combine 
def get_biases(): 

# Call routine to get darks and median combine (also determine if you're doing the "fast or slow method")
def get_darks():
    
# bad pixel masking/Cosmic ray correction? 

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
    sort_files
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
    