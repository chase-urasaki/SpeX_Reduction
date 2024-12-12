# Import modules 
import numpy as np 
from astropy.io import fits

#%% 
# Pre-processing
    # CAll routine to sort all files in a directory 


# Call routine to get flats and median combine
def get_flats(): 
    #for file in directory open object header and search for 'Flat' 

# Call routine to get biases and median combine 

# Call routine to get darks and median combine (also determine if you're doing the "fast or slow method")

# bad pixel masking/Cosmic ray correction? 

# Wavelength Calibration/solution

# Nodding subtraction 

# call routine to perform regular sky subtraction after nodding subtraction

# Create trace profile 
    # think astropy can do this 