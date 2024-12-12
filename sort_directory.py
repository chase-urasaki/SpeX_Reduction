#%%
import os
from astropy.io import fits 

#%%
def sort_directory():
    """ Function that sorts a directory of flats, biases, darks, wavelength solutions, and science frames"""
    # Create folders: 
    