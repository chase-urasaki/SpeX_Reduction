#%%
import os
from astropy.io import fits 

#%%
def make_folders():
    """ Function that sorts a directory of flats, biases, darks, wavelength solutions, and science frames"""
    # Create folders: 
    os.makedirs('flats')
    os.makedirs('biases')
    os.makedirs('darks')
    os.makedirs('wavelength_cal')
    os.makedirs('science')


def sort_directory(directory):
    for filename in os.listdir(directory): 
        with fits.open(str(filename)) as hdul:
            header = 