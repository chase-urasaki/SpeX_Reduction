# Code for sky subtraction via AB nod pair 
#%%
# Import modules 
import numpy as np 
from astropy.io import fits
import sys
import glob
from pathlib import Path
from matplotlib import pyplot as plt
#%%
class sky_subtraction: 
    def __init__(self, directory): 
        """
        Initialize the sky subtraction with a directory containing trace FITS files.
        """
        self.directory = directory

        self.file_list = glob.glob(os.path.join(directory, "*trace*.fits"))

        if len(self.file_list) == 0:
            raise ValueError("No Trace FITS files found in the specified directory.")
    
        print(self.file_list)

    def nod_pair(self, a_trace_fits, b_trace_fits):
        with fits.open(file) as hdul:
                    data = hdul[0].data
                    header = hdul[0].header
        
        with fits.open(a_trace_fits) as hdul_a: 
             a_data = hdul_a[0].data 
             a_header = hdul_a[0].header 

        with fits.open(b_trace_fits) as hdul_b: 
             b_data = hdul_b[0].data 
             b_header = hdul_b[0].header 

        
        

        
        return 
    
    # Write a function for the other way later

#%%
if __name__ == "__main__":
    reduced_science_dir = 'HW7/reduced_science'
    subtraction = sky_subtraction(reduced_science_dir)
# %%
