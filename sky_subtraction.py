# Code for sky subtraction via AB nod pair 
#%%
# Import modules 
import numpy as np 
from astropy.io import fits
import os
import glob
from pathlib import Path
from matplotlib import pyplot as plt
from astropy.stats import sigma_clip

#%%
def ds9(a):
    """ ds9 shrotcut from Mike Bottom"""
    fits.writeto('temp.fits', np.array(a), overwrite=True)
    os.system('ds9 temp.fits -zoom to fit &')
    return

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
        
        with fits.open(a_trace_fits) as hdul_a: 
             a_data = hdul_a[0].data 

             filtered_a_data = sigma_clip(a_data, sigma = 10)

             a_header = hdul_a[0].header 

             #collapsed_a_data = np.sum(filtered_a_data, axis=0)

        with fits.open(b_trace_fits) as hdul_b: 
             b_data = hdul_b[0].data 
             b_header = hdul_b[0].header 

             filtered_b_data = sigma_clip(b_data, sigma = 10)
             #collapsed_b_data = np.sum(filtered_b_data, axis=0)

        # Pairwise subtraction
        diff = np.subtract(filtered_a_data, filtered_b_data)
        ds9(diff)

        #cut the diff im in two spectra 
        split = diff.shape[0] // 2 

        a_spec = diff[split:, :] 
        b_spec = diff[:split, :]
       
        return a_spec, b_spec
    
    # Write a function for the other way later

#%%
if __name__ == "__main__":
    reduced_science_dir = 'HW7/reduced_science'
    subtraction = sky_subtraction(reduced_science_dir)
    spc0224a_trace, spc0225b_trace = subtraction.nod_pair('HW7/reduced_science/spc0224.a_trace.fits', 'HW7/reduced_science/spc0225.b_trace.fits')
    spc0227a_trace, spc0226b_trace = subtraction.nod_pair('HW7/reduced_science/spc0227.a_trace.fits', 'HW7/reduced_science/spc0226.b_trace.fits')
    
    def collapse(trace): 
        oneD = np.sum(trace, axis = 0) 
        return oneD

    spc0224a_spec = np.asarray(collapse(spc0224a_trace))
    spc0225b_spec = np.asarray(-1*collapse(spc0225b_trace))
    spc0227a_spec = np.asarray(collapse(spc0227a_trace))
    spc0226b_spec = np.asarray(-1*collapse(spc0226b_trace))
    plt.plot(spc0224a_spec, label = 'a')
    plt.plot(spc0225b_spec, label = 'b')
    plt.plot(spc0226b_spec, label = 'b2')
    plt.plot(spc0227a_spec, label = 'a2')
    plt.legend()

    #%%
    # Stack the spectra 
    stack = np.dstack((spc0224a_spec, spc0225b_spec, spc0226b_spec, spc0227a_spec))
    #%%
    np.shape(stack)
    #%%
    median_stack = np.median(stack, axis=2)[0]

    #shift all values to ensure they are above 0 
    #median_stack += np.abs(np.min(median_stack))
    #%%
    plt.plot(median_stack)


# %%

    reduced_standards_dir = 'HW7/reduced_standards'
    stand_subtraction = sky_subtraction(reduced_standards_dir)
        
    # %%
    spc0228a_trace, spc0229b_trace = subtraction.nod_pair('HW7/reduced_standards/spc0228.a_trace.fits',
                                                        'HW7/reduced_standards/spc0229.b_trace.fits')
    # %%
    spc0228a_spec = np.asarray(collapse(spc0228a_trace))
    spc0229b_spec = np.asarray(-1*collapse(spc0229b_trace))
    plt.plot(spc0228a_spec)
    plt.plot(spc0229b_spec)
    # %%
    standard_stac = np.dstack((spc0228a_spec, spc0229b_spec))
    standard_median_stack = np.median(standard_stac, axis = 2)[0]

    plt.plot(standard_median_stack)

#%% 
    # export the arrays 
    np.save('science_fluxes.npy', median_stack)
    np.save('standard_fluxes.npy', standard_median_stack)
# %%
