#%%# version using ccdproc 
from astropy.io import fits 
import numpy as np 
from pathlib import Path 
from ccdproc import ImageFileCollection

#%% 
flats_directory = 'HW7/flats'
science_directory = 'HW7/science'
standards_directory = 'HW7/standards'
wavelength_cal_directory = 'HW7/wavelength_cal'

# %%
flats_images = ImageFileCollection(flats_directory)
science_images = ImageFileCollection(science_directory)
standards_images = ImageFileCollection(standards_directory)
wavelength_cal_images = ImageFileCollection(wavelength_cal_directory) 
# %%
stuff = ImageFileCollection('HW7/flats')
# %%
image_name = 'flat0240.a.fits'
image_path = Path(flats_directory) / image_name
hdu_list = fits.open(image_path)
hdu_list =fits.open(image_path
            )
print(hdu_list[0].header)# %%

 
# %%
