# Code for pulling out order traces 
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

# Use the flat fields as a mask to find the traces 
