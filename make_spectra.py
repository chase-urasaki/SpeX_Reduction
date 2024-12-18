# Code to collapse the trace in 1-D 
#%%
import os
import numpy as np
from astropy.io import fits
import glob
from pathlib import Path
from matplotlib import pyplot as plt
#import order_traces

def collapse_trace_to_spec(file, y1, y2): 
    with fits.open(file) as hdul: 
        trace_data = hdul[0].data 
        spectrum = np.sum(trace_data, axis = 0) 
        pxls = np.linspace(y1, y2, len(spectrum))

        return spectrum, pxls
    
def plot_spec(spectrum, pxls, xlim=None, ylim=None): 
    fig, ax = plt.subplots()
    ax.plot(pxls, spectrum)
    ax.set_xlim(xlim, ylim)
    ax.set_xlabel('Pixel Number')
    ax.set_ylabel('Flux [arbitrary units]')
# %%
