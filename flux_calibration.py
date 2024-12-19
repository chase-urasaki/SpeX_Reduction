#%% 
import numpy as np 
from matplotlib import pyplot as plt
from scipy.constants import h, c, k

#%%
# Load in the spectra 
science_spec = np.load('science_spec.npy')
standard_spec = np.load('standard_spec.npy')

standard_spec[:,0]
#%%
def blackbody(T_eff, standard_spec): 
    """
    Computes the blackbody temperature of a star with T_eff  over 
    the wavelengths of the standard spectrum"""


def compute_response_function(standard_true, standard_observed):
    """
    Compute the response function using the observed and true spectrum of the standard star.

    Parameters:
    standard_observed (np.ndarray): Observed spectrum of the standard star.
    standard_true (np.ndarray): True flux values of the standard star (known).

    Returns:
    np.ndarray: Response function (true flux / observed flux).
    """
    #if np.any(standard_observed == 0):
    #    raise ValueError("Observed spectrum contains zero values.")
    # is this the correct thing to do?
    #standard_observed[standard_observed==0] = 1
    response = standard_true[:,1] / standard_observed

    return response 


def flux_calibrate(target_observed, response_function):
    """
    Apply the response function to flux calibrate the target spectrum.

    Parameters:
    target_observed (np.ndarray): Observed spectrum of the target star.
    response_function (np.ndarray): Response function (derived from standard star).

    Returns:
    np.ndarray: Flux-calibrated spectrum of the target star.
    """
    calibrated = target_observed[:,1] * response_function
    return calibrated

# flux calibrate system (eqns in homework) (maybe do this as a full function?)
#%% 
plt.plot(standard_spec[:,0], standard_spec[:,1])

# HD 0202025 from homework/simbad
flux_zp = 3.01e-9 # flux for Vega in AB in W/m^2/micron
J_band_filter_cut_on = 1.17 # microns 
J_band_fitler_cut_off = 1.33 # microns 

J_band_mag_std = 6.925 # J-band mag of the standard star 
mag_zp = 0.943 # ZP mag of vega in AB 

# Compute flux of standard star 
flux_std = flux_zp * 10**((mag_zp - J_band_mag_std) / 2.5) 

#%%
# Compute the flux through the wavelength range affected by the standard 
inst_flux_std = 0 

for idx in range(len(standard_spec)): 
    if J_band_filter_cut_on <= standard_spec[idx, 0] <= J_band_fitler_cut_off:
       inst_flux_std += standard_spec[idx,1]
       #print(standard_spec[idx,:])

#%%
inst_flux_std

#%%
inst_response = flux_std / inst_flux_std

science_spec[:,1] *= inst_response
#%
plt.plot(science_spec[:,0], science_spec[:,1])

# %%
