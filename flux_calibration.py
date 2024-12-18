#%% 
import numpy as np 
from sky_subtraction import sky_subtraction
from matplotlib import pyplot as plt
from scipy.constants import h, c, k

#%%
directory = 'HW7/reduced_standards'

# need sky subtracted standards as well 
reduced_standards_dir = 'HW7/reduced_standards'
subtraction = sky_subtraction(reduced_standards_dir)
diff, summed_diff = subtraction.nod_pair('HW7/reduced_standards/spc0228.a_trace.fits', 'HW7/reduced_standards/spc0229.b_trace.fits')
plt.plot(diff[:, 435] - np.mean(diff[:,435]))
# get the averaged spectra for the standard 

# %%
# divide by 10 seconds 



# zero_point_flux_J = 3.129e-13  # Zero-point flux for J-band in W/m^2/Hz (MKO system)
wavelength_range = (1.17, 1.33)  # J-band range in microns
J_magnitude = 6.92  # J-band magnitude of the standard star
T_eff = 9700  # Effective temperature of the standard star (in K)

# Wavelength grid in microns
wl = 
wavelengths = np.linspace(0.5, 2.9, 549)  # Example grid from 1.0 to 1.4 µm

# 1. Blackbody flux (spectral energy distribution of the standard star)
def blackbody(wavelength_micron, T):
    wavelength_m = wavelength_micron * 1e-6  # Convert microns to meters
    flux = (2 * h * c**2) / (wavelength_m**5) * (1 / (np.exp(h * c / (wavelength_m * k * T)) - 1))
    return flux  # In W/m^2/m

# Blackbody flux for the standard star
bb_flux = blackbody(wavelengths, T_eff)  # In W/m^2/m

# 2. Convert zero-point flux to flux density (F_lambda)
# F_nu to F_lambda: F_lambda = F_nu * c / lambda^2
j_band_central_wavelength = 1.25  # Central wavelength of J-band in microns
zero_point_flux_J_lambda = zero_point_flux_J * c / (j_band_central_wavelength * 1e-6)**2  # W/m^2/m

# Scale the blackbody flux in the J-band
j_band_indices = (wavelengths >= 1.17) & (wavelengths <= 1.33)
average_bb_flux_j_band = np.mean(bb_flux[j_band_indices])
scaling_factor = zero_point_flux_J_lambda / average_bb_flux_j_band

# Scale the entire blackbody flux
scaled_bb_flux = bb_flux * scaling_factor