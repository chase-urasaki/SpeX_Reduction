#%% 
import numpy as np
from matplotlib import pyplot as plt
from scipy.constants import h, c, k

def blackbody(T_eff, wavelengths):
    """
    Computes the blackbody spectrum for a given temperature and wavelengths.

    Parameters:
    T_eff (float): Effective temperature in Kelvin.
    wavelengths (np.ndarray): Array of wavelengths in microns.

    Returns:
    np.ndarray: Flux density in W/m^2/micron.
    """
    wavelengths_m = wavelengths * 1e-6  # Convert microns to meters
    bb_flux = (2 * h * c**2 / wavelengths_m**5) / (np.exp(h * c / (wavelengths_m * k * T_eff)) - 1)
    return bb_flux * 1e-6  # Convert from W/m^2/m to W/m^2/micron

# Load science spectrum
science_spec = np.load('science_spec.npy')  # Wavelength and flux

# Define the effective temperature of the standard star and its distance
T_eff = 9500  # Example value for the standard star in Kelvin
d_std = 120  # Distance to the standard star in parsecs

# Extract wavelengths and flux from the science spectrum
science_wavelengths = science_spec[:, 0]  # Wavelengths in microns
science_flux = science_spec[:, 1]  # Flux in arbitrary units

# Extract the observed flux of the standard star (from your data)
standard_spec = np.load('standard_spec.npy')  # Wavelength and observed flux for standard star
standard_wavelengths = standard_spec[:, 0]  # Wavelengths in microns
observed_flux = standard_spec[:, 1]  # Flux of the standard star

# Compute the blackbody flux over the standard star's wavelength range
blackbody_flux = blackbody(T_eff, standard_wavelengths)

# Flux Calibration: Calculate the calibration factor
calibration_factor = observed_flux / blackbody_flux  # Ratio of observed flux to blackbody flux

# Interpolate calibration factor to match the science spectrum wavelengths
from scipy.interpolate import interp1d
calibration_function = interp1d(standard_wavelengths, calibration_factor, kind='linear', fill_value="extrapolate")
calibration_factors = calibration_function(science_wavelengths)

# Apply flux calibration to the science spectrum
calibrated_flux = science_flux * calibration_factors

# Plot the science spectrum and calibrated spectrum
plt.figure(figsize=(10, 6))
#plt.plot(science_wavelengths, science_flux, label="Original Science Spectrum", alpha=0.7)
plt.plot(science_wavelengths, calibrated_flux, label="Flux Calibrated Spectrum", alpha=0.7)
plt.xlabel("Wavelength (microns)")
plt.ylabel("Flux (arbitrary units)")
plt.legend()
plt.title("Science Spectrum and Flux Calibrated Spectrum")
plt.show()

# Save the calibrated spectrum
np.save('flux_calibrated_science_spec.npy', np.column_stack((science_wavelengths, calibrated_flux)))


# %%
#%%
import numpy as np
from matplotlib import pyplot as plt
from scipy.constants import h, c, k

def blackbody(T_eff, wavelengths):
    """
    Computes the blackbody spectrum for a given temperature and wavelengths.

    Parameters:
    T_eff (float): Effective temperature in Kelvin.
    wavelengths (np.ndarray): Array of wavelengths in microns.

    Returns:
    np.ndarray: Blackbody flux density in W/m^2/micron.
    """
    wavelengths_m = wavelengths * 1e-6  # Convert microns to meters
    bb_flux = (2 * h * c**2 / wavelengths_m**5) / (np.exp(h * c / (wavelengths_m * k * T_eff)) - 1)
    return bb_flux * 1e-6  # Convert from W/m^2/m to W/m^2/micron

# Load spectra
science_spec = np.load('science_spec.npy')
standard_spec = np.load('standard_spec.npy')

# Define the effective temperature of the standard star
T_eff = 9500  # Example value in Kelvin

# Extract wavelengths and observed fluxes
wavelengths = standard_spec[:, 0]  # Wavelengths in microns
observed_flux = standard_spec[:, 1]  # Observed flux of the standard star

# Compute the blackbody flux over the full wavelength range
blackbody_flux = blackbody(T_eff, wavelengths)

# Compute the telluric correction factor
telluric_correction = blackbody_flux / observed_flux

# Apply telluric correction to the science spectrum
science_wavelengths = science_spec[:, 0]
science_flux = science_spec[:, 1]
corrected_flux = science_flux * telluric_correction

# Plot the results
plt.figure(figsize=(10, 6))

plt.subplot(2, 1, 1)
plt.plot(wavelengths, observed_flux, label="Standard Observed Spectrum", color="blue")
plt.plot(wavelengths, blackbody_flux, label=f"Blackbody (T_eff={T_eff} K)", color="orange")
plt.title("Standard Star Spectrum and Blackbody Curve")
plt.xlabel("Wavelength (microns)")
plt.ylabel("Flux (arbitrary units)")
plt.legend()

plt.subplot(2, 1, 2)
plt.plot(science_wavelengths, science_flux, label="Science Spectrum (Uncorrected)", color="red")
plt.plot(science_wavelengths, corrected_flux, label="Science Spectrum (Corrected)", color="green")
plt.title("Science Spectrum Before and After Telluric Correction")
plt.xlabel("Wavelength (microns)")
plt.ylabel("Flux (arbitrary units)")
plt.legend()

plt.tight_layout()
plt.show()
#%%
import numpy as np
from matplotlib import pyplot as plt
from scipy.constants import h, c, k

def blackbody(T_eff, wavelengths):
    """
    Computes the blackbody spectrum for a given temperature and wavelengths.

    Parameters:
    T_eff (float): Effective temperature in Kelvin.
    wavelengths (np.ndarray): Array of wavelengths in microns.

    Returns:
    np.ndarray: Flux density in W/m^2/micron.
    """
    wavelengths_m = wavelengths * 1e-6  # Convert microns to meters
    bb_flux = (2 * h * c**2 / wavelengths_m**5) / (np.exp(h * c / (wavelengths_m * k * T_eff)) - 1)
    return bb_flux * 1e-6  # Convert from W/m^2/m to W/m^2/micron

# Load science spectrum
science_spec = np.load('science_spec.npy')  # Wavelength and flux

# Define the effective temperature of the standard star and its distance
T_eff = 9500  # Example value for the standard star in Kelvin
d_std = 120  # Distance to the standard star in parsecs

# Extract wavelengths and flux from the science spectrum
science_wavelengths = science_spec[:, 0]  # Wavelengths in microns
science_flux = science_spec[:, 1]  # Flux in arbitrary units

# Extract the observed flux of the standard star (from your data)
standard_spec = np.load('standard_spec.npy')  # Wavelength and observed flux for standard star
standard_wavelengths = standard_spec[:, 0]  # Wavelengths in microns
observed_flux = standard_spec[:, 1]  # Flux of the standard star

# Compute the blackbody flux over the standard star's wavelength range
blackbody_flux = blackbody(T_eff, standard_wavelengths)

# Flux Calibration: Calculate the calibration factor
calibration_factor = observed_flux / blackbody_flux  # Ratio of observed flux to blackbody flux

# Interpolate calibration factor to match the science spectrum wavelengths
from scipy.interpolate import interp1d
calibration_function = interp1d(standard_wavelengths, calibration_factor, kind='linear', fill_value="extrapolate")
calibration_factors = calibration_function(science_wavelengths)

# Apply flux calibration to the science spectrum
calibrated_flux = science_flux * calibration_factors

# Plot the science spectrum and calibrated spectrum
plt.figure(figsize=(10, 6))
#plt.plot(science_wavelengths, science_flux, label="Original Science Spectrum", alpha=0.7)
plt.plot(science_wavelengths, calibrated_flux, label="Flux Calibrated Spectrum", alpha=0.7)
plt.xlabel("Wavelength (microns)")
plt.ylabel("Flux (arbitrary units)")
plt.legend()
plt.title("Science Spectrum and Flux Calibrated Spectrum")
plt.show()

# Save the calibrated spectrum
#np.save('flux_calibrated_science_spec.npy', np.column_stack((science_wavelengths, calibrated_flux)))

# %%
plt.plot(science_wavelengths, science_wavelengths*calibrated_flux)
# %%
