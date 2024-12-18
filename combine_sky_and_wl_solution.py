#%%
import numpy as np
from matplotlib import pyplot as plt 
#%%
science_fluxes = np.load('science_fluxes.npy')
standard_fluxes = np.load('standard_fluxes.npy')
wl_solution = np.load('wavelength_solution.npy')

assert len(wl_solution) == len(standard_fluxes) and len(wl_solution) == len(science_fluxes)


# %%
science_spec = np.column_stack((wl_solution, science_fluxes))
standard_spec = np.column_stack((wl_solution, standard_fluxes))

standard_spec
# %%
np.save('science_spec.npy', science_spec)
np.save('standard_spec.npy', standard_spec)
# %%
