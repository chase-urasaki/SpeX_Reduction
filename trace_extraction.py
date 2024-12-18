# Code for pulling out order traces 
#%%
import os
import numpy as np
from astropy.io import fits
import glob
from pathlib import Path
from matplotlib import pyplot as plt

def ds9(data):
    """
    Shortcut to display a FITS array in DS9 (external FITS viewer).
    """
    fits.writeto('temp.fits', np.array(data), overwrite=True)
    os.system('ds9 temp.fits -zoom to fit &')
    return

def write_fits(data, name, target_dir, header=None):
    """
    Write data to a FITS file in the specified directory.
    """
    full_path = os.path.join(target_dir, name)
    hdu = fits.PrimaryHDU(data, header=header)
    hdu.writeto(full_path, overwrite=True, output_verify='ignore')

class TraceExtraction_old_Spex:
    def __init__(self, directory):
        """
        Initialize the trace extraction class with a directory containing FITS files.
        """
        self.directory = directory

        trace_files = glob.glob(os.path.join(self.directory, "*trace*.fits"))
        for trace_file in trace_files:
            try:
                os.remove(trace_file)
                print(f"Removed existing trace file: {trace_file}")
            except Exception as e:
                print(f"Error removing file {trace_file}: {e}")
                
        self.file_list = glob.glob(os.path.join(directory, "*.fits"))

        if len(self.file_list) == 0:
            raise ValueError("No FITS files found in the specified directory.")

    def manual_extraction(self, x1, x2, y1, y2):
        """
        Manually extract a trace from a region defined by coordinates (x1, x2, y1, y2).
        Saves the extracted trace as a new FITS file with updated headers.
        Skips files that have already been processed.
        """


        for file in self.file_list:
            try:
                base_name = Path(file).stem
                output_name = f"{base_name}_trace.fits"
                output_path = os.path.join(self.directory, output_name)

                # Skip if the trace file already exists
                if os.path.exists(output_path):
                    print(f"Trace file already exists, skipping: {output_name}")
                    continue

                with fits.open(file) as hdul:
                    data = hdul[0].data
                    header = hdul[0].header

                    # Extract the trace from the specified region
                    cut_data = data[x1:x2, :]
                    trace = cut_data[:, y1:y2]

                    # Update header history
                    header.add_history('Order traced manually')

                    # Write the trace to a new file, overwriting if necessary
                    write_fits(trace, output_name, self.directory, header=header)
                    print(f"Extracted trace saved to {output_name}")
            except Exception as e:
                print(f"Error processing file {file}: {e}")


if __name__ == "__main__":
    x1=390; x2=482; y1=550; y2=950

    science_directory = "HW7/reduced_science"  # Update with the correct path\
    standard_directory = 'HW7/reduced_standards'
    sci_extractor = TraceExtraction_old_Spex(science_directory)
    sci_extractor.manual_extraction(x1, x2, y1, y2)

    stand_extractor = TraceExtraction_old_Spex(standard_directory)
    stand_extractor.manual_extraction(x1, x2, y1, y2)

    cal_extractor = TraceExtraction_old_Spex('HW7/wavelength_cal')
    cal_extractor.manual_extraction(x1,x2,y1,y2)
# %%
