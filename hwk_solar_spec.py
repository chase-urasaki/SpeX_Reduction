import os
import numpy as np
import ipdb
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

solar_mags = {#abs mag of sun in several filters
    #from
    #Willmer, Christopher NA. "The absolute magnitude of the sun in several filters." The Astrophysical Journal Supplement Series 236.2 (2018): 47.
    #Table 3
    "U": {
        "Vega": 5.61,
        "AB": 6.33,
    },
    "B": {
        "Vega": 5.44,
        "AB": 5.31,
    },
    "V": {
        "Vega": 4.81,
        "AB": 4.80,
    },
    "R": {
        "Vega": 4.43,
        "AB": 4.60,
    },
    "I": {
        "Vega": 4.10,
        "AB": 4.51,
    }
}


magnitude_system = {
    "U": {
        "lambda_c": 0.36,
        "dl_over_l": 0.15,
        "jy": 1810
    },
    "B": {
        "lambda_c": 0.44,
        "dl_over_l": 0.22,
        "jy": 4260
    },
    "V": {
        "lambda_c": 0.55,
        "dl_over_l": 0.16,
        "jy": 3640
    },
    "R": {
        "lambda_c": 0.64,
        "dl_over_l": 0.23,
        "jy": 3080
    },
    "I": {
        "lambda_c": 0.79,
        "dl_over_l": 0.19,
        "jy": 2550
    }
}

def mag_to_flux(total_mag=None, filter_choice=None):
    fc                     = filter_choice
    jansky                 = 1.51e7*magnitude_system[fc]['dl_over_l'] #photons sec^-1 m^-2
    flux                   = magnitude_system[fc]['jy']* 10**(-0.4 * total_mag)
    photons_per_m2_per_s   = flux * jansky
    return photons_per_m2_per_s


#constants
h = 6.62607015e-34  # joules seconds
c = 2.99792458e8  # meters per second
rsun = 6.957*10**10 #rad of sun, cm
pc   = 3.086*10**18 #1 pc in cm

def blackbodyspectrum(wave, temp):
    """wave: wavelengths in ANGSTROM
       temp: temperature in KELVIN
       returns the blackbody function in ergs/cm2/s/angstrom
       the  flux at the Earth can be found by multiplying by (R/D)**2"""
    w = wave / 1.E8                              # Angstroms to cm
    #constants appropriate to cgs units.
    c1 =  3.7417749*10**-5                # =2*!DPI*h*c*c
    c2 =  1.4387687                  # =h*c/k
    val =  c2/w/temp
    #mstr = machar(double = (size(val,/type) EQ 5) )  #Get machine precision
    #good = where( val LT alog(mstr.xmax), Ngood )    #Avoid floating underflow

    bbflux =  c1 / ( w**5 * ( np.exp( val)-1. ) )

    return bbflux*1.E-8              # Convert to ergs/cm2/s/A

class Blackbody:
    """a class used to get blackbody spectra of star-sized blackbodies
       teff in kelvin, solar rad in solar radius units, dist_pc in parsec
       wavelengths is in angstroms

       The output Watts/m2/nm and photons/m2/s/nm have been
       checked by comparing with independent calculations online,
       including using zeropoint magnitudes and found to be accurate"""
    def __init__(self, teff = 5772,
                 solar_rad = 1.0, dist_pc = 10,
                 wavelengths = np.arange(3000, 25000)):
        self.teff               = teff
        self.solar_rad          = solar_rad
        self.dist_pc            = dist_pc
        self.wavelengths        = wavelengths
        self.wavelengths_nm     = self.wavelengths/10
        #blackbodyspectrum is the equivalent blackbodyspectrum
        self.blackbodyspectrum  = blackbodyspectrum(self.wavelengths,
                                                       self.teff)
        self.cgsbbfluxatearth= (self.blackbodyspectrum*
                            (self.solar_rad*rsun)**2/
                            (self.dist_pc*pc)**2)
        self.watts_per_m2_per_nm = self.cgsbbfluxatearth/100
        #The weird thing below is hc/lambda with lambda in m
        self.photons_per_m2_per_s_per_nm = (self.watts_per_m2_per_nm*
                                      self.wavelengths_nm*(10**-9)/
                                      (6.626*3*10**(-34+8)) )

class SolarSpec:
    def __init__(self):
        self.refspec = np.genfromtxt('astm_solar_reference_spectrum.txt',
                                    skip_header=1)

        self.wavelengths_nm         = self.refspec[:,0]
        self.watts_per_m2_per_nm    = self.refspec[:,1]
        self.photons_per_m2_per_s_per_nm  = self.watts_per_m2_per_nm/(h * c / (self.wavelengths_nm / 1e9))

    def get_photons_per_m2_per_s(self, dist_pc=1/206265, filterfunc = None):
        #distance in parsecs
        if filterfunc is None:
            print("Need a filter function")
            return np.nan
        else:
            filt_trans = filterfunc(self.wavelengths_nm)
            wav_bin_width = np.diff(self.wavelengths_nm)
            #Should do trapezoidal rule, just addingthe last value
            #CLOSE ENOUGH
            wav_bin_width_eqsize = np.append(wav_bin_width, wav_bin_width[-1])
            return np.sum(filt_trans*
                          wav_bin_width_eqsize*
                          self.photons_per_m2_per_s_per_nm *
                          #scale to distance in pc
                          1/(206265*dist_pc)**2)

class Filter:
    def __init__(self, name, datafile):
        self.name = name
        self.data = np.loadtxt(datafile)
        self.wavelength_nm = self.data[:,0]/10
        self.transmission = self.data[:,1]
        self.filterfunc = interp1d(self.wavelength_nm, self.transmission, kind='linear', bounds_error=False, fill_value=0)

if __name__ == "__main__":
    solarspec = SolarSpec()
    bb_spec = Blackbody(teff = 5772, solar_rad=1.0, dist_pc = 1/206265)


    print("""a. Plot the measured solar spectrum against that expected from a solar-sized blackbody (Teff = 5772 K) at a distance of 1 AU. Confirm this is a decent approximation to the solar spectrum.""")

    plt.plot(bb_spec.wavelengths_nm, bb_spec.watts_per_m2_per_nm,
             label = 'blackbody',
             color = 'Red')
    plt.plot(solarspec.wavelengths_nm,
             solarspec.watts_per_m2_per_nm,
             color = 'Black', alpha = 0.5,
             label = 'ASTM E-490 solar irradiance spectrum')
    plt.xlim((300, 1000))
    plt.ylim((0, 2.5))
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Spectral flux density [W/m$^2$/nm]')
    plt.legend(prop={'size': 10})
    plt.tight_layout()
    plt.savefig('p3a.png', dpi = 300)
    os.system('open '+'p3a.png')

    print("Part 2")
    vfilt = Filter('V', "./filts/Generic_Bessell.V.dat")

    for filtname in ['U', 'B', 'V', 'R', 'I']:
        #create a filter
        filt = Filter(filtname, "./filts/Generic_Bessell."+filtname+".dat")
        print(filtname, " Filter:")
        #do the zeropoint calculation
        print("\tPhotons per second from filter zeropoint:",
              "%.2e"%mag_to_flux(total_mag=solar_mags[filtname]['Vega'],
                          filter_choice=filtname))
        #do the calculation from the spectrum
        print("\tPhotons per second from spectrum at 10 pc:",
              "%.2e"%solarspec.get_photons_per_m2_per_s(dist_pc=10,
                   filterfunc = filt.filterfunc))
