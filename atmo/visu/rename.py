from netCDF4 import *
import numpy as np

h = 6.62607015e-27   # erg s
c = 2.99792458e10    # cm/s


def plot_flux(file=""):
    
    ncdfFile = Dataset(file,'r')
    vars = ncdfFile.variables

    print(vars)
    
    


   

    wavelength_cm = np.asarray(wavelength_nm) * 1e-7

    photon_energy = h * c / wavelength_cm

def photon_erg(file=""):
    
    ncdfFile = Dataset(file,'r')
    vars = ncdfFile.variables

    Wavelength = vars["Wavelength"][:]
    Hnu = vars["Hnu"][:]
    
    wavelength_cm = np.asarray(Wavelength) * 1e-7
    Energy = h * c / wavelength_cm
    
    flux = np.asarray(Hnu) * Energy

    flux = flux * 1e7
    print(flux)
    




    