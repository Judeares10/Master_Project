from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

h = 6.62607015e-27 
c = 2.99792458e10


def plot_flux(t="",files=[]):
    vars=[0,0]

    plt.figure(figsize=(10,6))
    
    

    for i in files:

        if i.endswith(".ncdf"):
        
            ds = Dataset(i)
            vars = list(ds.variables.keys())     
        

        
            if vars[0] == "Wavelength":
                # Extracts wavelength with units of nm
                lambda_nm = ds.variables["Wavelength"][:]   
            
                # Extracts photon flux with units of photon s^-1 cm^-2 nm^-1
                phot_flux = ds.variables["Hnu"][:]
            
                # Converts wavelength units to cm
                lambda_cm = lambda_nm * 1e-7

                # Finds energy of photon at wavelength 
                E_photon = (h * c) / lambda_cm
            
                # Converts photon flux to energy flux per nm 
                F_lambda = phot_flux * E_photon

                #Converts to wave number
                nu_cm1 = 1.0 / lambda_cm
            

                F_wavenumber = F_lambda * (lambda_nm**2 / 1e7)
        
                I_wavenumber = F_wavenumber / np.pi

                idx = np.argsort(nu_cm1)

                nu_cm1 = nu_cm1[idx]
            
                hnu = I_wavenumber[idx]
            
                # Convert wavenumber to wavelength with cm to um unit conversion
                lam = 1e4 / nu_cm1
        
            elif vars[0] == "nu":
            
                ds = Dataset(i)

                nu = ds.variables["nu"][:]         
                hnu = ds.variables["hnu"][:] 

                # Convert wavenumber to wavelength with cm to um unit conversion
                lam = 1e4 / nu
            
            hnu = hnu*1e3*np.pi


        elif i.endswith(".star"):
            
            data = np.loadtxt(i, comments="!")
            c_m = 2.99792458e8
            lam = data[:, 0]   # wavelength [um]
            hnu = data[:, 1]   # Inu [erg s^-1 cm^-2 Hz^-1 sr^-1]

            lam_m = lam * 1e-6          # um -> m

            hnu = hnu * c_m / lam_m**2 * 1e-9
            hnu = hnu * np.pi

        
        # Converts erg s^-1 cm^-2 cm^-1 ster^-1 to W m^-2 um^-1
        
        plt.semilogy(lam, hnu, label=((i.split("/"))[-1]),linewidth=0.5)
    
    plt.xlabel(r"Wavelength ($\mu$m)")
    plt.ylabel(r"Specific Stellar Flux (W m$^{-2}$ $\mu$m$^{-1}$)")
    plt.title(t)
    plt.legend()
    plt.grid(True, which="both", alpha=0.3)
    plt.xlim(0,1)
    plt.tight_layout()
    plt.savefig((t) +".png", dpi=200,)
    plt.show
