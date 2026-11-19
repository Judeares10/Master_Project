from netCDF4 import Dataset
import numpy as np

# -------------------------------------------------------
# Convert:
#   photon s^-1 cm^-2 nm^-1
# to:
#   erg s^-1 cm^-1 ster^-1
#
# Input file:
#   hnu_uv_sun.ncdf
#
# Output file:
#   uv_sun_wasp39_units.ncdf
# -------------------------------------------------------

# Physical constants (cgs)
h = 6.62607015e-27      # Planck constant [erg s]
c = 2.99792458e10       # Speed of light [cm s^-1]

# -------------------------------------------------------
# Load UV Sun file
# -------------------------------------------------------
infile = "/home/jh1639/atmo/tutorial/Wasp39b/hnu_uv_sun.ncdf"
outfile = "/home/jh1639/atmo/tutorial/kurucz_sun_sun_units.ncdf"

ds = Dataset(infile)

# Input variables
lambda_nm = ds.variables["Wavelength"][:]   # nm
phot_flux = ds.variables["Hnu"][:]          # photon s^-1 cm^-2 nm^-1

ds.close()

# -------------------------------------------------------
# Convert wavelength to cm
# -------------------------------------------------------
lambda_cm = lambda_nm * 1e-7

# -------------------------------------------------------
# Photon energy:
# E = h*c/lambda
# -------------------------------------------------------
E_photon = (h * c) / lambda_cm   # erg

# -------------------------------------------------------
# Convert photon flux -> energy flux per nm
#
# erg s^-1 cm^-2 nm^-1
# -------------------------------------------------------
F_lambda = phot_flux * E_photon

# -------------------------------------------------------
# Convert wavelength -> wavenumber
#
# nu_tilde = 1/lambda(cm)
#
# output units need:
# erg s^-1 cm^-1 ster^-1
#
# Use:
# F_nu_tilde = F_lambda * |d(lambda_nm)/d(nu_tilde)|
#
# d(lambda_nm)/d(nu_tilde) = lambda_nm^2 / 1e7
# -------------------------------------------------------
nu_cm1 = 1.0 / lambda_cm

F_wavenumber = F_lambda * (lambda_nm**2 / 1e7)

# -------------------------------------------------------
# Convert flux -> intensity per steradian
#
# Assuming isotropic hemisphere:
# F = pi * I
# so:
# I = F/pi
# -------------------------------------------------------
I_wavenumber = F_wavenumber / np.pi

# -------------------------------------------------------
# Sort into increasing wavenumber
# -------------------------------------------------------
idx = np.argsort(nu_cm1)

nu_cm1 = nu_cm1[idx]
I_wavenumber = I_wavenumber[idx]

# -------------------------------------------------------
# Save new NetCDF file
# -------------------------------------------------------
out = Dataset(outfile, "w")

n = len(nu_cm1)

out.createDimension("nfreq", n)

nu_var = out.createVariable("nu", "f8", ("nfreq",))
hnu_var = out.createVariable("hnu", "f8", ("nfreq",))

nu_var[:] = nu_cm1
hnu_var[:] = I_wavenumber

nu_var.units = "cm-1"
hnu_var.units = "erg s-1 cm-1 ster-1"

out.close()

print(f"Saved converted file: {outfile}")