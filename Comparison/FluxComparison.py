
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf_file

# --------------------------------------------------
# File paths (EDIT THESE)
# --------------------------------------------------

star_file = r"/home/jh1639/pact-main/data/stars/wasp-39.star"

atmo_file = r"/home/jh1639/atmo/tutorial/Wasp39b_corrected/BT_Settke_Wasp39.ncdf"

third_file = r"/home/jh1639/atmo/tutorial/Wasp39b_corrected/hnu_uv_wasp39.ncdf"

combined_file = r"/home/jh1639/atmo/tutorial/Wasp39b_corrected/BT_Settke_Wasp39_combined.ncdf"

txt_file = r"/home/jh1639/VULCAN/atm/stellar_flux/sflux-wasp39-frances.txt"

# ==================================================
# 1. READ WASP-39b STELLAR SPECTRUM (.star)
# ==================================================

star_data = np.loadtxt(
    star_file,
    comments='!',
    skiprows=1
)

wlth = star_data[:, 0]
Inu = star_data[:, 1]


# ==================================================
# 2. READ ATMO NETCDF SPECTRUM
# ==================================================

with netcdf_file(atmo_file, 'r', mmap=False) as nc:

    nu = np.array(nc.variables['nu'][:])
    hnu = np.array(nc.variables['hnu'][:])

    # Save the original NetCDF structure and attributes
    dimensions = dict(nc.dimensions)
    global_attributes = dict(nc._attributes)

    variable_attributes = {}
    variable_dimensions = {}
    variable_types = {}

    for name, variable in nc.variables.items():
        variable_attributes[name] = dict(variable._attributes)
        variable_dimensions[name] = variable.dimensions
        variable_types[name] = variable.typecode()


# ==================================================
# 3. CONVERT ATMO SPECTRUM
# ==================================================

c = 2.99792458e10  # Speed of light [cm s-1]

# Convert wavenumber [cm-1] to wavelength [um]
wavelength = 1e4 / nu

# Convert spectral intensity per wavenumber
# to spectral intensity per frequency

Inu_atmo = hnu / c

# Convert specific intensity to hemispheric flux
Inu_atmo = Inu_atmo * np.pi * np.pi


# ==================================================
# 4. SORT ATMO WAVELENGTHS
# ==================================================

order = np.argsort(wavelength)

wavelength = wavelength[order]
Inu_atmo = Inu_atmo[order]


# ==================================================
# 5. READ THIRD NETCDF SPECTRUM
# ==================================================

with netcdf_file(third_file, 'r', mmap=False) as nc:

    third_wavelength = np.array(
        nc.variables['Wavelength'][:]
    )

    third_hnu = np.array(
        nc.variables['Hnu'][:]
    )


# ==================================================
# 6. CONVERT THIRD SPECTRUM
# ==================================================

h = 6.62607015e-27  # Planck constant [erg s]

# Convert wavelength from nm to microns
third_wavelength_um = third_wavelength / 1000.0

# Convert photon flux per wavelength to energy flux
# per frequency

third_flux = third_hnu * h * third_wavelength

# Keep your existing pi conversion
third_flux = third_flux * np.pi

# Sort wavelengths into ascending order

order = np.argsort(third_wavelength_um)

third_wavelength_um = third_wavelength_um[order]
third_flux = third_flux[order]


# ==================================================
# 7. READ ADDITIONAL TXT SPECTRUM
# ==================================================

txt_data = np.loadtxt(
    txt_file,
    comments="#",
    skiprows=1
)

txt_wavelength = txt_data[:, 0]
txt_flux = txt_data[:, 1]

# Convert wavelength from nm to microns
txt_wavelength_um = txt_wavelength / 1000.0

# Convert flux per nm to flux per micron
txt_flux_um = txt_flux * 1000.0

# Convert flux per micron to flux per frequency
txt_flux_nu = txt_flux_um * (txt_wavelength_um**2 / 2.99792458e14)

# Sort wavelengths
txt_order = np.argsort(txt_wavelength_um)

txt_wavelength_um = txt_wavelength_um[txt_order]
txt_flux_nu = txt_flux_nu[txt_order]
txt_flux_nu = txt_flux_nu

# ==================================================
# 7. COMBINE THE TWO NETCDF SPECTRA
# ==================================================

# Convert the third spectrum back to the original ATMO hnu convention
third_hnu_converted = third_flux * c / np.pi

# Read the original .star wavelength grid
star_data = np.loadtxt(
    "/home/jh1639/pact-main/data/stars/wasp-39.star",
    comments="!"
)

combined_wavelength = star_data[:, 0]
original_Inu = star_data[:, 1]

# Original ATMO wavelength grid
atmo_wavelength = 1e4 / nu

# Sort ATMO data by wavelength
atmo_order = np.argsort(atmo_wavelength)

atmo_wavelength = atmo_wavelength[atmo_order]
atmo_hnu = hnu[atmo_order]

# Sort UV data by wavelength
uv_order = np.argsort(third_wavelength_um)

uv_wavelength = third_wavelength_um[uv_order]
uv_hnu = third_hnu_converted[uv_order]

# Interpolate ATMO and UV spectra onto the original .star grid
atmo_interp = np.interp(
    combined_wavelength,
    atmo_wavelength,
    atmo_hnu,
    left=np.nan,
    right=np.nan
)

uv_interp = np.interp(
    combined_wavelength,
    uv_wavelength,
    uv_hnu,
    left=np.nan,
    right=np.nan
)

# Start with the original .star spectrum
combined_Inu = original_Inu.copy()

# Replace with ATMO values where available
mask_atmo = np.isfinite(atmo_interp)
combined_Inu[mask_atmo] = atmo_interp[mask_atmo] / c * np.pi

# Replace with UV values where available
mask_uv = np.isfinite(uv_interp)
combined_Inu[mask_uv] = uv_interp[mask_uv] / c * np.pi

# Convert combined spectrum back to hnu convention for plotting
combined_hnu = combined_Inu * c / np.pi

# Convert back to wavenumber
combined_nu = 1e4 / combined_wavelength

print("Spectrum combination complete")
print("Original .star points:", len(star_data))
print("ATMO points available:", np.count_nonzero(mask_atmo))
print("UV points available:", np.count_nonzero(mask_uv))
print("Total combined points:", len(combined_wavelength))


# ==================================================
# SAVE COMBINED SPECTRUM AS .STAR FILE
# ==================================================

with open(
    "/home/jh1639/pact-main/data/stars/wasp-39_atmo.star",
    "w"
) as f:

    f.write("!  wlth[um]      Inu[erg s-1 cm-2 Hz-1 sr-1]\n")

    for wl, inu in zip(combined_wavelength, combined_Inu):
        f.write(f"      {wl:.5f}    {inu:.4E}\n")

print("Combined .star file saved.")
# ==================================================
# 10. PLOT ALL FOUR SPECTRA
# ==================================================

Inu_combined = combined_Inu


plt.figure(figsize=(10, 6))

plt.plot(
    wlth,
    Inu,
    label='WASP-39b stellar spectrum',
    linewidth=0.5
)

plt.plot(
    wavelength,
    Inu_atmo,
    label='ATMO spectrum',
    linewidth=0.5
)

plt.plot(txt_wavelength_um,txt_flux_nu,label='Additional TXT spectrum',linewidth=0.5)

plt.plot(
    third_wavelength_um,
    third_flux,
    label='Third spectrum (converted from photons)',
    linewidth=0.5
)

plt.plot(
    combined_wavelength,
    Inu_combined,
    label='Combined ATMO + UV spectrum',
    linewidth=0.5
)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'Wavelength [$\mu$m]')

plt.ylabel(
    r'Spectral flux '
    r'[erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$]'
)

plt.title('WASP-39b Spectral Comparison')

plt.grid(True, which='both', alpha=0.3)
plt.legend()

plt.tight_layout()
plt.savefig("StellarFluxb.png", dpi=300)
plt.show()


# ==================================================
# 11. PLOT COMBINED SPECTRUM ON ITS OWN
# ==================================================

plt.figure(figsize=(10, 6))

plt.plot(
    combined_wavelength,
    Inu_combined,
    linewidth=0.7
)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'Wavelength [$\mu$m]')
plt.ylabel(
    r'Spectral flux '
    r'[erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$]'
)

plt.title('WASP-39b Combined Spectrum')

plt.grid(True, which='both', alpha=0.3)

plt.tight_layout()
plt.savefig("WASP39b_combined.png", dpi=300)
