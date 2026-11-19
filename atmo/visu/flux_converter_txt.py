import numpy as np
import matplotlib.pyplot as plt
def flux_convert(input_file="",output_file=""):

    wavelength, flux_1au = np.loadtxt(input_file, unpack=True)

    distance_au=0.04828
    flux_planet = flux_1au * (1 / distance_au)**2

    np.savetxt(output_file, np.column_stack([wavelength, flux_planet]),fmt="%.6e", header="Wavelength[nm]  Flux[erg/s/cm^2/nm at 0.04828 AU]")

    print(f"Converted spectrum saved to {output_file}")

