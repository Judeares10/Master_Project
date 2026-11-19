from netCDF4 import Dataset
import numpy as np

# Takes template profile for PACT and PT profile from ATMO and creates new input for PACT
def replace_pt_profile(Template_PACT="",PT_ATMO="",OUTPUT_NAME="",pressure_col="pressure[bar]",temperature_col="temperature[K]"):
	  
	  # Reads the content of the "template profile"
    with open(Template_PACT, "r") as f:
        lines = f.readlines()

    # Locates the coloumns for pressure and temperature
    header_idx = None
    for i, line in enumerate(lines):
        if pressure_col in line and temperature_col in line:
            header_idx = i
            break

    # Extract names and determine indices for pressure and temperature
    header = lines[header_idx].split()
    p_idx = header.index(pressure_col)
    t_idx = header.index(temperature_col)

		# Extracts data for each coloumn
    data_lines = lines[header_idx + 1:]

    # Extract PACT pressure profile
    pact_p = np.array([float(line.split()[p_idx]) for line in data_lines])

    # Reads the ATMO file and extracts pressure and temperature
    with Dataset(PT_ATMO) as nc:
        atmo_p = np.array(nc.variables["pressure"][:]) / 1e6  # bar
        atmo_T = np.array(nc.variables["temperature"][:])

    # Reorders profiles based on pressure to ensure low to high
    sort_idx = np.argsort(atmo_p)
    atmo_p = atmo_p[sort_idx]
    atmo_T = atmo_T[sort_idx]


    # Checks the profiles have same region
    if pact_p.min() < atmo_p.min() or pact_p.max() > atmo_p.max():
        print("WARNING: PACT pressure range exceeds ATMO range → edge extrapolation")

		# As pressure and temp spacing not the same, interpulates over the atmo profile to determine temps for PACT 
    new_T = np.interp(np.log10(pact_p),np.log10(atmo_p),atmo_T)

		# Writes the new pressure temp profile into .dat file with abundances from orginal PACT profile
    with open(OUTPUT_NAME, "w") as f:

        for line in lines[:header_idx]:
            f.write(line)

        f.write(lines[header_idx])

        for i, line in enumerate(data_lines):
            cols = line.split()
            cols[t_idx] = f"{new_T[i]:.8e}"
            f.write(" ".join(cols) + "\n")