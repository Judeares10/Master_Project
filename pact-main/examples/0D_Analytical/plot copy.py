
from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt

folder = Path("/home/jh1639/pact-main/examples/0D/0D")

# Species you want to plot
species_to_plot = ["CH4", "S2", "H2O", "NH3"]

# File prefixes, pressures and output folders
pressure_groups = {
    "0DLow_out_": ("1e-3bar", "4b"),
    "0DHigh_out_": ("100bar", "4a")
}

for prefix, (pressure, output_folder) in pressure_groups.items():

    # Create output folder
    save_folder = folder / output_folder
    save_folder.mkdir(parents=True, exist_ok=True)

    for species in species_to_plot:

        times = []
        abundances = []

        for file in sorted(folder.glob(f"{prefix}*.dat")):

            with open(file, "r") as f:
                lines = f.readlines()

            # Extract integration time
            time = None

            for line in lines:
                if "Continuity equation integrated up to" in line:
                    match = re.search(r"=\s*([0-9.Ee+-]+)\s*s", line)

                    if match:
                        time = float(match.group(1))

                    break

            if time is None:
                print(f"Could not find time in {file.name}")
                continue

            # Find species column names
            header_index = None

            for i, line in enumerate(lines):
                if "height[km]" in line:
                    header_index = i
                    break

            if header_index is None:
                print(f"Could not find column header in {file.name}")
                continue

            species_names = lines[header_index].split()

            if species not in species_names:
                print(f"{species} not found in {file.name}")
                continue

            # Read atmospheric data row
            data = lines[header_index + 1].split()

            # Get the correct column for the species
            species_index = species_names.index(species)

            abundance = float(data[species_index])

            times.append(time)
            abundances.append(abundance)

        if len(times) < 2:
            print(f"Not enough data for {species} at {pressure}")
            continue

        # Sort results by integration time
        times = np.array(times)
        abundances = np.array(abundances)

        sort_index = np.argsort(times)

        times = times[sort_index]
        abundances = abundances[sort_index]

        abundances[0] = 0.0
        times = times[1:]
        abundances = abundances[1:]

        # Convert log10 abundances into mole fractions
        mole_fractions = 10**abundances

        # Plot abundance against time
        plt.figure(figsize=(8, 5))

        plt.plot(times,mole_fractions,color="black",linestyle="-",linewidth=1.5)

        plt.xscale("log")
        plt.yscale("log")
        plt.xlim(1e-16,1e16)

        plt.xlabel("Integration time (s)")
        plt.ylabel(f"{species} mole fraction")
        plt.title(f"{species} abundance evolution ({pressure})")

        plt.tight_layout()

        # Save individual species plot
        filename = save_folder / f"Abundance_{species}_{pressure}.png"

        plt.savefig(filename, dpi=300)
        plt.close()

        print(f"Saved {filename}")

