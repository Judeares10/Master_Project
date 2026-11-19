
import numpy as np
import pathlib
import os
from pathlib import Path

input_dir = pathlib.Path("photo_cross")
output_dir = pathlib.Path("PACT_photo")
output_dir.mkdir(exist_ok=True)

# PACT reaction label
pact_label = "A"

# PACT photodissociation reactions
# Format: parent species: [(reaction ID, branching index, products)]

reactions = {
    "H2O": [
        (1071, 1, "H + OH"),
        (1073, 2, "H2 + O_1"),
        (1075, 3, "H + H + O"),
    ],
    "CH4": [
        (1077, 1, "CH3 + H"),
        (1079, 2, "CH2_1 + H2"),
        (1081, 3, "CH2_1 + H + H"),
        (1083, 4, "CH + H2 + H"),
    ],
    "CH3": [
        (1085, 1, "CH + H2"),
        (1087, 2, "CH2_1 + H"),
    ],
    "CO": [
        (1089, 1, "C + O"),
    ],
    "H2": [
        (1091, 1, "H + H"),
    ],
    "C2H2": [
        (1093, 1, "C2H + H"),
    ],
    "CO2": [
        (1095, 1, "CO + O"),
        (1097, 2, "CO + O_1"),
    ],
    "C2H4": [
        (1099, 1, "C2H2 + H2"),
        (1101, 2, "C2H2 + H + H"),
        (1103, 3, "C2H3 + H"),
    ],
    "C2H6": [
        (1105, 1, "C2H4 + H2"),
        (1107, 2, "C2H4 + H + H"),
        (1109, 3, "C2H2 + H2 + H2"),
        (1111, 4, "CH4 + CH2_1"),
        (1113, 5, "CH3 + CH3"),
    ],
    "OH": [
        (1115, 1, "H + O"),
    ],
    "HCO": [
        (1117, 1, "H + CO"),
    ],
    "H2CO": [
        (1119, 1, "H2 + CO"),
        (1121, 2, "H + HCO"),
    ],
    "CH3OH": [
        (1123, 1, "CH3 + OH"),
        (1125, 2, "H2CO + H2"),
        (1127, 3, "CH3O + H"),
    ],
    "CH3CHO": [
        (1129, 1, "CH4 + CO"),
        (1131, 2, "CH3 + HCO"),
    ],
    "N2": [
        (1133, 1, "N + N"),
    ],
    "NH3": [
        (1135, 1, "NH2 + H"),
        (1137, 2, "NH + H + H"),
    ],
    "HCN": [
        (1139, 1, "H + CN"),
    ],
    "NO": [
        (1141, 1, "N + O"),
    ],
    "NO2": [
        (1143, 1, "NO + O"),
    ],
    "N2H4": [
        (1145, 1, "N2H3 + H"),
    ],
    "HNCO": [
        (1147, 1, "NH + CO"),
        (1149, 2, "H + NCO"),
    ],
    "HNO2": [
        (1151, 1, "NO + OH"),
    ],
    "HO2": [
        (1153, 1, "O + OH"),
    ],
    "SO": [
        (1155, 1, "S + O"),
    ],
    "SO2": [
        (1157, 1, "SO + O"),
        (1159, 2, "S + O2"),
    ],
    "SH": [
        (1161, 1, "H + S"),
    ],
    "H2S": [
        (1163, 1, "SH + H"),
    ],
    "S2": [
        (1165, 1, "S + S"),
    ],
    "COS": [
        (1167, 1, "CO + S"),
    ],
    "CS": [
        (1169, 1, "C + S"),
    ],
    "CS2": [
        (1171, 1, "CS + S"),
    ],
    "O2": [
        (1173, 1, "O + O"),
        (1175, 2, "O + O_1"),
    ],
    "C4H2": [
        (1177, 1, "C2H2 + C2"),
    ],
    "CH2": [
        (1179, 1, "CH + H"),
    ],
    "C6H6": [
        (1181, 1, "C6H5 + H"),
        (1183, 2, "C3H3 + C3H3"),
    ],
    "CH3SH": [
        (1185, 1, "CH3S + H"),
        (1187, 2, "CH3 + SH"),
    ],
    "S4": [
        (1189, 1, "S2 + S2"),
    ],
    "H2CS": [
        (1191, 1, "H2 + CS"),
    ],
}


# Convert each parent species separately

for species, channels in reactions.items():

    cross_file = Path(input_dir) / species / f"{species}_cross.csv"
    branch_file = Path(input_dir) / species / f"{species}_branch.csv"

    if not cross_file.exists():
        print(f"Missing cross-section file: {cross_file}")
        continue

    if not branch_file.exists():
        print(f"Missing branching-ratio file: {branch_file}")
        continue

    print(f"\nProcessing {species}")

    # Read VULCAN cross-section data
    cross_data = np.loadtxt(cross_file, delimiter=",", comments="#", ndmin=2)
    
    wavelength = cross_data[:, 0]
    absorption = cross_data[:, 1]
    dissociation = cross_data[:, 2]
    
    unique_wavelengths = np.unique(wavelength)
    
    absorption_new = []
    dissociation_new = []
    
    for wl in unique_wavelengths:
        idx = np.where(wavelength == wl)[0]
        best_idx = idx[np.argmax(absorption[idx])]
    
        absorption_new.append(absorption[best_idx])
        dissociation_new.append(dissociation[best_idx])
    
    wavelength = unique_wavelengths
    absorption = np.array(absorption_new)
    dissociation = np.array(dissociation_new)

    # Read VULCAN branching-ratio data
    branch_data = np.loadtxt(
        branch_file,
        delimiter=",",
        comments="#",
        ndmin=2
    )

    branch_wavelength = branch_data[:, 0]
    branch_ratios = branch_data[:, 1:]

    # Check wavelength ordering
    if np.any(np.diff(wavelength) < 0):
        raise ValueError(f"{species}: cross-section wavelengths not sorted")

    if np.any(np.diff(branch_wavelength) < 0):
        raise ValueError(f"{species}: branching wavelengths not sorted")

    # Check that all required branching columns exist
    max_index = max(channel[1] for channel in channels)

    if max_index > branch_ratios.shape[1]:
        raise ValueError(
            f"{species}: requires branching column {max_index}, "
            f"but only {branch_ratios.shape[1]} columns exist"
        )

    # Calculate the PACT yields
    yield_data = []

    for reaction_id, br_index, products in channels:

        # Branching indices in the supplied list start at 1
        br = np.interp(
            wavelength,
            branch_wavelength,
            branch_ratios[:, br_index - 1]
        )

        # Yield relative to total photoabsorption
        yield_values = np.zeros_like(absorption)

        valid = absorption > 0

        yield_values[valid] = (
            dissociation[valid]
            / absorption[valid]
        ) * br[valid]

        yield_data.append(yield_values)

    yield_data = np.array(yield_data).T

    # Write PACT file
    output_file = output_dir / f"{species}.phot"

    with open(output_file, "w") as f:

        f.write(f"{len(channels)}  {pact_label}\n")

        for reaction_id, br_index, products in channels:
            f.write(
                f" {species} + photon = {products}\n"
            )

        f.write("!\n")
        f.write("! Converted from VULCAN photochemistry data\n")
        f.write(f"! Cross sections from {cross_file.name}\n")
        f.write(f"! Branching ratios from {branch_file.name}\n")
        f.write("!\n")
        f.write("! sigma is the total photoabsorption cross section\n")
        f.write("! yields are relative to total photoabsorption\n")
        f.write("!\n")
        f.write("! wlth[nm]  sigma[cm2]   yield(s)...\n")
        f.write("!---------  -----------  -----------\n")

        for i in range(len(wavelength)):
            line = f"{wavelength[i]:10.3f}   {absorption[i]:.3E}"

            for j in range(len(channels)):
                line += f"   {yield_data[i, j]:.3E}"

            f.write(line + "\n")

    print(f"  Reactions: {len(channels)}")
    print(f"  Wavelength points: {len(wavelength)}")
    print(f"  Output: {output_file}")


print("\nVULCAN to PACT photodissociation conversion complete.")

phot_files = sorted(output_dir.glob("*.phot"))

with open(output_dir / "photoreac.ipho", "w") as f:
    f.write("!\n")
    f.write("! photo cross section files available\n")
    f.write("!\n")

    for phot_file in phot_files:
        f.write(phot_file.name + "\n")

print(f"Created photoreac.ipho with {len(phot_files)} files")