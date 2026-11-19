import re
from collections import OrderedDict
import sys




#!/usr/bin/env python3

import re
from pathlib import Path


# ============================================================
# FILES
# ============================================================

NASA_FILE = "/home/jh1639/pact-main/data/chemistryvulcan/nasa9_combined.txt"
OUTPUT_FILE = "/home/jh1639/pact-main/data/chemistryvulcan/therm_vulcan8.therm"
NETWORK_FILE = "/home/jh1639/pact-main/data/chemistryvulcan/SNCHO_photo_network_2025.txt"


NUMBER_RE = re.compile(
    r"""
    [+-]?
    (?:
        \d+(?:\.\d*)?
        |
        \.\d+
    )
    (?:[EeDd][+-]?\d+)?
    """,
    re.VERBOSE,
)


def extract_numbers(line):
    """Extract floating point numbers from a line."""

    return [
        float(x.replace("D", "E").replace("d", "e"))
        for x in NUMBER_RE.findall(line)
    ]


# ============================================================
# SPECIES NAME -> CHEMICAL FORMULA
# ============================================================

def species_to_formula(species):
    """
    Convert a NASA species name into the chemical formula
    used to determine elemental abundances.

    The original species name is NOT changed.
    """

    formula = species

    # Remove phase/state suffixes only for the formula.
    formula = re.sub(
        r"_(?:p|m|s|l|g|l_s)$",
        "",
        formula
    )

    # Remove numerical suffixes such as _1 or _2.
    formula = re.sub(
        r"_\d+$",
        "",
        formula
    )

    # Remove labels such as _2D.
    formula = re.sub(
        r"_\d+[A-Za-z]+$",
        "",
        formula
    )

    return formula


# ============================================================
# CHEMICAL FORMULA PARSER
# ============================================================

ELEMENT_RE = re.compile(
    r"([A-Z][a-z]?)(\d*)"
)


def parse_formula(formula):
    """
    Parse a chemical formula and combine repeated elements.

    Examples:

        H       -> H 1
        H2O     -> H 2, O 1
        CO2     -> C 1, O 2
        C2H3CN  -> C 3, H 3, N 1
        CH3NH2  -> C 1, H 5, N 1
    """

    element_counts = {}

    position = 0

    while position < len(formula):

        match = ELEMENT_RE.match(
            formula,
            position
        )

        if match is None:

            raise ValueError(
                f"Cannot parse chemical formula "
                f"'{formula}' near "
                f"'{formula[position:]}'"
            )

        element = match.group(1)

        number = match.group(2)

        if number == "":
            number = 1
        else:
            number = int(number)

        # Combine repeated elements
        if element in element_counts:
            element_counts[element] += number
        else:
            element_counts[element] = number

        position = match.end()

    return list(element_counts.items())

# ============================================================
# NASA-9 COEFFICIENT CONVERSION
# ============================================================

def convert_nasa9(values, species):
    """
    NASA-9 format:

        a1 a2 a3 a4 a5 a6 a7 0 a8 a9

    ATMO DataType=9 expects:

        a1 a2 a3 a4 a5 a6 a7 a8 a9

    Therefore remove the zero from every 10-value block.
    """

    if len(values) == 0:

        raise ValueError(
            f"{species}: no thermochemical coefficients found"
        )

    if len(values) % 10 != 0:

        raise ValueError(
            f"{species}: found {len(values)} coefficient "
            f"values; expected a multiple of 10."
        )

    number_of_ranges = len(values) // 10

    coefficients = []

    for i in range(number_of_ranges):

        block = values[
            i * 10:
            (i + 1) * 10
        ]

        # NASA-9 eighth value should be zero.
        placeholder = block[7]

        if abs(placeholder) > 1.0e-10:

            raise ValueError(
                f"{species}: NASA-9 placeholder is "
                f"{placeholder}, not zero"
            )

        # Remove the zero.
        coefficients.extend(
            block[:7] + block[8:]
        )

        if number_of_ranges == 3:

            blocks = [
                coefficients[0:9],
                coefficients[9:18],
                coefficients[18:27]
            ]

            coefficients = (
                blocks[2] +
                blocks[0] +
                blocks[1]
            )

    return coefficients, number_of_ranges


# ============================================================
# TEMPERATURE RANGES
# ============================================================

def get_temperatures(number_of_ranges):

    if number_of_ranges == 2:

        return [
            200.0,
            1000.0,
            6000.0
        ]

    if number_of_ranges == 3:

        return [
            50.0,
            200.0,
            1000.0,
            6000.0
        ]

    raise ValueError(
        f"No temperature definition implemented for "
        f"{number_of_ranges} ranges"
    )


# ============================================================
# READ NASA-9 FILE
# ============================================================

def read_nasa9_file(filename):
    """
    Read nasa9_combined.txt.

    A species starts with:

        species.txt

    and continues until the next species.txt.

    Only the FIRST thermochemical dataset for each species
    is used.
    """

    with open(
        filename,
        "r",
        encoding="utf-8"
    ) as f:

        lines = f.readlines()

    species_data = []

    i = 0

    while i < len(lines):

        line = lines[i].strip()

        # ----------------------------------------------------
        # Find species name
        # ----------------------------------------------------

        if not line.endswith(".txt"):

            i += 1
            continue

        species = line[:-4]

        i += 1

        coefficients = []

        comment = ""

        # ----------------------------------------------------
        # Read first dataset
        # ----------------------------------------------------

        while i < len(lines):

            current = lines[i].strip()

            # Next species
            if current.endswith(".txt"):

                break

            # Comments
            if current.startswith("#"):

                if coefficients:

                    comment = current[1:].strip()

                    i += 1

                    break

                i += 1
                continue

            # Numerical coefficient line
            if current:

                values = extract_numbers(
                    current
                )

                if values:

                    coefficients.extend(
                        values
                    )

            i += 1

        species_data.append(
            {
                "name": species,
                "coefficients": coefficients,
                "comment": comment,
            }
        )

    return species_data


# ============================================================
# FORMAT SPECIES FOR ATMO
# ============================================================

def format_species(data):
    """
    Produce exactly the sequence expected by the Fortran reader.

    Line 1:
        species nattot nelem charge type

    Line 2:
        element number element number ...

    Line 3:
        ntemp temperatures...

    Following lines:
        nine coefficients per temperature range
    """

    species = data["name"]

    raw_coefficients = data["coefficients"]

    comment = data["comment"]

    # --------------------------------------------------------
    # Chemical composition
    # --------------------------------------------------------

    formula = species_to_formula(
        species
    )

    elements = parse_formula(
        formula
    )

    nattot = sum(
        number
        for element, number in elements
    )

    nelem = len(elements)

    # --------------------------------------------------------
    # NASA-9 coefficients
    # --------------------------------------------------------

    coefficients, ntemp = convert_nasa9(
        raw_coefficients,
        species
    )

    temperatures = get_temperatures(
        ntemp
    )

    # --------------------------------------------------------
    # LINE 1
    # --------------------------------------------------------

    line1 = (
        f"{species} "
        f"{nattot} "
        f"{nelem} "
        f"0 "
        f"9"
    )

    if comment:

        line1 += (
            f" ! {comment}"
        )

    # --------------------------------------------------------
    # LINE 2
    # --------------------------------------------------------

    element_parts = []

    for element, number in elements:

        element_parts.append(
            element
        )

        element_parts.append(
            str(number)
        )

    line2 = " ".join(
        element_parts
    )

    # --------------------------------------------------------
    # LINE 3
    # --------------------------------------------------------

    line3 = str(ntemp)

    for temperature in temperatures:

        line3 += (
            f" {temperature:.3f}"
        )

    # --------------------------------------------------------
    # COEFFICIENT LINES
    # --------------------------------------------------------

    coefficient_lines = []

    for k in range(ntemp):

        start = k * 9
        end = start + 9

        block = coefficients[
            start:end
        ]

        if len(block) != 9:

            raise ValueError(
                f"{species}: expected 9 coefficients "
                f"for temperature range {k + 1}"
            )

        coefficient_lines.append(
            " ".join(
                f"{value:.9E}"
                for value in block
            )
        )

    return [
        line1,
        line2,
        line3,
        *coefficient_lines
    ]


# ============================================================
# VALIDATE SPECIES
# ============================================================

def validate_species(data):

    species = data["name"]

    formula = species_to_formula(
        species
    )

    elements = parse_formula(
        formula
    )

    if len(elements) == 0:

        raise ValueError(
            f"{species}: no elements found"
        )

    coefficients, ntemp = convert_nasa9(
        data["coefficients"],
        species
    )

    temperatures = get_temperatures(
        ntemp
    )

    if len(temperatures) != ntemp + 1:

        raise ValueError(
            f"{species}: temperature count mismatch"
        )

    if len(coefficients) != 9 * ntemp:

        raise ValueError(
            f"{species}: expected "
            f"{9 * ntemp} coefficients, got "
            f"{len(coefficients)}"
        )


# ============================================================
# MAIN
# ============================================================

def main():

    print()
    print(
        "Reading NASA-9 data:"
    )
    print(
        f"  {NASA_FILE}"
    )

    data = read_nasa9_file(
        NASA_FILE
    )

    print(
        f"Found {len(data)} NASA species."
    )

    # --------------------------------------------------------
    # Validate
    # --------------------------------------------------------

    valid_data = []

    failed = []

    for entry in data:

        try:

            validate_species(
                entry
            )

            valid_data.append(
                entry
            )

        except Exception as error:

            failed.append(
                (
                    entry["name"],
                    str(error)
                )
            )

    # --------------------------------------------------------
    # Check duplicate species
    # --------------------------------------------------------

    names = [
        entry["name"]
        for entry in valid_data
    ]

    duplicates = set(
        name
        for name in names
        if names.count(name) > 1
    )

    if duplicates:

        print()
        print(
            "ERROR: duplicate species:"
        )

        for name in sorted(duplicates):

            print(
                f"  {name}"
            )

        sys.exit(1)

    # --------------------------------------------------------
    # Write output
    # --------------------------------------------------------

    with open(
        OUTPUT_FILE,
        "w",
        encoding="utf-8"
    ) as f:

        # Header
        f.write(
            "!\n"
            "! NASA-9 thermochemical data\n"
            "! Converted for ATMO DataType = 9\n"
            "!\n"
        )

        for entry in valid_data:

            lines = format_species(
                entry
            )

            for line in lines:

                f.write(
                    line + "\n"
                )

            # NO BLANK LINE HERE

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    print()
    print(
        "Conversion complete."
    )

    print(
        f"  Output:   {OUTPUT_FILE}"
    )

    print(
        f"  Written:  {len(valid_data)} species"
    )

    print(
        f"  Failed:   {len(failed)} species"
    )

    if failed:

        print()
        print(
            "Species that could not be converted:"
        )

        for name, error in failed:

            print(
                f"  {name}: {error}"
            )

    print()


if __name__ == "__main__":

    main()