from pathlib import Path
import re

SOURCE = Path("/home/jh1639/pact-main/data/chemistryvulcan/SNCHO_photo_network_2025.txt")
TEMPLATE = Path("/home/jh1639/pact-main/data/chemistryvulcan/network.reac")
OUTPUT = Path("/home/jh1639/pact-main/data/chemistryvulcan/SNCHO_vulcan_network.reac")

def convert(A, B, C):
    """A T^B exp(-C/T) -> alpha (T/300)^beta exp(-gamma/T)."""
    return A * 300.0**B, B, C


def fmt(x):
    return f"{x:.3E}"


NUM = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?'


def parse_coefficients(text, n):
    pattern = r'^\s*' + r'\s+'.join(f'({NUM})' for _ in range(n)) + r'(.*)$'
    m = re.match(pattern, text)
    if not m:
        return None
    values = [float(m.group(i + 1)) for i in range(n)]
    return values, m.group(n + 1).strip()


def extract_ranges(text):
    """
    Return two temperature ranges.

    Missing ranges default to 0-2500 K.
    If only one range is present, the second also defaults to 0-2500 K.
    """
    pairs = re.findall(
        r'(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)',
        text
    )

    if len(pairs) >= 2:
        #return pairs[0], pairs[1]
        return ("0", "2500"), ("0", "2500")
    
    if len(pairs) == 1:
        print(pairs[0])
        #return pairs[0], ("0", "2500")
        return ("0", "2500"), ("0", "2500")

    return ("0", "2500"), ("0", "2500")


def strip_ranges(text):
    return re.sub(
        r'(?<![\w.])\d+(?:\.\d+)?\s*-\s*\d+(?:\.\d+)?(?![\w.])',
        '',
        text
    ).strip(" \t,/")


def extract_lerror(text):
    """
    Try to find an integer/float lerror in the source tail.

    SCHOAr references often contain numbers, so only use a standalone
    numeric token at the end if it is clearly separated. Otherwise use 0.
    """
    # Explicitly labelled forms such as lerror=1 or error 1.
    m = re.search(
        r'(?:lerror|error)\s*[=:]\s*([+-]?\d+(?:\.\d+)?)\s*$',
        text,
        re.IGNORECASE
    )
    if m:
        return m.group(1), text[:m.start()].strip()

    # If no explicit lerror is present, use 0.
    return "A", text


# ----------------------------------------------------------------------
# Preserve the supplied network.reac header.
# ----------------------------------------------------------------------

template_lines = TEMPLATE.read_text(errors="replace").splitlines()

first_reaction = next(
    i for i, line in enumerate(template_lines)
    if re.match(r'^\s*[A-Za-z0-9_].*:', line)
    and not line.lstrip().startswith(("!", "#"))
)

output_lines = template_lines[:first_reaction] + [""]

source_lines = SOURCE.read_text(errors="replace").splitlines()

section = "two_body"
two_count = 0
pressure_count = 0
skipped = []


for line in source_lines:

    stripped = line.strip()

    if "3-body and Disscoiation Reactions" in stripped:
        section = "pressure"
        continue

    if not stripped or stripped.startswith("#"):
        continue

    # Accept both:
    # 1 [ reaction ] ...
    # and:
    # 1153 1153 [ reaction ] ...
    m = re.match(
        r'^\s*(\d+)\s+(?:(?:\d+)\s+)?\[\s*(.*?)\s*\]\s+(.*)$',
        line
    )

    if not m:
        continue

    rid, reaction, remainder = m.groups()
    reaction = reaction.strip().replace("->", "=")


    # ==================================================================
    # TYPE 1: ordinary reaction
    #
    # reaction : alpha beta gamma : Tmin Tmax : lerror : comment
    # ==================================================================

    if section == "two_body":

        parsed = parse_coefficients(remainder, 3)

        if parsed is None:
            skipped.append((rid, line))
            continue

        (A, B, C), tail = parsed

        alpha, beta, gamma = convert(A, B, C)

        (tmin, tmax), _ = extract_ranges(tail)

        # For ordinary reactions the source does not necessarily provide
        # an lerror value. Default to 0.
        lerror, reference = extract_lerror(strip_ranges(tail))

        output_lines.append(
            f"{reaction:<48} : "
            f"{fmt(alpha):>11} {beta:>9.3f} {gamma:>10.3f} : "
            f"{tmin:>7} {tmax:>7} : "
            f"{lerror:>3} : "
            f"{reference}"
        )

        two_count += 1


    # ==================================================================
    # PRESSURE-DEPENDENT REACTIONS
    #
    # reaction :
    # alpha beta gamma :
    # alphainf betainf gammainf :
    # afcent :
    # bfcent :
    # mcoll :
    # Tmin Tmax :
    # Tmininf Tmaxinf :
    # Tminfc Tmaxfc :
    # lerror :
    # comment
    #
    # A/B/M and falloff data are intentionally blank.
    # ==================================================================

    else:

        parsed6 = parse_coefficients(remainder, 6)

        if parsed6 is not None:

            (A0, B0, C0, Ainf, Binf, Cinf), tail = parsed6

            alpha0, beta0, gamma0 = convert(A0, B0, C0)
            alphainf, betainf, gammainf = convert(Ainf, Binf, Cinf)

            (range0, rangeinf) = extract_ranges(tail)

            lerror, reference = extract_lerror(
                strip_ranges(tail)
            )

            output_lines.append(
                f"{reaction:<48} : "
                f"{fmt(alpha0):>11} {beta0:>9.3f} {gamma0:>10.3f} : "
                f"{fmt(alphainf):>11} {betainf:>9.3f} {gammainf:>10.3f} : "
                f" : "
                f" : "
                f" : "
                f"{range0[0]:>7} {range0[1]:>7} : "
                f"{rangeinf[0]:>7} {rangeinf[1]:>7} : "
                f" : "
                f"{lerror:>3} : "
                f"{reference}"
            )

            pressure_count += 1
            continue


        # --------------------------------------------------------------
        # One-sided pressure-dependent / dissociation reaction.
        # --------------------------------------------------------------

        parsed3 = parse_coefficients(remainder, 3)

        if parsed3 is not None:

            (A0, B0, C0), tail = parsed3

            alpha0, beta0, gamma0 = convert(A0, B0, C0)

            (range0, _) = extract_ranges(tail)

            lerror, reference = extract_lerror(
                strip_ranges(tail)
            )

            output_lines.append(
                f"{reaction:<48} : "
                f"{fmt(alpha0):>11} {beta0:>9.3f} {gamma0:>10.3f} : "
                f"{'':>11} {'':>9} {'':>10} : "
                f" : "
                f" : "
                f" : "
                f"{range0[0]:>7} {range0[1]:>7} : "
                f"{'0':>7} {'2500':>7} : "
                f" : "
                f"{lerror:>3} : "
                f"{reference}"
            )

            pressure_count += 1
            continue


        # --------------------------------------------------------------
        # Include unparseable pressure-dependent reactions rather than
        # dropping them. All unavailable numerical fields are blank.
        # lerror is set to A so the Fortran reader can still read it.
        # --------------------------------------------------------------

        (range0, rangeinf) = extract_ranges(remainder)

        lerror, reference = extract_lerror(
            strip_ranges(remainder)
        )

        output_lines.append(
            f"{reaction:<48} : "
            f"{'':>11} {'':>9} {'':>10} : "
            f"{'':>11} {'':>9} {'':>10} : "
            f" : "
            f" : "
            f" : "
            f"{range0[0]:>7} {range0[1]:>7} : "
            f"{rangeinf[0]:>7} {rangeinf[1]:>7} : "
            f" : "
            f"{lerror:>3} : "
            f"{reference}"
        )

        pressure_count += 1


OUTPUT.write_text("\n".join(output_lines) + "\n")

print(f"Created: {OUTPUT}")
print(f"Two-body reactions converted: {two_count}")
print(f"Pressure-dependent/dissociation reactions included: {pressure_count}")
print(f"Total reactions written: {two_count + pressure_count}")
print(f"Skipped malformed source lines: {len(skipped)}")

for rid, line in skipped[:20]:
    print(f"SKIPPED {rid}: {line}")