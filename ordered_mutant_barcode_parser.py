#!/usr/bin/env python3
# I frequently have to determine which ordered mutant is found in which well of our ordered mutant library. The following takes output files from Effie's pipeline and assigns the most likely mutant

# cd /uufs/chpc.utah.edu/common/home/karasov-group1/lab_members/esymeonidi/charact_sep26/final_comparison/dc3000
"""
Assign barcodes to plate x well coordinates from combinatorial pools.

Reads:
  dc3000_plate001_top_ids.txt ... dc3000_plate025_top_ids.txt
  dc3000_all_positions.txt

The script:
  1. Calls a barcode's plate from the 25 plate pools.
  2. Calls its well from the 96 position pools.
  3. Corrects position-pool barcodes that are one base away from a
     uniquely matching, well-supported barcode in the plate pools.
  4. Reports confident, ambiguous, and incomplete calls.
"""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

PLATE_FILE_RE = re.compile(r"plate(\d+)_top_ids\.txt$", re.IGNORECASE)
COUNT_RE = re.compile(r":(\d+)$")
# Matches A01P02id01, H12P25id60, or A01id01
POSITION_RE = re.compile(r"^([A-H]\d{2})(?:P\d+)?id\d+$", re.IGNORECASE)
BASES = "ACGT"


def parse_line(line, source):
    """Return (pool_name, count, barcode); return None for blank lines."""
    fields = line.strip().split()

    if not fields:
        return None

    if len(fields) != 2:
        raise ValueError(f"Unexpected line in {source}: {line.rstrip()}")

    label, barcode = fields

    count_match = COUNT_RE.search(label)
    if not count_match:
        raise ValueError(
            f"Could not parse ':count' from {source}: {line.rstrip()}"
        )

    pool_name = label[:count_match.start()]
    count = int(count_match.group(1))

    return pool_name, count, barcode.upper()


def make_one_error_lookup(valid_barcodes):
    """
    For each valid barcode, generate all one-substitution neighbors.

    Returns:
        dict:
          erroneous sequence -> set of possible true parent barcode(s)
    """
    lookup = defaultdict(set)

    for barcode in valid_barcodes:
        for i, original_base in enumerate(barcode):
            for new_base in BASES:
                if new_base != original_base:
                    neighbor = (
                        barcode[:i] + new_base + barcode[i + 1:]
                    )
                    lookup[neighbor].add(barcode)

    return lookup


def rank_scores(scores):
    """
    Given {coordinate: read_count}, return:
      best_coordinate, best_count,
      second_coordinate, second_count,
      best/second ratio
    """
    ranked = sorted(scores.items(), key=lambda item: (-item[1], item[0]))

    if not ranked:
        return "NA", 0, "NA", 0, "NA"

    best_name, best_count = ranked[0]

    if len(ranked) == 1:
        return best_name, best_count, "NA", 0, "inf"

    second_name, second_count = ranked[1]

    if second_count == 0:
        ratio = "inf"
    else:
        ratio = round(best_count / second_count, 3)

    return best_name, best_count, second_name, second_count, ratio


def ratio_is_confident(ratio, cutoff):
    """Return True if the top call exceeds the chosen confidence cutoff."""
    return ratio == "inf" or ratio >= cutoff


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Match barcodes to plate x well coordinates, with "
            "one-base barcode error correction."
        )
    )

    parser.add_argument(
        "--directory",
        default=".",
        help="Directory containing input files [default: current directory].",
    )

    parser.add_argument(
        "--positions-file",
        default="dc3000_all_positions.txt",
        help="Combined position-pool file.",
    )

    parser.add_argument(
        "--output",
        default="barcode_plate_well_errorcorrected.tsv",
        help="Output TSV filename.",
    )

    parser.add_argument(
        "--min-whitelist-reads",
        type=int,
        default=250,
        help=(
            "Minimum summed plate-pool reads for a barcode to be used "
            "as a true barcode for one-base error correction [default: 250]."
        ),
    )

    parser.add_argument(
        "--min-plate-ratio",
        type=float,
        default=1.25,
        help=(
            "Minimum top/second plate-count ratio for a confident plate call "
            "[default: 1.25]."
        ),
    )

    parser.add_argument(
        "--min-well-ratio",
        type=float,
        default=1.25,
        help=(
            "Minimum top/second well-count ratio for a confident well call "
            "[default: 1.25]."
        ),
    )

    args = parser.parse_args()

    data_dir = Path(args.directory)
    position_file = data_dir / args.positions_file
    output_file = data_dir / args.output

    if not position_file.is_file():
        raise FileNotFoundError(
            f"Position-pool file not found: {position_file}"
        )

    plate_files = []

    for path in data_dir.glob("*plate*_top_ids.txt"):
        match = PLATE_FILE_RE.search(path.name)

        if match:
            plate_number = int(match.group(1))
            plate_files.append((plate_number, path))

    plate_files.sort()

    if not plate_files:
        raise FileNotFoundError(
            "No files matching '*plate*_top_ids.txt' were found."
        )

    # barcode -> {plate: total count}
    plate_counts = defaultdict(lambda: defaultdict(int))

    print(f"Reading {len(plate_files)} plate-pool files.")

    for plate_number, path in plate_files:
        plate = f"plate{plate_number:03d}"

        with path.open() as handle:
            for line in handle:
                parsed = parse_line(line, path.name)

                if parsed is None:
                    continue

                _, count, barcode = parsed
                plate_counts[barcode][plate] += count

    # Only well-supported plate-pool barcodes can be parents for correcting
    # sequencing/PCR variants in the position pools.
    valid_barcodes = {
        barcode
        for barcode, counts in plate_counts.items()
        if sum(counts.values()) >= args.min_whitelist_reads
    }

    one_error_lookup = make_one_error_lookup(valid_barcodes)

    print(
        f"Plate-barcode whitelist: {len(valid_barcodes)} barcodes "
        f"(minimum {args.min_whitelist_reads} reads)."
    )

    # barcode -> {well: total count}
    well_counts = defaultdict(lambda: defaultdict(int))

    exact_position_reads = 0
    corrected_position_reads = 0
    ambiguous_error_reads = 0
    unmatched_position_reads = 0
    unknown_position_labels = set()

    print(f"Reading position pools from {position_file.name}")

    with position_file.open() as handle:
        for line in handle:
            parsed = parse_line(line, position_file.name)

            if parsed is None:
                continue

            raw_pool, count, observed_barcode = parsed

            position_match = POSITION_RE.match(raw_pool)

            if not position_match:
                unknown_position_labels.add(raw_pool)
                continue

            # For example, A01P02id01 is assigned to well A01.
            well = position_match.group(1).upper()

            # Exact barcode call.
            if observed_barcode in valid_barcodes:
                barcode = observed_barcode
                exact_position_reads += count

            # Correct a unique one-base sequence variant.
            else:
                possible_parents = one_error_lookup.get(
                    observed_barcode,
                    set()
                )

                if len(possible_parents) == 1:
                    barcode = next(iter(possible_parents))
                    corrected_position_reads += count

                elif len(possible_parents) > 1:
                    # Do not arbitrarily choose between two possible parents.
                    ambiguous_error_reads += count
                    continue

                else:
                    unmatched_position_reads += count
                    continue

            well_counts[barcode][well] += count

    if unknown_position_labels:
        examples = ", ".join(sorted(unknown_position_labels)[:10])

        print(
            f"WARNING: skipped {len(unknown_position_labels)} "
            f"unrecognized position labels. Examples: {examples}"
        )

    # Report only plate-supported barcode identities, rather than treating
    # low-count sequence errors as independent mutants.
    barcodes_to_report = sorted(valid_barcodes)

    confident = 0
    incomplete = 0
    ambiguous = 0

    with output_file.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")

        writer.writerow([
            "barcode",
            "assigned_plate",
            "plate_count",
            "second_plate",
            "second_plate_count",
            "plate_top_to_second_ratio",
            "assigned_well",
            "well_count",
            "second_well",
            "second_well_count",
            "well_top_to_second_ratio",
            "status",
        ])

        for barcode in barcodes_to_report:
            (
                plate,
                plate_count,
                second_plate,
                second_plate_count,
                plate_ratio,
            ) = rank_scores(plate_counts[barcode])

            (
                well,
                well_count,
                second_well,
                second_well_count,
                well_ratio,
            ) = rank_scores(well_counts[barcode])

            flags = []

            if plate == "NA":
                flags.append("missing_plate")

            elif not ratio_is_confident(
                plate_ratio,
                args.min_plate_ratio,
            ):
                flags.append("ambiguous_plate")

            if well == "NA":
                flags.append("missing_well")

            elif not ratio_is_confident(
                well_ratio,
                args.min_well_ratio,
            ):
                flags.append("ambiguous_well")

            status = "confident" if not flags else ";".join(flags)

            if status == "confident":
                confident += 1

            elif "missing_" in status:
                incomplete += 1

            else:
                ambiguous += 1

            writer.writerow([
                barcode,
                plate,
                plate_count,
                second_plate,
                second_plate_count,
                plate_ratio,
                well,
                well_count,
                second_well,
                second_well_count,
                well_ratio,
                status,
            ])

    print("\nPosition-pool read handling:")
    print(f"  exact whitelist barcode reads: {exact_position_reads}")
    print(
        f"  corrected one-base variant reads: "
        f"{corrected_position_reads}"
    )
    print(
        f"  ambiguous one-base variant reads skipped: "
        f"{ambiguous_error_reads}"
    )
    print(f"  unmatched reads skipped: {unmatched_position_reads}")

    print(f"\nWrote: {output_file}")
    print(f"Whitelist barcodes assessed: {len(barcodes_to_report)}")
    print(f"Confident assignments: {confident}")
    print(f"Incomplete assignments: {incomplete}")
    print(f"Ambiguous assignments: {ambiguous}")


if __name__ == "__main__":
    main()
