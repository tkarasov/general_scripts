#!/usr/bin/env python3
# I frequently have to determine which ordered mutant is found in which well of our ordered mutant library. The following takes output files from Effie's pipeline and assigns the most likely mutant

# cd /uufs/chpc.utah.edu/common/home/karasov-group1/lab_members/esymeonidi/charact_sep26/final_comparison/dc3000

"""
Assign each barcode to its original plate x well coordinate.

Input:
  dc3000_plate001_top_ids.txt ... dc3000_plate025_top_ids.txt
  dc3000_all_positions.txt

Logic:
  - Take the top 96 barcodes from each of the 25 plate pools:
      25 plates x 96 wells = 2,400 expected mutants.
  - Use these as the true-barcode whitelist.
  - Aggregate position-pool data by well (A01-H12).
  - Correct unique one-base barcode variants to the whitelist barcode.
  - Call the best plate and well for each barcode.
"""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

PLATE_FILE_RE = re.compile(r"plate(\d+)_top_ids\.txt$", re.IGNORECASE)
COUNT_RE = re.compile(r":(\d+)$")
POSITION_RE = re.compile(r"^([A-H]\d{2})(?:P\d+)?id\d+$", re.IGNORECASE)
BASES = "ACGT"


def parse_line(line, source):
    """Return (pool_name, count, barcode), or None for blank lines."""
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
    Make a lookup table of all one-base substitutions.

    Example:
      true barcode = ACGT
      possible observed one-base variants = CCGT, GCGT, TCGT, AAGT, etc.

    The value is a set because a variant could theoretically be one base
    away from more than one true barcode.
    """
    lookup = defaultdict(set)

    for barcode in valid_barcodes:
        for i, original_base in enumerate(barcode):
            for new_base in BASES:
                if new_base != original_base:
                    variant = (
                        barcode[:i] + new_base + barcode[i + 1:]
                    )
                    lookup[variant].add(barcode)

    return lookup


def rank_scores(scores):
    """
    Given {coordinate: read_count}, return:
      best_coordinate, best_count,
      second_coordinate, second_count,
      top/second ratio
    """
    ranked = sorted(scores.items(), key=lambda x: (-x[1], x[0]))

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


def is_confident(ratio, cutoff):
    return ratio == "inf" or ratio >= cutoff


def main():
    parser = argparse.ArgumentParser()

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
        default="barcode_plate_well_assignments.tsv",
        help="Output TSV filename.",
    )

    parser.add_argument(
        "--barcodes-per-plate",
        type=int,
        default=96,
        help="Number of expected true barcodes per plate [default: 96].",
    )

    parser.add_argument(
        "--min-plate-ratio",
        type=float,
        default=1.25,
        help="Top/second plate-count ratio for confident calls [default: 1.25].",
    )

    parser.add_argument(
        "--min-well-ratio",
        type=float,
        default=1.25,
        help="Top/second well-count ratio for confident calls [default: 1.25].",
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

    print(f"Reading {len(plate_files)} plate-pool files.")

    # barcode -> {plate: total reads}
    plate_counts = defaultdict(lambda: defaultdict(int))

    # plate -> the top 96 expected barcode sequences for that plate
    top_barcodes_by_plate = {}

    for plate_number, path in plate_files:
        plate = f"plate{plate_number:03d}"
        top_barcodes = []

        with path.open() as handle:
            for line in handle:
                parsed = parse_line(line, path.name)

                if parsed is None:
                    continue

                _, count, barcode = parsed
                plate_counts[barcode][plate] += count

                # Files are already in descending rank order: id001, id002, etc.
                if len(top_barcodes) < args.barcodes_per_plate:
                    top_barcodes.append(barcode)

        top_barcodes_by_plate[plate] = top_barcodes

        if len(top_barcodes) < args.barcodes_per_plate:
            print(
                f"WARNING: {plate} contained only {len(top_barcodes)} "
                f"barcodes, not {args.barcodes_per_plate}."
            )

    # The true mutant-barcode whitelist:
    # top 96 calls from each of 25 plates = expected 2,400 mutants.
    valid_barcodes = set()

    for plate, barcodes in top_barcodes_by_plate.items():
        valid_barcodes.update(barcodes)

    expected_total = len(plate_files) * args.barcodes_per_plate

    print(
        f"Expected barcode entries: {expected_total} "
        f"({len(plate_files)} plates x {args.barcodes_per_plate})."
    )
    print(f"Unique whitelist barcodes: {len(valid_barcodes)}.")

    if len(valid_barcodes) < expected_total:
        print(
            "WARNING: fewer unique whitelist barcodes than expected. "
            "This means at least one barcode appeared among the top calls "
            "for more than one plate."
        )

    one_error_lookup = make_one_error_lookup(valid_barcodes)

    # barcode -> {well: total reads}
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

            # A label such as A01P02id01 is assigned to well A01.
            well = position_match.group(1).upper()

            # Exact match to a true barcode.
            if observed_barcode in valid_barcodes:
                barcode = observed_barcode
                exact_position_reads += count

            else:
                possible_parents = one_error_lookup.get(
                    observed_barcode,
                    set()
                )

                # Correct only if this variant has one unambiguous parent.
                if len(possible_parents) == 1:
                    barcode = next(iter(possible_parents))
                    corrected_position_reads += count

                elif len(possible_parents) > 1:
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

        for barcode in sorted(valid_barcodes):
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

            elif not is_confident(plate_ratio, args.min_plate_ratio):
                flags.append("ambiguous_plate")

            if well == "NA":
                flags.append("missing_well")

            elif not is_confident(well_ratio, args.min_well_ratio):
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
    print(f"  corrected one-base variant reads: {corrected_position_reads}")
    print(
        f"  ambiguous one-base variants skipped: {ambiguous_error_reads}"
    )
    print(f"  unmatched reads skipped: {unmatched_position_reads}")

    print(f"\nWrote: {output_file}")
    print(f"Whitelist barcodes assessed: {len(valid_barcodes)}")
    print(f"Confident assignments: {confident}")
    print(f"Incomplete assignments: {incomplete}")
    print(f"Ambiguous assignments: {ambiguous}")


if __name__ == "__main__":
    main()
