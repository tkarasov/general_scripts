#!/usr/bin/env python3
# I frequently have to determine which ordered mutant is found in which well of our ordered mutant library. The following takes output files from Effie's pipeline and assigns the most likely mutant



# cd /uufs/chpc.utah.edu/common/home/karasov-group1/lab_members/esymeonidi/charact_sep26/final_comparison/dc3000

"""
Match barcodes to original mutant plate x well coordinates.

Expected inputs:
  dc3000_plate001_top_ids.txt ... dc3000_plate025_top_ids.txt
  dc3000_all_positions.txt

Output:
  barcode_plate_well_assignments.tsv
"""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

PLATE_FILE_RE = re.compile(r"plate(\d+)_top_ids\.txt$")
COUNT_RE = re.compile(r":(\d+)$")
# Examples matched: A01P02id01, H12P25id60, A01id01
POSITION_RE = re.compile(r"^([A-H]\d{2})(?:P\d+)?id\d+$", re.IGNORECASE)


def parse_line(line, source):
    """Return (pool_label, count, barcode), or None for blank lines."""
    fields = line.strip().split()
    if not fields:
        return None
    if len(fields) != 2:
        raise ValueError(f"Unexpected line in {source}: {line.rstrip()}")

    label, barcode = fields
    match = COUNT_RE.search(label)
    if not match:
        raise ValueError(f"Could not find ':count' in {source}: {line.rstrip()}")

    return label[:match.start()], int(match.group(1)), barcode.upper()


def best_call(scores):
    """
    scores: dict of coordinate -> aggregated read count
    returns best coordinate, best count, second count, ratio
    """
    if not scores:
        return "NA", 0, 0, "NA"

    ranked = sorted(scores.items(), key=lambda x: (-x[1], x[0]))
    best_name, best_count = ranked[0]
    second_count = ranked[1][1] if len(ranked) > 1 else 0
    ratio = "inf" if second_count == 0 else round(best_count / second_count, 3)

    return best_name, best_count, second_count, ratio


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--directory",
        default=".",
        help="Directory containing the top-ID files [default: current directory].",
    )
    parser.add_argument(
        "--positions-file",
        default="dc3000_all_positions.txt",
        help="Filename of the concatenated position-pool file.",
    )
    parser.add_argument(
        "--output",
        default="barcode_plate_well_assignments.tsv",
        help="Output TSV filename.",
    )
    parser.add_argument(
        "--min-plate-ratio",
        type=float,
        default=3.0,
        help="Minimum top/second plate count ratio for a confident plate call.",
    )
    parser.add_argument(
        "--min-well-ratio",
        type=float,
        default=3.0,
        help="Minimum top/second well count ratio for a confident well call.",
    )
    args = parser.parse_args()

    data_dir = Path(args.directory)
    position_file = data_dir / args.positions_file

    if not position_file.exists():
        raise FileNotFoundError(f"Cannot find {position_file}")

    plate_files = []
    for path in data_dir.glob("*plate*_top_ids.txt"):
        match = PLATE_FILE_RE.search(path.name)
        if match:
            plate_files.append((int(match.group(1)), path))

    plate_files.sort()

    if not plate_files:
        raise FileNotFoundError("No files matching *plate*_top_ids.txt were found.")

    print(f"Reading {len(plate_files)} plate-pool files.")
    print(f"Reading position pools from {position_file.name}")

    # barcode -> {plate: counts}
    plate_counts = defaultdict(lambda: defaultdict(int))

    # barcode -> {well: counts}
    well_counts = defaultdict(lambda: defaultdict(int))

    # Read 25 plate pools
    for plate_number, path in plate_files:
        plate_name = f"plate{plate_number:03d}"

        with path.open() as handle:
            for line in handle:
                parsed = parse_line(line, path.name)
                if parsed is None:
                    continue

                _, count, barcode = parsed
                plate_counts[barcode][plate_name] += count

    # Read 96 well/position pools.
    # A label such as A01P02id01 is assigned to well A01.
    unmatched_position_labels = set()

    with position_file.open() as handle:
        for line in handle:
            parsed = parse_line(line, position_file.name)
            if parsed is None:
                continue

            raw_pool, count, barcode = parsed
            match = POSITION_RE.match(raw_pool)

            if not match:
                unmatched_position_labels.add(raw_pool)
                continue

            well = match.group(1).upper()
            well_counts[barcode][well] += count

    if unmatched_position_labels:
        print(
            f"WARNING: {len(unmatched_position_labels)} unrecognized position labels "
            "were skipped. First few:",
            sorted(unmatched_position_labels)[:10],
        )

    all_barcodes = sorted(set(plate_counts) | set(well_counts))

    output_file = data_dir / args.output
    confident = 0
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

        for barcode in all_barcodes:
            plate, plate_n, plate_second_n, plate_ratio = best_call(
                plate_counts[barcode]
            )
            well, well_n, well_second_n, well_ratio = best_call(
                well_counts[barcode]
            )

            second_plate = "NA"
            if len(plate_counts[barcode]) > 1:
                second_plate = sorted(
                    plate_counts[barcode].items(),
                    key=lambda x: (-x[1], x[0])
                )[1][0]

            second_well = "NA"
            if len(well_counts[barcode]) > 1:
                second_well = sorted(
                    well_counts[barcode].items(),
                    key=lambda x: (-x[1], x[0])
                )[1][0]

            problems = []
            if plate == "NA":
                problems.append("missing_plate")
            elif plate_ratio != "inf" and plate_ratio < args.min_plate_ratio:
                problems.append("ambiguous_plate")

            if well == "NA":
                problems.append("missing_well")
            elif well_ratio != "inf" and well_ratio < args.min_well_ratio:
                problems.append("ambiguous_well")

            status = "confident" if not problems else ";".join(problems)

            if status == "confident":
                confident += 1
            else:
                ambiguous += 1

            writer.writerow([
                barcode,
                plate,
                plate_n,
                second_plate,
                plate_second_n,
                plate_ratio,
                well,
                well_n,
                second_well,
                well_second_n,
                well_ratio,
                status,
            ])

    print(f"\nWrote: {output_file}")
    print(f"Barcodes assessed: {len(all_barcodes)}")
    print(f"Confident assignments: {confident}")
    print(f"Ambiguous or incomplete assignments: {ambiguous}")


if __name__ == "__main__":
    main()
