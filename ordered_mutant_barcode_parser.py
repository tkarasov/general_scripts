#!/usr/bin/env python3
# I frequently have to determine which ordered mutant is found in which well of our ordered mutant library. The following takes output files from Effie's pipeline and assigns the most likely mutant

# cd /uufs/chpc.utah.edu/common/home/karasov-group1/lab_members/esymeonidi/charact_sep26/final_comparison/dc3000

"""Assign barcodes to plate x well coordinates from combinatorial pools.

The script reads 25 plate-pool ``*_top_ids.txt`` files plus a combined
position-pool file.  It error-corrects position-pool barcode sequences that
are exactly one nucleotide away from one uniquely identifiable, well-supported
barcode in the plate pools.
"""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


PLATE_FILE_RE = re.compile(r"plate(\d+)_top_ids\.txt$", re.IGNORECASE)
COUNT_RE = re.compile(r":(\d+)$")
# Examples: A01P02id01, H12P25id60, A01id01
POSITION_RE = re.compile(r"^([A-H]\d{2})(?:P\d+)?id\d+$", re.IGNORECASE)
BASES = "ACGT"


def parse_line(line, source):
    """Return (pool_name, read_count, barcode), or None for blank lines."""
    fields = line.strip().split()
    if not fields:
        return None
    if len(fields) != 2:
        raise ValueError(f"Unexpected line in {source}: {line.rstrip()}")

    label, barcode = fields
    count_match = COUNT_RE.search(label)
    if not count_match:
        raise ValueError(f"Could not parse ':count' in {source}: {line.rstrip()}")

    return label[:count_match.start()], int(count_match.group(1)), barcode.upper()


def make_one_error_lookup(valid_barcodes):
    """Map every one-substitution neighbor to all possible parent barcodes."""
    lookup = defaultdict(set)
    for barcode in valid_barcodes:
        for i, original_base in enumerate(barcode):
            for new_base in BASES:
                if new_base != original_base:
                    lookup[barcode[:i] + new_base + barcode[i + 1:]].add(barcode)
    return lookup


def rank_scores(scores):
    """Return highest- and second-highest-scoring labels and their counts."""
    ranked = sorted(scores.items(), key=lambda item: (-item[1], item[0]))
    if not ranked:
        return "NA", 0, "NA", 0, "NA"

    best_name, best_count = ranked[0]
    if len(ranked) == 1:
        return best_name, best_count, "NA", 0, "inf"

    second_name, second_count = ranked[1]
    ratio = round(best_count / second_count, 3) if second_count else "inf"
    return best_name, best_count, second_name, second_count, ratio


def ratio_is_confident(ratio, cutoff):
    return ratio == "inf" or ratio >= cutoff


def main():
    parser = argparse.ArgumentParser(
        description="Match barcodes to plate x well coordinates with one-base error correction."
    )
    parser.add_argument(
        "--directory", default=".", help="Directory containing the input files [default: current directory]."
    )
    parser.add_argument(
        "--positions-file", default="dc3000_all_positions.txt", help="Combined position-pool file."
    )
    parser.add_argument(
        "--output", default="barcode_plate_well_errorcorrected.tsv", help="Output TSV filename."
    )
    parser.add_argument(
        "--barcodes-per-plate", type=int, default=96,
        help="Number of expected mutant barcodes on each plate [default: 96]."
    )
    parser.add_argument(
        "--min-plate-ratio", type=float, default=1.25,
        help="Minimum top/second plate read ratio for a confident call [default: 1.25]."
    )
    parser.add_argument(
        "--min-well-ratio", type=float, default=1.25,
        help="Minimum top/second well read ratio for a confident call [default: 1.25]."
    )
    args = parser.parse_args()

    data_dir = Path(args.directory)
    position_file = data_dir / args.positions_file
    output_file = data_dir / args.output
    if not position_file.is_file():
        raise FileNotFoundError(f"Position-pool file not found: {position_file}")

    plate_files = []
    for path in data_dir.glob("*plate*_top_ids.txt"):
        match = PLATE_FILE_RE.search(path.name)
        if match:
            plate_files.append((int(match.group(1)), path))
    plate_files.sort()
    if not plate_files:
        raise FileNotFoundError("No files matching '*plate*_top_ids.txt' were found.")

    # barcode -> {plate: total read count}
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

    # Assign each observed sequence to its most enriched plate.  Recurrent
    # high-count contaminants are down-ranked; a barcode is never considered
    # as a candidate for a non-best plate.
    candidates_by_plate = defaultdict(list)
    for barcode, counts in plate_counts.items():
        ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
        best_plate, best_count = ranked[0]
        second_count = ranked[1][1] if len(ranked) > 1 else 0

        # Do not break exact ties between plate pools arbitrarily.
        if best_count == second_count:
            continue

        other_count = sum(counts.values()) - best_count
        # This gives high scores to abundant, plate-specific barcodes, while
        # penalizing sequences that recur broadly across the other 24 pools.
        specificity_score = (best_count * best_count) / (other_count + 1)
        candidates_by_plate[best_plate].append(
            (specificity_score, best_count, barcode)
        )

    # Choose 96 distinct, plate-enriched barcodes per plate.
    valid_barcodes = set()
    for plate_number, _ in plate_files:
        plate = f"plate{plate_number:03d}"
        candidates = sorted(
            candidates_by_plate[plate],
            key=lambda item: (-item[0], -item[1], item[2]),
        )
        selected = candidates[:args.barcodes_per_plate]
        valid_barcodes.update(barcode for _, _, barcode in selected)
        if len(selected) < args.barcodes_per_plate:
            print(
                f"WARNING: {plate} has only {len(selected)} "
                "plate-enriched candidate barcodes."
            )

    expected_total = len(plate_files) * args.barcodes_per_plate
    one_error_lookup = make_one_error_lookup(valid_barcodes)
    print(f"Expected barcode entries: {expected_total} "
          f"({len(plate_files)} plates x {args.barcodes_per_plate}).")
    print(f"Unique plate-enriched whitelist barcodes: {len(valid_barcodes)}.")
    if len(valid_barcodes) < expected_total:
        print("WARNING: fewer than 96 plate-enriched candidates were available for at least one plate.")

    # barcode -> {well: total read count}; well is A01 through H12.
    well_counts = defaultdict(lambda: defaultdict(int))
    exact_position_reads = corrected_position_reads = 0
    ambiguous_error_reads = unmatched_position_reads = 0
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

            # A label such as A01P02id01 is aggregated into well A01.
            well = position_match.group(1).upper()

            if observed_barcode in valid_barcodes:
                barcode = observed_barcode
                exact_position_reads += count
            else:
                possible_parents = one_error_lookup.get(observed_barcode, set())
                if len(possible_parents) == 1:
                    barcode = next(iter(possible_parents))
                    corrected_position_reads += count
                elif len(possible_parents) > 1:
                    # A one-base sequence that has multiple possible parents is
                    # deliberately not forced into an arbitrary barcode.
                    ambiguous_error_reads += count
                    continue
                else:
                    unmatched_position_reads += count
                    continue

            well_counts[barcode][well] += count

    if unknown_position_labels:
        examples = ", ".join(sorted(unknown_position_labels)[:10])
        print(f"WARNING: skipped {len(unknown_position_labels)} unknown position labels: {examples}")

    # Report the 25 x 96 plate-derived candidate barcodes, including any that
    # lack a position-pool call so missing wells remain visible.
    barcodes_to_report = sorted(valid_barcodes)
    confident = 0
    partial = 0
    ambiguous = 0

    with output_file.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "barcode", "assigned_plate", "plate_count", "second_plate",
            "second_plate_count", "plate_top_to_second_ratio", "assigned_well",
            "well_count", "second_well", "second_well_count",
            "well_top_to_second_ratio", "status"
        ])

        for barcode in barcodes_to_report:
            plate, plate_count, second_plate, second_plate_count, plate_ratio = rank_scores(plate_counts[barcode])
            well, well_count, second_well, second_well_count, well_ratio = rank_scores(well_counts[barcode])

            flags = []
            if plate == "NA":
                flags.append("missing_plate")
            elif not ratio_is_confident(plate_ratio, args.min_plate_ratio):
                flags.append("ambiguous_plate")

            if well == "NA":
                flags.append("missing_well")
            elif not ratio_is_confident(well_ratio, args.min_well_ratio):
                flags.append("ambiguous_well")

            status = "confident" if not flags else ";".join(flags)
            if status == "confident":
                confident += 1
            elif "missing_" in status:
                partial += 1
            else:
                ambiguous += 1

            writer.writerow([
                barcode, plate, plate_count, second_plate, second_plate_count,
                plate_ratio, well, well_count, second_well, second_well_count,
                well_ratio, status
            ])

    print("\nPosition-pool read handling:")
    print(f"  exact whitelist barcode reads: {exact_position_reads}")
    print(f"  corrected one-base variant reads: {corrected_position_reads}")
    print(f"  ambiguous one-base variant reads skipped: {ambiguous_error_reads}")
    print(f"  unmatched reads skipped: {unmatched_position_reads}")
    print(f"\nWrote: {output_file}")
    print(f"Whitelist barcodes assessed: {len(barcodes_to_report)}")
    print(f"Confident assignments: {confident}")
    print(f"Incomplete assignments: {partial}")
    print(f"Ambiguous assignments: {ambiguous}")


if __name__ == "__main__":
    main()
