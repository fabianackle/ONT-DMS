#!/usr/bin/env python3
"""
Reads TSV from stdin from bcftools query, expects CHROM, QUAL, FILTER, and TBCSQ columns.
Outputs CSV to stdout; CSV with header and the following columns:
barcode_id (CHROM), quality, filter, consequence, aa_change, dna_change.
"""

import csv
import sys

header = [
    "barcode_id",
    "quality",
    "filter",
    "consequence",
    "aa_change",
    "dna_change"
]


def main():
    reader = csv.reader(sys.stdin, delimiter='\t')
    writer = csv.writer(sys.stdout, delimiter=',')
    writer.writerow(header)

    for row in reader:
        if len(row) != 4:
            sys.stderr.write(f"Skipping row; not 4 columns: {row}\n")
            continue

        barcode_id, quality, filt, tbsqc = row
        parts = tbsqc.split('|')

        if len(parts) == 7:
            consequence, _, _, _, _, aa_change, dna_change = parts

            writer.writerow([
                barcode_id,
                quality,
                filt,
                consequence,
                aa_change,
                dna_change
            ])
        else:
            sys.stderr.write(f"Skipping row; tbsqc has {len(parts)} parts (expected 7): {tbsqc}\n")


if __name__ == '__main__':
    main()
