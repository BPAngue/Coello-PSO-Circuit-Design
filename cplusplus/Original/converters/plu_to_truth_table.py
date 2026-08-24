#!/usr/bin/env python3

"""
plu_to_truth_table.py
=====================

Converts a CGP++ .plu benchmark file into a normal CSV truth table.

The .plu format stores each input/output column as a packed integer.
For example:

    .i 3
    .o 2
    .p 1
    240 204 170   232 150
    .e

represents a truth table with:

    2^3 = 8 rows

where:

    240 -> input 0
    204 -> input 1
    170 -> input 2
    232 -> output 0
    150 -> output 1

Each bit of each integer corresponds to one truth-table row.
"""

import argparse
import csv
import sys


CHUNK_BITS_DEFAULT = 32


def error(msg):
    print(f"Error: {msg}", file=sys.stderr)
    sys.exit(1)


def read_plu(path):
    """
    Read a .plu file.

    Returns:
        n              number of inputs
        m              number of outputs
        instances      list of (input_values, output_values)
    """

    n = None
    m = None
    p = None
    instances = []

    with open(path, "r") as f:
        for line_num, raw_line in enumerate(f, start=1):

            line = raw_line.strip()

            # Ignore empty lines
            if not line:
                continue

            # End of file marker
            if line == ".e":
                break

            # Number of inputs
            if line.startswith(".i"):
                parts = line.split()

                if len(parts) != 2:
                    error(f"Line {line_num}: invalid .i declaration")

                n = int(parts[1])

            # Number of outputs
            elif line.startswith(".o"):
                parts = line.split()

                if len(parts) != 2:
                    error(f"Line {line_num}: invalid .o declaration")

                m = int(parts[1])

            # Number of instances
            elif line.startswith(".p"):
                parts = line.split()

                if len(parts) != 2:
                    error(f"Line {line_num}: invalid .p declaration")

                p = int(parts[1])

            # Data line
            else:
                if n is None or m is None:
                    error(
                        f"Line {line_num}: encountered data before "
                        f".i/.o declarations"
                    )

                values = line.split()

                expected = n + m

                if len(values) != expected:
                    error(
                        f"Line {line_num}: expected {expected} values "
                        f"({n} inputs + {m} outputs), "
                        f"but found {len(values)}"
                    )

                try:
                    values = [int(v) for v in values]
                except ValueError:
                    error(f"Line {line_num}: non-integer value found")

                input_values = values[:n]
                output_values = values[n:]

                instances.append((input_values, output_values))

    if n is None:
        error("Missing .i declaration")

    if m is None:
        error("Missing .o declaration")

    if p is None:
        error("Missing .p declaration")

    if len(instances) != p:
        error(
            f".p says {p} instances, but {len(instances)} "
            f"data lines were found"
        )

    return n, m, instances


def unpack_plu(n, m, instances, chunk_bits=32):
    """
    Convert packed .plu integers into individual truth-table rows.

    Returns:
        rows = [
            [input0, input1, ..., output0, output1, ...],
            ...
        ]
    """

    # Total number of rows represented by the input variables
    total_rows = 1 << n

    rows = []

    for instance_index, (input_values, output_values) in enumerate(instances):

        # Number of rows represented by this instance
        start_row = instance_index * chunk_bits

        remaining_rows = total_rows - start_row

        if remaining_rows <= 0:
            break

        rows_in_instance = min(chunk_bits, remaining_rows)

        for bit_index in range(rows_in_instance):

            row = []

            # ---------------------------------------------------------
            # Inputs
            # ---------------------------------------------------------

            for value in input_values:

                bit = (value >> bit_index) & 1

                row.append(bit)

            # ---------------------------------------------------------
            # Outputs
            # ---------------------------------------------------------

            for value in output_values:

                bit = (value >> bit_index) & 1

                row.append(bit)

            rows.append(row)

    if len(rows) != total_rows:
        error(
            f"Expected {total_rows} truth-table rows, "
            f"but reconstructed {len(rows)}"
        )

    return rows


def write_csv(path, n, m, rows, include_header=True):
    """
    Write reconstructed truth table to CSV.
    """

    with open(path, "w", newline="") as f:

        writer = csv.writer(f)

        if include_header:

            header = []

            for i in range(n):
                header.append(f"in{i}")

            for i in range(m):
                header.append(f"out{i}")

            writer.writerow(header)

        writer.writerows(rows)


def print_truth_table(n, m, rows):
    """
    Print the truth table to the terminal.
    """

    headers = []

    for i in range(n):
        headers.append(f"I{i}")

    for i in range(m):
        headers.append(f"O{i}")

    # Header
    print(" | ".join(headers))

    print("-" * (len(" | ".join(headers))))

    # Rows
    for row in rows:
        print(" | ".join(str(x) for x in row))


def main():

    parser = argparse.ArgumentParser(
        description="Convert a CGP++ .plu file into a truth table."
    )

    parser.add_argument(
        "plu_path",
        help="Input .plu file"
    )

    parser.add_argument(
        "csv_path",
        help="Output CSV file"
    )

    parser.add_argument(
        "--chunk-bits",
        type=int,
        default=CHUNK_BITS_DEFAULT,
        help="Number of rows represented by each .plu instance "
             "(default: 32)"
    )

    parser.add_argument(
        "--no-header",
        action="store_true",
        help="Do not write CSV column headers"
    )

    parser.add_argument(
        "--print",
        action="store_true",
        help="Also print the reconstructed truth table"
    )

    args = parser.parse_args()

    # -------------------------------------------------------------
    # Read PLU
    # -------------------------------------------------------------

    n, m, instances = read_plu(args.plu_path)

    # -------------------------------------------------------------
    # Convert packed integers to individual bits
    # -------------------------------------------------------------

    rows = unpack_plu(
        n,
        m,
        instances,
        args.chunk_bits
    )

    # -------------------------------------------------------------
    # Write CSV
    # -------------------------------------------------------------

    write_csv(
        args.csv_path,
        n,
        m,
        rows,
        include_header=not args.no_header
    )

    # -------------------------------------------------------------
    # Information
    # -------------------------------------------------------------

    print(f"OK: converted {args.plu_path}")
    print(f"    inputs={n}")
    print(f"    outputs={m}")
    print(f"    rows={len(rows)}")
    print(f"    instances={len(instances)}")
    print(f"    output={args.csv_path}")

    if args.print:
        print()
        print_truth_table(n, m, rows)


if __name__ == "__main__":
    main()