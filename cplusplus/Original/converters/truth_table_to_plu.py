#!/usr/bin/env python3
"""
truth_table_to_plu.py
======================

Converts a truth table (CSV) into the .plu benchmark format used by CGP++
(cgp-plusplus) for logic synthesis.

-----------------------------------------------------------------------
WHY THIS EXISTS
-----------------------------------------------------------------------
CGP++'s LogicSynthesisProblem does NOT read a truth table row-by-row.
Instead, it reads whole *columns* of the truth table, each packed bit-for-bit
into a single integer (this is what makes its fitness evaluation so fast --
one bitwise XOR + popcount evaluates every row at once). A row in a .plu
file therefore looks like:

    240 204 170   232 150

where 240, 204, 170 are the three INPUT columns (each a packed bit pattern
across every row of the truth table) and 232, 150 are the two OUTPUT
columns, for a problem with 3 inputs and 2 outputs.

For truth tables with 2^n rows > CHUNK_BITS (default 32, matching CGP++'s
hardcoded MAX_BITS), the table is split into multiple "instances" -- each
instance is one 32-row chunk, and the .plu ".p" field records how many
chunks there are.

-----------------------------------------------------------------------
BIT-ORDER CONVENTION (IMPORTANT)
-----------------------------------------------------------------------
Row r of the truth table (r = 0 .. 2^n - 1) is expected to correspond to
the standard binary-counting order of the inputs, with the FIRST input
column being the most-significant bit (slowest changing) and the LAST
input column being the least-significant bit (fastest changing):

    row 0: in0=0 in1=0 in2=0 ...
    row 1: in0=0 in1=0 in2=1 ...
    row 2: in0=0 in1=0 in2=... 1 0
    ...

This matches CGP++'s own example benchmarks (verified against add1c.plu:
column values 240/204/170 for a 3-input adder match exactly this
convention). If your CSV rows are NOT already in this canonical order,
this script will re-sort them for you automatically, as long as every
combination of inputs appears exactly once.

-----------------------------------------------------------------------
CSV INPUT FORMAT
-----------------------------------------------------------------------
- Plain CSV, comma-separated.
- Optional header row (auto-detected: if the first row contains any
  non-{0,1} token, it's treated as a header and column names are used
  for --input-cols/--output-cols matching).
- Every cell must be 0 or 1 (no don't-cares) -- see --dontcare-fill below
  if you need to resolve don't-cares before conversion.
- By default, the first N columns are treated as inputs and the
  remaining M columns as outputs, where N and M are given by
  --num-inputs and --num-outputs. Alternatively, pass --input-cols and
  --output-cols with explicit header names.

-----------------------------------------------------------------------
EXAMPLE
-----------------------------------------------------------------------
    python3 truth_table_to_plu.py my_table.csv my_circuit.plu \\
        --num-inputs 4 --num-outputs 1

    # or with explicit named columns:
    python3 truth_table_to_plu.py my_table.csv my_circuit.plu \\
        --input-cols A,B,C,D --output-cols F
"""

import argparse
import csv
import sys
from itertools import product


CHUNK_BITS_DEFAULT = 32  # matches CGP++'s hardcoded MAX_BITS in LogicSynthesisProblem.h


def error(msg):
    print(f"Error: {msg}", file=sys.stderr)
    sys.exit(1)


def looks_like_header(row):
    """A row is treated as a header if any cell isn't a plain 0/1 token."""
    for cell in row:
        cell = cell.strip()
        if cell not in ("0", "1"):
            return True
    return False


def read_csv_rows(path):
    with open(path, newline="") as f:
        reader = csv.reader(f)
        rows = [r for r in reader if r and any(c.strip() != "" for c in r)]
    if not rows:
        error("CSV file is empty.")
    return rows


def parse_table(path, num_inputs, num_outputs, input_cols, output_cols,
                 dontcare_fill):
    raw_rows = read_csv_rows(path)

    header = None
    if looks_like_header(raw_rows[0]):
        header = [h.strip() for h in raw_rows[0]]
        data_rows = raw_rows[1:]
    else:
        data_rows = raw_rows

    if not data_rows:
        error("No data rows found (only a header row was present).")

    total_cols = len(data_rows[0])

    # Determine which columns are inputs vs outputs
    if input_cols or output_cols:
        if header is None:
            error("--input-cols/--output-cols given, but no header row was "
                  "detected in the CSV to match them against.")
        try:
            in_idx = [header.index(c) for c in input_cols]
            out_idx = [header.index(c) for c in output_cols]
        except ValueError as e:
            error(f"Column name not found in CSV header: {e}")
        n = len(in_idx)
        m = len(out_idx)
    else:
        if num_inputs is None or num_outputs is None:
            error("Provide either --input-cols/--output-cols (with a header "
                  "row) or both --num-inputs and --num-outputs.")
        n, m = num_inputs, num_outputs
        if n + m != total_cols:
            error(f"--num-inputs ({n}) + --num-outputs ({m}) = {n+m}, but "
                  f"each data row has {total_cols} columns.")
        in_idx = list(range(n))
        out_idx = list(range(n, n + m))

    expected_rows = 1 << n  # 2^n
    table = {}  # tuple(input bits) -> tuple(output bits)

    for line_num, row in enumerate(data_rows, start=(2 if header else 1)):
        if len(row) != total_cols:
            error(f"Row {line_num}: expected {total_cols} columns, got {len(row)}.")

        def get_bit(raw):
            raw = raw.strip().lower()
            if raw in ("0", "1"):
                return int(raw)
            if raw in ("x", "-", "d", "dc", ""):
                if dontcare_fill is None:
                    error(
                        f"Row {line_num}: found a don't-care value ('{raw}') "
                        f"but no --dontcare-fill was given. Pass "
                        f"--dontcare-fill 0 or --dontcare-fill 1 to resolve "
                        f"don't-cares before conversion (the .plu format "
                        f"requires a fully specified table)."
                    )
                return dontcare_fill
            error(f"Row {line_num}: value '{raw}' is not 0, 1, or a "
                  f"recognized don't-care marker (x, -, d).")

        in_bits = tuple(get_bit(row[i]) for i in in_idx)
        out_bits = tuple(get_bit(row[i]) for i in out_idx)

        if in_bits in table:
            error(f"Row {line_num}: duplicate input combination {in_bits} "
                  f"(already seen earlier in the file).")
        table[in_bits] = out_bits

    if len(table) != expected_rows:
        missing = expected_rows - len(table)
        # Show a few missing combinations to help debugging
        all_combos = set(product([0, 1], repeat=n))
        missing_combos = list(all_combos - set(table.keys()))[:5]
        error(
            f"Truth table is incomplete: {n} inputs require {expected_rows} "
            f"rows (2^{n}), but only {len(table)} unique rows were found "
            f"({missing} missing). Example missing input combination(s): "
            f"{missing_combos}. The .plu format requires a fully specified "
            f"table -- use --dontcare-fill if these are don't-cares."
        )

    return n, m, table


def build_columns(n, m, table):
    """
    Returns:
      input_columns:  list of n integers, one full-width bit pattern per
                       input variable, ordered canonically by row index.
      output_columns: list of m integers, one full-width bit pattern per
                       output, ordered canonically by row index.
    Row index r corresponds to input combination given by the binary
    representation of r, with input 0 as the MOST significant bit.
    """
    num_rows = 1 << n
    input_columns = [0] * n
    output_columns = [0] * m

    for r in range(num_rows):
        # Reconstruct the input combination for row r (input 0 = MSB)
        in_bits = tuple((r >> (n - 1 - k)) & 1 for k in range(n))
        out_bits = table[in_bits]

        for k in range(n):
            if in_bits[k]:
                input_columns[k] |= (1 << r)
        for k in range(m):
            if out_bits[k]:
                output_columns[k] |= (1 << r)

    return input_columns, output_columns


def chunk_columns(n, input_columns, output_columns, chunk_bits):
    """
    Splits full-width columns into CGP++-style chunks of `chunk_bits` rows
    each (matching the ".p" instance mechanism in the .plu format). Returns
    a list of (input_chunk_values, output_chunk_values) tuples, one per
    instance.
    """
    num_rows = 1 << n
    num_chunks = max(1, (num_rows + chunk_bits - 1) // chunk_bits)
    mask = (1 << chunk_bits) - 1 if num_rows > chunk_bits else (1 << num_rows) - 1

    instances = []
    for c in range(num_chunks):
        shift = c * chunk_bits
        in_vals = [(col >> shift) & mask for col in input_columns]
        out_vals = [(col >> shift) & mask for col in output_columns]
        instances.append((in_vals, out_vals))
    return instances


def write_plu(path, n, m, instances):
    with open(path, "w", newline="\n") as f:
        f.write(f".i {n}\n")
        f.write(f".o {m}\n")
        f.write(f".p {len(instances)}\n")
        for in_vals, out_vals in instances:
            in_str = " ".join(str(v) for v in in_vals)
            out_str = " ".join(str(v) for v in out_vals)
            f.write(f"{in_str}   {out_str}\n")
        f.write(".e\n")


def main():
    parser = argparse.ArgumentParser(
        description="Convert a CSV truth table into a CGP++ .plu benchmark file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("csv_path", help="Path to the input CSV truth table.")
    parser.add_argument("plu_path", help="Path to write the output .plu file.")
    parser.add_argument("--num-inputs", type=int, default=None,
                         help="Number of input columns (if CSV has no usable header).")
    parser.add_argument("--num-outputs", type=int, default=None,
                         help="Number of output columns (if CSV has no usable header).")
    parser.add_argument("--input-cols", type=str, default=None,
                         help="Comma-separated input column names (requires a CSV header).")
    parser.add_argument("--output-cols", type=str, default=None,
                         help="Comma-separated output column names (requires a CSV header).")
    parser.add_argument("--dontcare-fill", type=int, choices=[0, 1], default=None,
                         help="Resolve don't-care cells (x, -, d) to this fixed value.")
    parser.add_argument("--chunk-bits", type=int, default=CHUNK_BITS_DEFAULT,
                         help=f"Rows per .plu instance chunk (default {CHUNK_BITS_DEFAULT}, "
                              f"matching CGP++'s hardcoded MAX_BITS -- only change this if "
                              f"you've also changed MAX_BITS in LogicSynthesisProblem.h).")
    args = parser.parse_args()

    input_cols = args.input_cols.split(",") if args.input_cols else None
    output_cols = args.output_cols.split(",") if args.output_cols else None

    n, m, table = parse_table(
        args.csv_path, args.num_inputs, args.num_outputs,
        input_cols, output_cols, args.dontcare_fill,
    )

    input_columns, output_columns = build_columns(n, m, table)
    instances = chunk_columns(n, input_columns, output_columns, args.chunk_bits)
    write_plu(args.plu_path, n, m, instances)

    num_rows = 1 << n
    print(f"OK: wrote {args.plu_path}")
    print(f"    inputs={n}  outputs={m}  rows=2^{n}={num_rows}  "
          f"instances(.p)={len(instances)}")


if __name__ == "__main__":
    main()
