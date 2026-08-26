#!/usr/bin/env python3
"""
Truth-table simulator for the supplied .i/.o/.p structure.

Inputs:
    A, B, C, D, E, F, G

Outputs:
    Y1 = AND1(F, A)
    Y2 = XOR1(D, AND1(OR1(A, B), F))
    Y3 = AND1(OR1(A, B), XOR1(E, G))

The .p lines in the supplied file are interpreted as packed bit-vectors:
each integer is a 16-bit block, with bit 0 corresponding to the first
row of that block. For this example, the five pattern rows represent
5 combinations of A/B/C; the remaining D-G inputs are fixed by the
bit-pattern columns.

Run:
    python truth_table_simulator.py
"""

from dataclasses import dataclass


@dataclass
class Row:
    A: int
    B: int
    C: int
    D: int
    E: int
    F: int
    G: int
    expected: tuple[int, int, int]


def AND1(x: int, y: int) -> int:
    return x & y


def OR1(x: int, y: int) -> int:
    return x | y


def XOR1(x: int, y: int) -> int:
    return x ^ y


def circuit(A, B, C, D, E, F, G):
    # Structure from the question.
    and1 = AND1(F, A)
    or1 = OR1(A, B)
    and2 = AND1(or1, F)
    xor1 = XOR1(D, and2)
    xor2 = XOR1(E, G)
    y1 = and1
    y2 = xor1
    y3 = AND1(or1, xor2)
    return y1, y2, y3


def bits16(value: int):
    """Return 16 bits, least-significant bit first."""
    return [(value >> i) & 1 for i in range(16)]


# Exact data from the supplied .p section.
patterns = [
    [0, 0, 0, 65280, 61680, 52428, 43690, 0, 52224, 41120],
    [0, 0, 65535, 65280, 61680, 52428, 43690, 0, 65484, 64250],
    [0, 65535, 0, 65280, 61680, 52428, 43690, 0, 13260, 23130],
    [0, 65535, 65535, 65280, 61680, 52428, 43690, 60544, 37740, 23130],
    [65535, 0, 0, 65280, 61680, 52428, 43690, 2254, 14790, 23130],
]


def simulate_packed_table():
    """
    Expand the five .p rows into individual truth-table rows.

    Columns 0..6 are inputs A..G.
    Columns 7..9 are expected outputs.
    """
    rows = []

    for packed in patterns:
        input_bits = [bits16(v) for v in packed[:7]]
        output_bits = [bits16(v) for v in packed[7:10]]

        for i in range(16):
            inputs = [input_bits[j][i] for j in range(7)]
            expected = tuple(output_bits[j][i] for j in range(3))
            rows.append(Row(*inputs, expected=expected))

    return rows


def print_results(rows):
    header = "A B C D E F G | Y1 Y2 Y3 | E1 E2 E3 | OK"
    print(header)
    print("-" * len(header))

    correct = 0

    for r in rows:
        actual = circuit(r.A, r.B, r.C, r.D, r.E, r.F, r.G)
        ok = actual == r.expected
        correct += ok

        print(
            f"{r.A} {r.B} {r.C} {r.D} {r.E} {r.F} {r.G} | "
            f"{actual[0]}  {actual[1]}  {actual[2]} | "
            f"{r.expected[0]}  {r.expected[1]}  {r.expected[2]} | "
            f"{'OK' if ok else 'FAIL'}"
        )

    print()
    print(f"Correct rows: {correct}/{len(rows)}")
    print(f"Fitness: {100.0 * correct / len(rows):.2f}%")


def test_exhaustive_7_input_truth_table():
    """
    Optional: simulate all 2^7 = 128 possible input combinations.
    This is useful for checking the Boolean structure independently
    from the packed .p benchmark representation.
    """
    print("\nExhaustive 7-input simulation:")
    print("A B C D E F G | Y1 Y2 Y3")
    print("-" * 25)

    for n in range(128):
        A = (n >> 6) & 1
        B = (n >> 5) & 1
        C = (n >> 4) & 1
        D = (n >> 3) & 1
        E = (n >> 2) & 1
        F = (n >> 1) & 1
        G = n & 1

        y1, y2, y3 = circuit(A, B, C, D, E, F, G)
        print(f"{A} {B} {C} {D} {E} {F} {G} | {y1}  {y2}  {y3}")


if __name__ == "__main__":
    rows = simulate_packed_table()
    print_results(rows)

    # Uncomment this if you want every one of the 128 input combinations.
    # test_exhaustive_7_input_truth_table()
