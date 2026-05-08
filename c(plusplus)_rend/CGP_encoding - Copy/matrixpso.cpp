/********************************************************/
/*                  CINVESTAV - IPN                     */
/*        Department of Electrical Engineering          */
/*                 Computing Section                    */
/*                                                      */
/*               Evolutionary Computation               */
/*                                                      */
/*                Erika Hernandez Luna                  */
/*         eluna@computacion.cs.cinvestav.mx            */
/*                   August 2, 2003                     */
/*         Converted to C++ - April 2026                */
/*                                                      */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* File: matrixpso.cpp                                  */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Library for handling the matrix that    */
/* represents a combinational logic circuit.            */
/********************************************************/

#include "matrixpso.h"
#include "psomatrixcircuit.h"   /* Full Swarm definition needed here */
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <string>
#include <algorithm>

/* Portable 64-bit popcount.
   __builtin_popcountll is supported by GCC, Clang and MinGW. */
#if defined(__GNUC__) || defined(__clang__)
  #define PSO_POPCOUNT64(x) __builtin_popcountll(x)
#else
  static inline unsigned pso_popcount64(std::uint64_t x)
  {
      x = x - ((x >> 1) & 0x5555555555555555ULL);
      x = (x & 0x3333333333333333ULL) +
          ((x >> 2) & 0x3333333333333333ULL);
      x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0fULL;
      return static_cast<unsigned>((x * 0x0101010101010101ULL) >> 56);
  }
  #define PSO_POPCOUNT64(x) pso_popcount64(x)
#endif

namespace PSwarm {

/* -------------------------------------------------------
 * reserve
 * Allocates working vectors sized for the given circuit.
 * inTT / output are flat buffers of numRows wires x
 * numWords 64-bit words each (bit-packed, 64 truth-table
 * rows per word).
 * ------------------------------------------------------- */
void MatrixDecoder::reserve(const CircuitData& circuit)
{
    input1.assign(tMat, 0u);
    input2.assign(tMat, 0u);
    gateType.assign(tMat, 0u);
    gateCount.assign(tMat, 0u);

    const std::size_t buf = static_cast<std::size_t>(numRows) * circuit.numWords;
    inTT.assign(buf, 0ULL);
    output.assign(buf, 0ULL);
}

/* -------------------------------------------------------
 * free
 * Releases all working vectors.
 * ------------------------------------------------------- */
void MatrixDecoder::free()
{
    input1.clear();
    input2.clear();
    gateType.clear();
    output.clear();
    gateCount.clear();
    inTT.clear();
}

/* -------------------------------------------------------
 * decode
 * Converts the flat chromosome vector M into per-cell
 * input1, input2, gateType arrays.
 * ------------------------------------------------------- */
void MatrixDecoder::decode(const std::vector<unsigned>& M,
                           const Swarm& swarm,
                           const CircuitData& circuit)
{
    unsigned j = 0;   /* bit index used only for binary representation */
    unsigned in = 0;  /* gate cell index */

    for (unsigned i = 0; i < swarm.nVar; ++i) {
        unsigned num = 0;

        if (swarm.representation != BINARY) {
            num = M[i];
        } else {
            for (int k = static_cast<int>(swarm.bitVariable[i]) - 1; k >= 0; --k, ++j)
                if (M[j]) num += static_cast<unsigned>(std::pow(2.0, static_cast<double>(k)));
        }

        switch (i % 3) {
            case 0: input1[in]       = num % numRows;           break;
            case 1: input2[in]       = num % numRows;           break;
            case 2: gateType[in++]   = num % circuit.numGates;  break;
        }
    }
}

/* -------------------------------------------------------
 * evaluate  (bit-parallel)
 *
 * Each wire's value across all 2^numInputs truth-table
 * rows is packed into numWords 64-bit words.  A single
 * bitwise AND/OR/XOR/NOT on one word evaluates 64 rows
 * at once, turning the inner loop from
 *     numRowsTT * tMat
 * into
 *     numWords  * tMat
 * operations - typically a 32x-64x speed-up for circuits
 * with 9+ inputs.
 *
 * The observable fitness (numEqual) is identical to the
 * row-at-a-time implementation: it counts the number of
 * (row, output) pairs whose simulated bit matches the
 * target, ignoring DONTCARE entries via the precomputed
 * care mask in CircuitData.
 * ------------------------------------------------------- */
void MatrixDecoder::evaluate(const std::vector<unsigned>& M,
                             unsigned& numEqual,
                             const Swarm& swarm,
                             const CircuitData& circuit)
{
    numEqual = 0;
    decode(M, swarm, circuit);

    const unsigned W = circuit.numWords;
    if (W == 0) return;   /* No rows -> nothing to evaluate */

    /* --------- Initialise the rolling wire buffer ---------
       Matches the original behaviour: wire j is seeded with
       primary input (j mod numInputs).  In packed form this
       is just a copy of numWords words per wire. */
    for (unsigned wire = 0; wire < numRows; ++wire) {
        const unsigned src = wire % circuit.numInputs;
        const uint64_t* srcPtr = circuit.inputTTPacked[src].data();
        uint64_t*       dstPtr = inTT.data() + static_cast<std::size_t>(wire) * W;
        for (unsigned w = 0; w < W; ++w) dstPtr[w] = srcPtr[w];
    }

    /* --------- Simulate the gate matrix column by column ---------
       Every cell writes W words; at the end of each column the
       freshly computed outputs are fed back into inTT as the inputs
       to the next column. */
    unsigned in = 0;
    for (unsigned j = 0; j < tMat; ++j) {
        const uint64_t* a = inTT.data() + static_cast<std::size_t>(input1[j]) * W;
        const uint64_t* b = inTT.data() + static_cast<std::size_t>(input2[j]) * W;
        uint64_t*       o = output.data() + static_cast<std::size_t>(in)        * W;

        switch (gateType[j]) {
            case AND:
                for (unsigned w = 0; w < W; ++w) o[w] = a[w] & b[w];
                break;
            case OR:
                for (unsigned w = 0; w < W; ++w) o[w] = a[w] | b[w];
                break;
            case NOT:
            case NOT1:
                for (unsigned w = 0; w < W; ++w) o[w] = ~a[w];
                break;
            case WIRE:
            case WIRE1:
                for (unsigned w = 0; w < W; ++w) o[w] = a[w];
                break;
            case XOR:
            case XOR1:
                for (unsigned w = 0; w < W; ++w) o[w] = a[w] ^ b[w];
                break;
        }

        if (!((j + 1) % numRows)) {
            /* End of column: copy the whole output buffer into inTT */
            const std::size_t total = static_cast<std::size_t>(numRows) * W;
            for (std::size_t k = 0; k < total; ++k) inTT[k] = output[k];
            in = 0;
        } else {
            ++in;
        }
    }

    /* --------- Count matching output bits ---------
       For each circuit output j, output[j] holds the simulated
       bits produced by the final column.  A row matches when
           (~(sim ^ target)) & care
       has the bit set.  popcount over all words gives numEqual.
       Bits beyond numRowsTT are zero in outputTTCare, so they
       never contribute. */
    for (unsigned j = 0; j < circuit.numOutputs; ++j) {
        const uint64_t* sim    = output.data() + static_cast<std::size_t>(j) * W;
        const uint64_t* target = circuit.outputTTTarget[j].data();
        const uint64_t* care   = circuit.outputTTCare[j].data();
        for (unsigned w = 0; w < W; ++w) {
            const uint64_t match = (~(sim[w] ^ target[w])) & care[w];
            numEqual += PSO_POPCOUNT64(match);
        }
    }
}

/* -------------------------------------------------------
 * countGates
 * Returns the number of non-wire gates actually used by
 * the last decoded solution.
 * ------------------------------------------------------- */
unsigned MatrixDecoder::countGates(const CircuitData& circuit)
{
    std::fill(gateCount.begin(), gateCount.end(), 0u);

    unsigned startOutput = numRows * (numCols - 1);
    for (unsigned i = startOutput; i < startOutput + circuit.numOutputs; ++i)
        countGate(i, static_cast<int>(numCols) - 1, circuit);

    unsigned gates = 0;
    for (unsigned i = 0; i < tMat; ++i)
        if (gateCount[i]) ++gates;

    return gates;
}

/* -------------------------------------------------------
 * countGate (private, recursive)
 * Marks each gate cell that contributes to an output.
 * ------------------------------------------------------- */
void MatrixDecoder::countGate(unsigned out, int col, const CircuitData& circuit)
{
    if (col < 0) return;

    /* Mark non-wire, non-trivial gates */
    if (gateType[out] != WIRE && gateType[out] != WIRE1) {
        if (gateType[out] != AND && gateType[out] != OR) {
            gateCount[out] = 1;
        } else if (input1[out] != input2[out]) {
            gateCount[out] = 1;
        }
    }

    unsigned prev = static_cast<unsigned>(numRows * (col - 1));
    countGate(prev + input1[out], col - 1, circuit);
    countGate(prev + input2[out], col - 1, circuit);
}

/* -------------------------------------------------------
 * expression
 * Returns the Boolean expression string for output 'out'
 * of the last decoded matrix.
 * ------------------------------------------------------- */
std::string MatrixDecoder::expression(const std::vector<unsigned>& M,
                                      unsigned out,
                                      const Swarm& swarm,
                                      const CircuitData& circuit)
{
    decode(M, swarm, circuit);
    return booleanString(numRows * (numCols - 1) + out,
                         static_cast<int>(numCols) - 1,
                         circuit);
}

/* -------------------------------------------------------
 * booleanString (private, recursive)
 * Builds the Boolean expression string for one cell.
 * Fixed bug: original had `!= NOT && != NOT` (duplicate);
 * corrected to `!= NOT && != NOT1`.
 * ------------------------------------------------------- */
std::string MatrixDecoder::booleanString(unsigned cell, int col,
                                         const CircuitData& circuit)
{
    std::string str;

    switch (gateType[cell]) {
        case AND:          str = "(AND1 ";  break;
        case OR:           str = "(OR1 ";   break;
        case NOT:
        case NOT1:         str = "(NOT1 ";  break;
        case XOR:
        case XOR1:         str = "(XOR1 ";  break;
        case WIRE:
        case WIRE1:        str = "(WIRE ";  break;
        default:           str = "(??? ";   break;
    }

    bool isBinary = (gateType[cell] != WIRE  && gateType[cell] != WIRE1 &&
                     gateType[cell] != NOT   && gateType[cell] != NOT1);

    if (col <= 0) {
        /* Base case: leaf node refers to a primary input */
        str += static_cast<char>(input1[cell] % circuit.numInputs + 'A');
        if (isBinary) {
            str += ' ';
            str += static_cast<char>(input2[cell] % circuit.numInputs + 'A');
        }
    } else {
        unsigned new1 = static_cast<unsigned>(numRows * (col - 1)) + input1[cell];
        unsigned new2 = static_cast<unsigned>(numRows * (col - 1)) + input2[cell];

        str += booleanString(new1, col - 1, circuit);
        if (isBinary) {
            str += ' ';
            str += booleanString(new2, col - 1, circuit);
        }
    }

    str += ')';
    return str;
}

/* -------------------------------------------------------
 * printMatrix
 * Prints the decoded gate matrix to stdout.
 * Replaces fcvt() (POSIX) with std::to_string().
 * ------------------------------------------------------- */
void MatrixDecoder::printMatrix(const std::vector<unsigned>& M,
                                const Swarm& swarm,
                                const CircuitData& circuit)
{
    std::printf("\n");
    decode(M, swarm, circuit);

    for (unsigned i = 0; i < numRows; ++i) {
        std::string row;
        for (unsigned j = 0; j < numCols; ++j) {
            unsigned idx = j * numRows + i;
            switch (gateType[idx]) {
                case AND:   row += "AND(";  break;
                case OR:    row += "OR(";   break;
                case NOT:
                case NOT1:  row += "NOT(";  break;
                case WIRE:
                case WIRE1: row += "WIRE("; break;
                case XOR:
                case XOR1:  row += "XOR(";  break;
            }
            row += std::to_string(input1[idx]);
            row += ' ';
            row += std::to_string(input2[idx]);
            row += ')';
            if (j != numCols - 1) row += ',';
        }
        std::printf("%s\n", row.c_str());
    }
}

} // namespace PSwarm
