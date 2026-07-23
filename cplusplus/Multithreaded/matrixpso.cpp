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
#include <string>
#include <algorithm>

namespace PSwarm {

/* -------------------------------------------------------
 * reserve
 * Allocates working vectors sized for the given circuit.
 * Replaces the original malloc-based reserveMatrixMemory.
 * ------------------------------------------------------- */
void MatrixDecoder::reserve(const CircuitData& circuit)
{
    input1.assign(tMat, 0u);
    input2.assign(tMat, 0u);
    gateType.assign(tMat, 0u);
    output.assign(numRows, 0u);
    gateCount.assign(tMat, 0u);
    inTT.assign(circuit.numRowsTT, std::vector<unsigned>(numRows, 0u));
}

/* -------------------------------------------------------
 * free
 * Releases all working vectors.
 * Original freeMatrixMemory had a bug (malloc instead of
 * free in the inTT loop) — fixed here by simply clearing.
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
 * evaluate
 * Evaluates the circuit for every row of the truth table
 * and counts how many outputs match the target.
 * ------------------------------------------------------- */
void MatrixDecoder::evaluate(const std::vector<unsigned>& M,
                             unsigned& numEqual,
                             const Swarm& swarm,
                             const CircuitData& circuit)
{
    numEqual = 0;
    decode(M, swarm, circuit);

    /* Initialise the rolling input/output table */
    for (unsigned i = 0; i < circuit.numRowsTT; ++i)
        for (unsigned j = 0; j < numRows; ++j)
            inTT[i][j] = (j < circuit.numInputs)
                            ? circuit.inputTT[i][j]
                            : inTT[i][j - circuit.numInputs];

    /* Simulate every column of the gate matrix */
    for (unsigned i = 0; i < circuit.numRowsTT; ++i) {
        unsigned in = 0;
        for (unsigned j = 0; j < tMat; ++j) {
            switch (gateType[j]) {
                case AND:
                    output[in] = inTT[i][input1[j]] & inTT[i][input2[j]]; break;
                case OR:
                    output[in] = inTT[i][input1[j]] | inTT[i][input2[j]]; break;
                case NOT:
                case NOT1:
                    output[in] = inTT[i][input1[j]] ? 0u : 1u; break;
                case WIRE:
                case WIRE1:
                    output[in] = inTT[i][input1[j]]; break;
                case XOR:
                case XOR1:
                    output[in] = inTT[i][input1[j]] ^ inTT[i][input2[j]]; break;
            }

            /* At the end of each column, feed outputs back as new inputs */
            if (!((j + 1) % numRows)) {
                in = 0;
                for (unsigned k = 0; k < numRows; ++k)
                    inTT[i][k] = output[k];
            } else {
                ++in;
            }
        }

        /* Count matching outputs */
        for (unsigned j = 0; j < circuit.numOutputs; ++j)
            if (output[j] == circuit.outputTT[i][j]) ++numEqual;
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
