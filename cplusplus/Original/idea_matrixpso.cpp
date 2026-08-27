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
#include "psomatrixcircuit.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <algorithm>
#include <bit>

namespace PSwarm {

/* =======================================================
 * reserve
 *
 * Allocate all matrix/evaluation working memory.
 *
 * IMPORTANT:
 * inTT is now:
 *
 *     numInputs + tMat
 *
 * rather than just numRows.
 *
 * The first numInputs entries contain primary inputs.
 * The remaining tMat entries contain gate outputs.
 * ======================================================= */
void MatrixDecoder::reserve(const CircuitData& circuit)
{
    input1.assign(tMat, 0u);
    input2.assign(tMat, 0u);
    gateType.assign(tMat, 0u);

    output.assign(numRows, 0u);
    gateCount.assign(tMat, 0u);

    /*
     * Cumulative source pool:
     *
     * primary inputs + every matrix cell output.
     */
    inTT.assign(circuit.numInputs + tMat, 0u);
}


/* =======================================================
 * free
 * ======================================================= */
void MatrixDecoder::free()
{
    input1.clear();
    input2.clear();
    gateType.clear();

    output.clear();
    gateCount.clear();
    inTT.clear();
}


/* =======================================================
 * decode
 *
 * Converts the flat chromosome into:
 *
 *     input1[cell]
 *     input2[cell]
 *     gateType[cell]
 *
 * SOURCE INDEXING
 *
 * Primary inputs:
 *
 *     0 ... numInputs-1
 *
 * Gate outputs:
 *
 *     numInputs + cell
 *
 * For a gate in column C:
 *
 *     legal sources =
 *
 *     primary inputs
 *     +
 *     all gates from columns 0 ... C-1
 *
 * Therefore:
 *
 *     sourceCount = numInputs + C * numRows
 *
 * A gate cannot reference:
 *
 *     - itself
 *     - another gate in the same column
 *     - a future column
 *
 * This preserves the DAG property.
 * ======================================================= */
void MatrixDecoder::decode(const std::vector<unsigned>& M,
                           const Swarm& swarm,
                           const CircuitData& circuit)
{
    unsigned bitIndex = 0;

    for (unsigned cell = 0; cell < tMat; ++cell) {

        /*
         * Determine which matrix column this cell belongs to.
         *
         * Matrix cells are stored column-major:
         *
         *     column * numRows + row
         */
        const unsigned column = cell / numRows;

        /*
         * Number of source values that already exist before
         * this column.
         *
         * Example:
         *
         * numInputs = 5
         * numRows   = 5
         *
         * Column 0:
         *     sourceCount = 5
         *
         * Column 1:
         *     sourceCount = 10
         *
         * Column 2:
         *     sourceCount = 15
         */
        const unsigned sourceCount =
            circuit.numInputs + column * numRows;

        /*
         * Safety check.
         *
         * Column 0 must always have at least the primary
         * inputs available.
         */
        if (sourceCount == 0)
            continue;


        /* -----------------------------------------------
         * input1
         * ----------------------------------------------- */

        unsigned num = 0;

        if (swarm.representation != BINARY) {

            num = M[cell * 3];

        } else {

            const unsigned variable = cell * 3;

            for (int k =
                     static_cast<int>(swarm.bitVariable[variable]) - 1;
                 k >= 0;
                 --k, ++bitIndex)
            {
                if (M[bitIndex])
                    num +=
                        static_cast<unsigned>(
                            std::pow(
                                2.0,
                                static_cast<double>(k)
                            )
                        );
            }
        }

        input1[cell] = num % sourceCount;


        /* -----------------------------------------------
         * input2
         * ----------------------------------------------- */

        num = 0;

        if (swarm.representation != BINARY) {

            num = M[cell * 3 + 1];

        } else {

            const unsigned variable = cell * 3 + 1;

            for (int k =
                     static_cast<int>(swarm.bitVariable[variable]) - 1;
                 k >= 0;
                 --k, ++bitIndex)
            {
                if (M[bitIndex])
                    num +=
                        static_cast<unsigned>(
                            std::pow(
                                2.0,
                                static_cast<double>(k)
                            )
                        );
            }
        }

        input2[cell] = num % sourceCount;


        /* -----------------------------------------------
         * gate type
         * ----------------------------------------------- */

        num = 0;

        if (swarm.representation != BINARY) {

            num = M[cell * 3 + 2];

        } else {

            const unsigned variable = cell * 3 + 2;

            for (int k =
                     static_cast<int>(swarm.bitVariable[variable]) - 1;
                 k >= 0;
                 --k, ++bitIndex)
            {
                if (M[bitIndex])
                    num +=
                        static_cast<unsigned>(
                            std::pow(
                                2.0,
                                static_cast<double>(k)
                            )
                        );
            }
        }

        gateType[cell] = num % circuit.numGates;
    }
}


/* =======================================================
 * evaluate
 *
 * Bit-parallel evaluation of the circuit.
 *
 * OLD BEHAVIOUR:
 *
 *     inTT = current column only
 *
 * NEW BEHAVIOUR:
 *
 *     inTT contains:
 *
 *       [primary inputs]
 *       [column 0 outputs]
 *       [column 1 outputs]
 *       [column 2 outputs]
 *       ...
 *
 * Every output is written once and retained.
 * ======================================================= */
void MatrixDecoder::evaluate(const std::vector<unsigned>& M,
                             unsigned& numEqual,
                             const Swarm& swarm,
                             const CircuitData& circuit)
{
    numEqual = 0;

    decode(M, swarm, circuit);


    /* ===================================================
     * Evaluate every packed 32-row chunk.
     * =================================================== */
    for (unsigned c = 0; c < circuit.numInstances; ++c) {

        /* -----------------------------------------------
         * Load primary inputs.
         *
         * These occupy:
         *
         *     inTT[0 ... numInputs-1]
         *
         * They remain available for the entire evaluation
         * of this packed chunk.
         * ----------------------------------------------- */
        for (unsigned j = 0; j < circuit.numInputs; ++j) {
            inTT[j] = circuit.inputTT[c][j];
        }


        /* -----------------------------------------------
         * Evaluate every matrix cell in column-major order.
         * ----------------------------------------------- */
        for (unsigned cell = 0; cell < tMat; ++cell) {

            /*
             * Row within the current matrix column.
             */
            const unsigned row = cell % numRows;


            /*
             * Evaluate the gate using GLOBAL source
             * indices.
             */
            switch (gateType[cell]) {

                case AND:
                    output[row] =
                        inTT[input1[cell]]
                        &
                        inTT[input2[cell]];
                    break;


                case OR:
                    output[row] =
                        inTT[input1[cell]]
                        |
                        inTT[input2[cell]];
                    break;


                case NOT:
                case NOT1:
                    /*
                     * Bitwise complement operates on all
                     * packed truth-table rows simultaneously.
                     */
                    output[row] =
                        ~inTT[input1[cell]];
                    break;


                case WIRE:
                case WIRE1:
                    output[row] =
                        inTT[input1[cell]];
                    break;


                case XOR:
                case XOR1:
                    output[row] =
                        inTT[input1[cell]]
                        ^
                        inTT[input2[cell]];
                    break;


                default:
                    /*
                     * Should never happen because gateType
                     * is bounded by circuit.numGates.
                     */
                    output[row] = 0u;
                    break;
            }


            /*
             * Store this gate's result permanently in the
             * cumulative source pool.
             *
             * Primary inputs occupy:
             *
             *     0 ... numInputs-1
             *
             * Gate cell 0 occupies:
             *
             *     numInputs
             *
             * Gate cell 1:
             *
             *     numInputs + 1
             *
             * etc.
             */
            inTT[circuit.numInputs + cell] = output[row];
        }


        /* -----------------------------------------------
         * Compare final-column outputs against target.
         *
         * The circuit outputs are expected to be located
         * in the last matrix column, exactly as in the
         * original representation.
         * ----------------------------------------------- */

        const unsigned validBits =
            circuit.chunkValidBits[c];

        const unsigned validMask =
            (validBits >= 32u)
                ? 0xFFFFFFFFu
                : ((1u << validBits) - 1u);


        const unsigned finalColumnStart =
            numRows * (numCols - 1);


        for (unsigned j = 0;
             j < circuit.numOutputs;
             ++j)
        {
            const unsigned cell =
                finalColumnStart + j;

            const unsigned diff =
                (inTT[circuit.numInputs + cell]
                 ^
                 circuit.outputTT[c][j])
                &
                validMask;

            numEqual +=
                validBits
                -
                static_cast<unsigned>(
                    __builtin_popcount(diff)
                );
        }
    }
}


/* =======================================================
 * countGates
 *
 * Find all matrix gates that contribute to the final
 * outputs.
 *
 * This version understands GLOBAL source indices.
 * ======================================================= */
unsigned MatrixDecoder::countGates(const CircuitData& circuit)
{
    std::fill(
        gateCount.begin(),
        gateCount.end(),
        0u
    );


    /*
     * Final matrix column.
     */
    const unsigned startOutput =
        numRows * (numCols - 1);


    /*
     * Trace every circuit output backwards.
     */
    for (unsigned i = 0;
         i < circuit.numOutputs;
         ++i)
    {
        const unsigned cell =
            startOutput + i;

        if (cell < tMat)
            countGate(cell, circuit);
    }


    unsigned gates = 0;

    for (unsigned i = 0; i < tMat; ++i) {
        if (gateCount[i])
            ++gates;
    }

    return gates;
}


/* =======================================================
 * countGate
 *
 * Recursively traces a matrix cell's dependencies.
 *
 * GLOBAL SOURCE INDEX:
 *
 *     source < numInputs
 *
 * means:
 *
 *     primary input
 *
 * Otherwise:
 *
 *     source >= numInputs
 *
 * means:
 *
 *     matrix gate
 *
 * and therefore:
 *
 *     gateCell = source - numInputs
 * ======================================================= */
void MatrixDecoder::countGate(unsigned cell,
                              const CircuitData& circuit)
{
    if (cell >= tMat)
        return;


    /*
     * Prevent repeatedly walking the same gate.
     *
     * This is particularly useful now that a gate may be
     * referenced by many later columns.
     */
    if (gateCount[cell] == 2u)
        return;


    /*
     * Mark as visited.
     *
     * 1 = visited
     * 2 = counted gate
     *
     * We temporarily use 1 while tracing.
     */
    gateCount[cell] = 1u;


    /*
     * Count this cell if it represents a real gate.
     *
     * Preserve the original optimization behaviour:
     *
     * - WIRE does not count
     * - AND/OR with identical inputs do not count
     * - NOT/XOR count
     */
    bool countsAsGate = false;

    switch (gateType[cell]) {

        case WIRE:
        case WIRE1:
            countsAsGate = false;
            break;

        case AND:
        case OR:
            countsAsGate =
                (input1[cell] != input2[cell]);
            break;

        default:
            countsAsGate = true;
            break;
    }


    if (countsAsGate)
        gateCount[cell] = 2u;


    /*
     * Follow input 1.
     */
    const unsigned source1 = input1[cell];

    if (source1 >= circuit.numInputs) {

        const unsigned sourceCell =
            source1 - circuit.numInputs;

        countGate(sourceCell, circuit);
    }


    /*
     * Follow input 2 for binary gates.
     */
    const bool isBinary =
        gateType[cell] != WIRE
        &&
        gateType[cell] != WIRE1
        &&
        gateType[cell] != NOT
        &&
        gateType[cell] != NOT1;


    if (isBinary) {

        const unsigned source2 = input2[cell];

        if (source2 >= circuit.numInputs) {

            const unsigned sourceCell =
                source2 - circuit.numInputs;

            countGate(sourceCell, circuit);
        }
    }
}


/* =======================================================
 * expression
 * ======================================================= */
std::string MatrixDecoder::expression(
    const std::vector<unsigned>& M,
    unsigned out,
    const Swarm& swarm,
    const CircuitData& circuit)
{
    decode(M, swarm, circuit);


    const unsigned finalColumnStart =
        numRows * (numCols - 1);


    const unsigned cell =
        finalColumnStart + out;


    if (cell >= tMat)
        return "(INVALID)";


    return booleanString(cell, circuit);
}


/* =======================================================
 * primaryInputName
 *
 * Converts:
 *
 *     0 -> A
 *     1 -> B
 *     ...
 *
 * For inputs beyond Z, use a more explicit form.
 * ======================================================= */
std::string MatrixDecoder::primaryInputName(unsigned index) const
{
    if (index < 26) {

        std::string result;

        result +=
            static_cast<char>('A' + index);

        return result;
    }


    return "I" + std::to_string(index);
}


/* =======================================================
 * booleanString
 *
 * Recursively follows GLOBAL source indices.
 * ======================================================= */
std::string MatrixDecoder::booleanString(
    unsigned cell,
    const CircuitData& circuit)
{
    if (cell >= tMat)
        return "(INVALID)";


    std::string str;


    /* -----------------------------------------------
     * Gate name
     * ----------------------------------------------- */
    switch (gateType[cell]) {

        case AND:
            str = "(AND1 ";
            break;

        case OR:
            str = "(OR1 ";
            break;

        case NOT:
        case NOT1:
            str = "(NOT1 ";
            break;

        case XOR:
        case XOR1:
            str = "(XOR1 ";
            break;

        case WIRE:
        case WIRE1:
            str = "(WIRE ";
            break;

        default:
            str = "(??? ";
            break;
    }


    const bool isBinary =
        gateType[cell] != WIRE
        &&
        gateType[cell] != WIRE1
        &&
        gateType[cell] != NOT
        &&
        gateType[cell] != NOT1;


    /* -----------------------------------------------
     * Input 1
     * ----------------------------------------------- */
    const unsigned source1 = input1[cell];

    if (source1 < circuit.numInputs) {

        /*
         * Primary input.
         */
        str += primaryInputName(source1);

    } else {

        /*
         * Earlier matrix gate.
         */
        const unsigned sourceCell =
            source1 - circuit.numInputs;

        str += booleanString(sourceCell, circuit);
    }


    /* -----------------------------------------------
     * Input 2
     * ----------------------------------------------- */
    if (isBinary) {

        str += ' ';

        const unsigned source2 = input2[cell];

        if (source2 < circuit.numInputs) {

            str += primaryInputName(source2);

        } else {

            const unsigned sourceCell =
                source2 - circuit.numInputs;

            str += booleanString(sourceCell, circuit);
        }
    }


    str += ')';

    return str;
}


/* =======================================================
 * printMatrix
 *
 * Prints GLOBAL source indices.
 *
 * Example:
 *
 *     AND(0 1)
 *     XOR(5 2)
 *
 * Here:
 *
 *     0,1,2,... = primary inputs
 *     numInputs+ = matrix gate outputs
 * ======================================================= */
void MatrixDecoder::printMatrix(
    const std::vector<unsigned>& M,
    const Swarm& swarm,
    const CircuitData& circuit)
{
    std::printf("\n");

    decode(M, swarm, circuit);


    for (unsigned row = 0;
         row < numRows;
         ++row)
    {
        std::string matrixRow;


        for (unsigned column = 0;
             column < numCols;
             ++column)
        {
            const unsigned cell =
                column * numRows + row;


            switch (gateType[cell]) {

                case AND:
                    matrixRow += "AND(";
                    break;

                case OR:
                    matrixRow += "OR(";
                    break;

                case NOT:
                case NOT1:
                    matrixRow += "NOT(";
                    break;

                case WIRE:
                case WIRE1:
                    matrixRow += "WIRE(";
                    break;

                case XOR:
                case XOR1:
                    matrixRow += "XOR(";
                    break;

                default:
                    matrixRow += "???(";
                    break;
            }


            matrixRow +=
                std::to_string(input1[cell]);

            matrixRow += ' ';

            matrixRow +=
                std::to_string(input2[cell]);

            matrixRow += ')';


            if (column != numCols - 1)
                matrixRow += ',';
        }


        std::printf(
            "%s\n",
            matrixRow.c_str()
        );
    }
}

} // namespace PSwarm