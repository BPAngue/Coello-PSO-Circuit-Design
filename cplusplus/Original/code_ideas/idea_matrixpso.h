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
/* File: matrixpso.h                                    */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Header for the matrix circuit           */
/* representation and evaluation library.               */
/********************************************************/

#pragma once

#include <vector>
#include <string>

#include "circuits.h"

namespace PSwarm {

/* Forward-declare Swarm to avoid circular dependency. */
class Swarm;

constexpr unsigned maxBoolExprSize = 4000;

/*------------------------------------------------------*/
/* MatrixDecoder                                         */
/*                                                      */
/* Matrix representation with cumulative source access. */
/*                                                      */
/* Source index layout:                                  */
/*                                                      */
/*   0 ... numInputs-1                                  */
/*       Primary inputs                                  */
/*                                                      */
/*   numInputs ... numInputs+tMat-1                    */
/*       Gate outputs, stored column-major              */
/*                                                      */
/* A gate in column C may reference:                     */
/*                                                      */
/*   - any primary input                                */
/*   - any gate in a column before C                    */
/*                                                      */
/* It may NOT reference a gate in its own column or a    */
/* later column. This preserves the DAG structure.      */
/*------------------------------------------------------*/
struct MatrixDecoder {

    /* ---- Matrix dimensions ---- */

    unsigned numRows = 0;
    unsigned numCols = 0;
    unsigned tMat    = 0;   /* numRows * numCols */


    /* ---- Decoded chromosome ---- */

    /*
     * One entry per matrix cell.
     *
     * input1[cell] and input2[cell] are GLOBAL source
     * indices into the cumulative source pool.
     */
    std::vector<unsigned> input1;
    std::vector<unsigned> input2;
    std::vector<unsigned> gateType;


    /* ---- Evaluation temporaries ---- */

    /*
     * output contains the outputs of the current matrix
     * column.
     *
     * Size = numRows.
     */
    std::vector<unsigned> output;

    /*
     * gateCount marks matrix cells that contribute to
     * the requested circuit outputs.
     */
    std::vector<unsigned> gateCount;

    /*
     * Cumulative packed source pool.
     *
     * Layout:
     *
     *   [0 ... numInputs-1]
     *       Primary inputs
     *
     *   [numInputs ... numInputs+tMat-1]
     *       Matrix gate outputs
     *
     * Unlike the old implementation, this is NEVER
     * overwritten at a column boundary.
     */
    std::vector<unsigned> inTT;


    /* ---- Lifecycle ---- */

    void reserve(const CircuitData& circuit);

    void free();


    /* ---- Core operations ---- */

    /*
     * Decode chromosome M into the matrix representation.
     *
     * Source indices are converted to GLOBAL indices.
     */
    void decode(const std::vector<unsigned>& M,
                const Swarm& swarm,
                const CircuitData& circuit);


    /*
     * Evaluate the circuit against the packed truth table.
     *
     * Each packed word represents up to 32 truth-table rows.
     */
    void evaluate(const std::vector<unsigned>& M,
                  unsigned& numEqual,
                  const Swarm& swarm,
                  const CircuitData& circuit);


    /*
     * Count the non-wire gates that contribute to the
     * requested outputs.
     */
    unsigned countGates(const CircuitData& circuit);


    /*
     * Generate the Boolean expression for an output.
     */
    std::string expression(const std::vector<unsigned>& M,
                           unsigned out,
                           const Swarm& swarm,
                           const CircuitData& circuit);


    /*
     * Print the decoded matrix.
     */
    void printMatrix(const std::vector<unsigned>& M,
                     const Swarm& swarm,
                     const CircuitData& circuit);


private:

    /*
     * Recursively count gates using GLOBAL source indices.
     */
    void countGate(unsigned cell,
                   const CircuitData& circuit);


    /*
     * Recursively construct a Boolean expression from a
     * matrix cell.
     */
    std::string booleanString(unsigned cell,
                              const CircuitData& circuit);


    /*
     * Convert a primary-input index into a printable name.
     *
     * 0 -> A
     * 1 -> B
     * ...
     */
    std::string primaryInputName(unsigned index) const;
};

} // namespace PSwarm