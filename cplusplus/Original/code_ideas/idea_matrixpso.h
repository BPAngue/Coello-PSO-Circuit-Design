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

/* Forward-declare Swarm so MatrixDecoder can reference it
   in function signatures without a circular include. */
class Swarm;

constexpr unsigned maxBoolExprSize = 4000;

/*------------------------------------------------------*/
/* MatrixDecoder: owns all matrix-decoding state that   */
/* was previously scattered across globals in           */
/* matrixpso.h.  Requires a CircuitData reference and   */
/* a Swarm reference for functions that need PSO or     */
/* circuit parameters.                                  */
/*                                                      */
/* Source Index Layout:                                 */              
/*     0 ... numInputs-1: primary inputs                */                                   
/*     numInputs ... numInputs+tMat-1: gate outputs,    */
/*                                   stored column-major*/
/*                                                      */
/* A gate in column C may reference:                    */                                
/*    - any primary input                               */                   
/*    - any gate in a column before C                   */
/* It may NOT reference a gate in its own column        */
/*     or a later column. This preserves the            */                                      
/*------------------------------------------------------*/
struct MatrixDecoder {

    /* Circuit matrix dimensions */
    unsigned numRows = 0;
    unsigned numCols = 0;
    unsigned tMat    = 0;   /* numRows * numCols */

    /* Per-cell decoded values (rebuilt each evaluation) */
    std::vector<unsigned> input1;
    std::vector<unsigned> input2;
    std::vector<unsigned> gateType;

    /* Evaluation temporaries */
    std::vector<unsigned> output; // contains outputs of the current matrix column; size = numRows
    std::vector<unsigned> gateCount; // marks matrix cells that contribute to the requested circuit outputs
   //  std::vector<std::vector<unsigned>> inTT;   /* [numRowsTT][numRows] */

   /* Rolling register file: one packed 32-bit word per matrix row,
       reused for every chunk of 32 truth-table rows during evaluate().
       Replaces the old [numRowsTT][numRows] table -- since evaluation
       now processes one packed chunk at a time (mirroring CGP++), we
       only ever need one row's worth of registers live at once. 
       
    EDIT: 
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
     * overwritten at a column boundary.*/
    std::vector<unsigned> inTT;

    /* ---- Lifecycle ---- */

    /* Allocate working vectors sized for the given circuit */
    void reserve(const CircuitData& circuit);

    /* Release all working vectors */
    void free();

    /* ---- Core operations ---- */

    /* Decode chromosome M into input1, input2, gateType arrays.
       Needs nVar / representation / bitVariable from Swarm,
       and numGates from CircuitData. */
    void decode(const std::vector<unsigned>& M,
                const Swarm& swarm,
                const CircuitData& circuit);

   //  /* Evaluate circuit against all truth-table rows.
   //     Sets numEqual to the number of matching outputs. */
   //  void evaluate(const std::vector<unsigned>& M,
   //                unsigned& numEqual,
   //                const Swarm& swarm,
   //                const CircuitData& circuit);

   /* Evaluate the circuit against the packed truth table (circuit.inputTT /
       circuit.outputTT, one chunk of 32 rows per iteration). Sets numEqual
       to the number of matching (row, output) pairs across all logical
       rows -- NOT number of chunks, and padding bits in a partial last
       chunk are masked out via circuit.chunkValidBits so they never count
       as spurious matches. */
   void evaluate(const std::vector<unsigned>& M,
                  unsigned& numEqual,
                  const Swarm& swarm,
                  const CircuitData& circuit);

    /* Count the gates actually used in the last decoded solution */
    unsigned countGates(const CircuitData& circuit);

    /* Get the Boolean expression for output 'out' as a string */
    std::string expression(const std::vector<unsigned>& M,
                           unsigned out,
                           const Swarm& swarm,
                           const CircuitData& circuit);

    /* Print the decoded gate matrix to stdout */
    void printMatrix(const std::vector<unsigned>& M,
                     const Swarm& swarm,
                     const CircuitData& circuit);

private:
    /* Recursive helpers */
    // void countGate(unsigned out, int col, const CircuitData& circuit);
    void countGate(unsigned cell, const CircuitData& circuit);
    // std::string booleanString(unsigned cell, int col, const CircuitData& circuit);
    std::string booleanString(unsigned cell, const CircuitData& circuit);
    std::string primaryInputName(unsigned index) const; // Returns the name of a primary input by its index 0 -> A, 1 -> B, etc.
};

} // namespace PSwarm