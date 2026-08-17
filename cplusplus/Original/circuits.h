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
/* File: circuits.h                                     */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Header for the circuit truth-table      */
/* loading and I/O library.                             */
/********************************************************/

#pragma once

#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <string>
#include <vector>

namespace PSwarm {

/* -------------------------------------------------------
 * Character sentinels (replaced #define with constexpr)
 * ------------------------------------------------------- */
constexpr char EndLine   = '\n';
constexpr char EndString = '\x0';

/* -------------------------------------------------------
 * Gate type constants (replaced #define with constexpr)
 * ------------------------------------------------------- */
constexpr unsigned AND   = 0;
constexpr unsigned OR    = 1;
constexpr unsigned NOT   = 2;
constexpr unsigned WIRE  = 3;
constexpr unsigned XOR   = 4;
constexpr unsigned NOT1  = 5;
constexpr unsigned WIRE1 = 6;
constexpr unsigned XOR1  = 7;

constexpr unsigned DONTCARE = 2;  /* "Don't care" in the output truth table */

/* Number of truth-table rows packed into a single word by the .plu format. */
constexpr unsigned PLU_CHUNK_BITS = 32u;

/*------------------------------------------------------*/
/* CircuitData: owns all truth-table state that was     */
/* previously scattered across global variables in      */
/* circuits.h.                                          */
/*------------------------------------------------------*/
struct CircuitData {
    unsigned numGates        = 0;   /* Number of available gate types */
    unsigned numInputs       = 0;   /* Inputs in the truth table (.i) */
    unsigned numOutputs      = 0;   /* Outputs in the truth table (.o) */
    unsigned numInstances    = 0;   /* Packed 32-row chunks (.p) */
    unsigned numRowsTT       = 0;   /* Rows in the truth table (2^numInputs) */
    unsigned numTotalOutputs = 0;   /* numOutputs * numRowsTT */

    /* Packed truth-table data: inputTT[instance][inputCol] / outputTT[instance][outputCol]. 
       Each value is a 32-bit word; bit k encodes the truth-table row
       (instance * 32 + k).*/
    std::vector<std::vector<unsigned>> inputTT;   /* Truth-table inputs  */
    std::vector<std::vector<unsigned>> outputTT;  /* Truth-table outputs */

    // /* Build truth-table input combinations */
    // void generateTT();

    /* Number of valid bits in each chunk (32 for all but possibly the last
       chunk, when numRowsTT isn't a multiple of 32). Used to mask off
       padding bits during evaluation so they never count as matches. */
    std::vector<unsigned> chunkValidBits;

    /* Load a .plu-format truth-table block from an already-open file
       stream, positioned at the start of ".i ...". Leaves the stream
       positioned right after ".e" so the caller can keep reading
       further parameters from the same file. */
    void loadTT(std::FILE* pf);

    /* Display the truth table on stdout */
    void printTT() const;

    /* Release truth-table memory (vectors clear themselves;
       provided for explicit lifetime management) */
    void freeMemory();

    /* ---- File-parsing helpers (static: no circuit state needed) ---- */
    static double      readNumber(std::FILE* pf);
    static std::string readString(std::FILE* pf);
    static std::string readToken(std::FILE* pf);
    static unsigned readUnsigned(std::FILE* pf);
    static bool        isDigit(char c);
};

} // namespace PSwarm
