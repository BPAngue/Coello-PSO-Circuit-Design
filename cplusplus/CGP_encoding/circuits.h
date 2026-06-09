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
#include <cstdint>
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

/*------------------------------------------------------*/
/* CircuitData: owns all truth-table state that was     */
/* previously scattered across global variables in      */
/* circuits.h.                                          */
/*------------------------------------------------------*/
struct CircuitData {
    unsigned numGates        = 0;   /* Number of available gate types */
    unsigned numInputs       = 0;   /* Inputs in the truth table */
    unsigned numOutputs      = 0;   /* Outputs in the truth table */
    unsigned numRowsTT       = 0;   /* Rows in the truth table (2^numInputs) */
    unsigned numTotalOutputs = 0;   /* numOutputs * numRowsTT */

    std::vector<std::vector<unsigned>> inputTT;   /* Truth-table inputs  */
    std::vector<std::vector<unsigned>> outputTT;  /* Truth-table outputs */

    /* ---- Bit-packed truth table (64 rows per word) ----
       Built once by packTT() and used by the bit-parallel
       evaluator in MatrixDecoder.  Each word holds the bit
       for 64 consecutive truth-table rows. */
    unsigned numWords = 0;   /* (numRowsTT + 63) / 64 */
    std::vector<std::vector<uint64_t>> inputTTPacked;   /* [numInputs][numWords]  */
    std::vector<std::vector<uint64_t>> outputTTTarget;  /* [numOutputs][numWords] */
    std::vector<std::vector<uint64_t>> outputTTCare;    /* [numOutputs][numWords] */

    /* Build truth-table input combinations */
    void generateTT();

    /* Load truth-table and circuit meta-data from file */
    void loadTT(std::FILE* pf);

    /* Build the bit-packed representation from inputTT / outputTT.
       Called automatically by loadTT(). */
    void packTT();

    /* Display the truth table on stdout */
    void printTT() const;

    /* Release truth-table memory (vectors clear themselves;
       provided for explicit lifetime management) */
    void freeMemory();

    /* ---- File-parsing helpers (static: no circuit state needed) ---- */
    static double      readNumber(std::FILE* pf);
    static std::string readString(std::FILE* pf);
    static bool        isDigit(char c);
};

} // namespace PSwarm
