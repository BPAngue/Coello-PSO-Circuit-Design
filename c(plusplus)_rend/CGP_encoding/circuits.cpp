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
/* File: circuits.cpp                                   */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Library to load into memory a file with */
/* the design of the logic circuit to optimise.         */
/********************************************************/

#include "circuits.h"
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cstdint>

namespace PSwarm {

/* -------------------------------------------------------
 * generateTT
 * Builds the input combinations for the truth table.
 * Replaces the original malloc-based approach with nested
 * vectors that manage their own memory.
 * ------------------------------------------------------- */
void CircuitData::generateTT()
{
    inputTT.assign(numRowsTT, std::vector<unsigned>(numInputs, 0u));

    for (unsigned i = 0; i < numRowsTT; ++i)
        for (unsigned j = 0; j < numInputs; ++j)
            inputTT[i][j] = (i >> j) & 1u;
}

/* -------------------------------------------------------
 * loadTT
 * Reads circuit dimensions and expected outputs from file.
 * ------------------------------------------------------- */
void CircuitData::loadTT(std::FILE* pf)
{
    numRowsTT = 1u;
    numInputs = static_cast<unsigned>(readNumber(pf));
    numRowsTT <<= numInputs;
    numOutputs      = static_cast<unsigned>(readNumber(pf));
    numTotalOutputs = numOutputs * numRowsTT;

    generateTT();

    outputTT.assign(numRowsTT, std::vector<unsigned>(numOutputs, 0u));

    for (unsigned i = 0; i < numRowsTT; ++i)
        for (unsigned j = 0; j < numOutputs; ++j)
            outputTT[i][j] = static_cast<unsigned>(readNumber(pf));

    packTT();
}

/* -------------------------------------------------------
 * packTT
 * Builds the bit-packed truth-table representation used
 * by the bit-parallel MatrixDecoder::evaluate.
 *
 *   - inputTTPacked[in][w]  : bit 'b' of word 'w' is the
 *                             value of input 'in' on row
 *                             (w*64 + b).
 *   - outputTTTarget[out][w]: expected output bit.
 *   - outputTTCare[out][w]  : 1 iff the row is NOT a
 *                             DONTCARE.  Bits beyond
 *                             numRowsTT are left zero so
 *                             they never contribute to
 *                             the match count.
 * ------------------------------------------------------- */
void CircuitData::packTT()
{
    numWords = (numRowsTT + 63u) / 64u;

    inputTTPacked.assign(numInputs,
                         std::vector<uint64_t>(numWords, 0ULL));
    outputTTTarget.assign(numOutputs,
                          std::vector<uint64_t>(numWords, 0ULL));
    outputTTCare.assign(numOutputs,
                        std::vector<uint64_t>(numWords, 0ULL));

    for (unsigned row = 0; row < numRowsTT; ++row) {
        const unsigned w = row >> 6;            /* row / 64 */
        const uint64_t b = 1ULL << (row & 63u); /* row % 64 */

        for (unsigned in = 0; in < numInputs; ++in)
            if (inputTT[row][in]) inputTTPacked[in][w] |= b;

        for (unsigned out = 0; out < numOutputs; ++out) {
            const unsigned v = outputTT[row][out];
            if (v != DONTCARE) {
                outputTTCare[out][w] |= b;
                if (v) outputTTTarget[out][w] |= b;
            }
        }
    }
}

/* -------------------------------------------------------
 * printTT
 * Displays the full truth table on stdout.
 * ------------------------------------------------------- */
void CircuitData::printTT() const
{
    std::printf("      Number of inputs in the truth table: %u\n",  numInputs);
    std::printf("     Number of outputs in the truth table: %u\n",  numOutputs);
    std::printf("       Number of rows in the truth table: %u\n\n", numRowsTT);

    std::printf("Truth Table\n\n");
    for (unsigned i = 0; i < numInputs;  ++i) std::printf("E%u ", i);
    for (unsigned i = 0; i < numOutputs; ++i) std::printf("S%u ", i);

    for (unsigned i = 0; i < numRowsTT; ++i) {
        std::printf("\n");
        for (int j = static_cast<int>(numInputs) - 1; j >= 0; --j)
            std::printf(" %u ", inputTT[i][static_cast<unsigned>(j)]);
        for (unsigned j = 0; j < numOutputs; ++j)
            std::printf(outputTT[i][j] == DONTCARE ? " * " : " %u ", outputTT[i][j]);
    }
    std::printf("\n\n");
}

/* -------------------------------------------------------
 * freeMemory
 * Vectors free themselves when Swarm is destroyed, but
 * this allows explicit early release if desired.
 * ------------------------------------------------------- */
void CircuitData::freeMemory()
{
    inputTT.clear();
    outputTT.clear();
    inputTTPacked.clear();
    outputTTTarget.clear();
    outputTTCare.clear();
}

/* -------------------------------------------------------
 * readNumber (static)
 * Parses the next numeric token from a PSO data file,
 * skipping comments (lines starting with ';').
 * ------------------------------------------------------- */
double CircuitData::readNumber(std::FILE* pf)
{
    char car = '\0';
    std::string num;

    /* Skip non-numeric characters; handle comment lines starting with ';' */
    while (!std::feof(pf)) {
        car = static_cast<char>(std::fgetc(pf));
        if (isDigit(car) || car == '.') break;
        if (car == ';') {
            while (car != EndLine && !std::feof(pf))
                car = static_cast<char>(std::fgetc(pf));
            car = static_cast<char>(std::fgetc(pf));
            if (isDigit(car) || car == '.') break;
        }
    }

    /* Consume the numeric token */
    while (!std::feof(pf) && car != EndLine && car != ' ' && car != ';') {
        if (isDigit(car) || car == '.') num += car;
        car = static_cast<char>(std::fgetc(pf));
    }

    /* Skip rest of inline comment if present */
    if (car == ';')
        while ((car = static_cast<char>(std::fgetc(pf))) != EndLine && !std::feof(pf));

    return std::atof(num.c_str());
}

/* -------------------------------------------------------
 * readString (static)
 * Parses the next string token from a PSO data file.
 * Returns std::string instead of writing into a char*.
 * ------------------------------------------------------- */
std::string CircuitData::readString(std::FILE* pf)
{
    char car = '\0';
    std::string str;

    /* Skip until the start of an alphanumeric token */
    while (!std::feof(pf)) {
        car = static_cast<char>(std::fgetc(pf));
        if (car == ';') {
            while (car != EndLine && !std::feof(pf))
                car = static_cast<char>(std::fgetc(pf));
            car = static_cast<char>(std::fgetc(pf));
        }
        if (std::isalpha(static_cast<unsigned char>(car))) break;
    }

    /* Consume the token */
    while (!std::feof(pf) && car != EndLine && car != ' ' && car != ';' && car != '\t') {
        str += car;
        car = static_cast<char>(std::fgetc(pf));
    }

    if (car == ';')
        while ((car = static_cast<char>(std::fgetc(pf))) != EndLine && !std::feof(pf));

    return str;
}

/* -------------------------------------------------------
 * isDigit (static)
 * Returns true if c is an ASCII decimal digit.
 * ------------------------------------------------------- */
bool CircuitData::isDigit(char c)
{
    return c >= '0' && c <= '9';
}

} // namespace PSwarm
