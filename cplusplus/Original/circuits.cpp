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

namespace PSwarm {

/* -------------------------------------------------------
 * generateTT
 * Builds the input combinations for the truth table.
 * Replaces the original malloc-based approach with nested
 * vectors that manage their own memory.
 * ------------------------------------------------------- */
// void CircuitData::generateTT()
// {
//     inputTT.assign(numRowsTT, std::vector<unsigned>(numInputs, 0u));

//     for (unsigned i = 0; i < numRowsTT; ++i)
//         for (unsigned j = 0; j < numInputs; ++j)
//             inputTT[i][j] = (i >> j) & 1u;
// }

/* -------------------------------------------------------
 * loadTT
 * Reads circuit dimensions and expected outputs from file.
 * ------------------------------------------------------- */
// void CircuitData::loadTT(std::FILE* pf)
// {
//     numRowsTT = 1u;
//     numInputs = static_cast<unsigned>(readNumber(pf));
//     numRowsTT <<= numInputs;
//     numOutputs      = static_cast<unsigned>(readNumber(pf));
//     numTotalOutputs = numOutputs * numRowsTT;

//     generateTT();

//     outputTT.assign(numRowsTT, std::vector<unsigned>(numOutputs, 0u));

//     for (unsigned i = 0; i < numRowsTT; ++i)
//         for (unsigned j = 0; j < numOutputs; ++j)
//             outputTT[i][j] = static_cast<unsigned>(readNumber(pf));
// }

/* -------------------------------------------------------
 * loadTT
 * Reads a .plu-format truth-table block from an
 * already-open file stream:
 * 
 *     .i <num_inputs>
 *     .o <num_outputs>
 *     .p <num_instances>
 *     <numInputs packed words> <numOutputs packed words> x numInstances
 *     .e
 * 
 * Each packed word's bit k encodes one logical truth-table
 * row. Data is stored directly in packed form (no
 * expansion to one-bit-per-row), so MatrixDecoder can
 * evaluate 32 rows per gate operation.
 * ------------------------------------------------------- */
void CircuitData::loadTT(std::FILE* pf)
{
    readToken(pf);
    numInputs = readUnsigned(pf);   // ".i"

    readToken(pf);
    numOutputs = readUnsigned(pf);  // ".o"

    readToken(pf);
    numInstances = readUnsigned(pf); // ".p"

    /* Total logical rows. Assumes numInputs <= 31, true for any
       realistic benchmark circuit (2^31 rows would already be
       infeasible to fully enumerate). */
    numRowsTT = 1u << numInputs;
    numTotalOutputs = numOutputs * numRowsTT;

    inputTT.assign(numInstances, std::vector<unsigned>(numInputs, 0u));
    outputTT.assign(numInstances, std::vector<unsigned>(numOutputs, 0u));
    chunkValidBits.assign(numInstances, PLU_CHUNK_BITS);

    const unsigned remainder = numRowsTT % PLU_CHUNK_BITS;
    if (remainder != 0 && numInstances > 0) {
        chunkValidBits[numInstances - 1] = remainder;
    }

    for (unsigned i = 0; i < numInstances; ++i) {
        for (unsigned j = 0; j < numInputs; ++j) {
            inputTT[i][j] = readUnsigned(pf);
        }

        for (unsigned j = 0; j < numOutputs; ++j) {
            outputTT[i][j] = readUnsigned(pf);
        }
    }

    readToken(pf); // ".e"
}

// /* -------------------------------------------------------
//  * printTT
//  * Displays the full truth table on stdout.
//  * ------------------------------------------------------- */
// void CircuitData::printTT() const
// {
//     std::printf("      Number of inputs in the truth table: %u\n",  numInputs);
//     std::printf("     Number of outputs in the truth table: %u\n",  numOutputs);
//     std::printf("       Number of rows in the truth table: %u\n\n", numRowsTT);

//     std::printf("Truth Table\n\n");
//     for (unsigned i = 0; i < numInputs;  ++i) std::printf("E%u ", i);
//     for (unsigned i = 0; i < numOutputs; ++i) std::printf("S%u ", i);

//     for (unsigned i = 0; i < numRowsTT; ++i) {
//         std::printf("\n");
//         for (int j = static_cast<int>(numInputs) - 1; j >= 0; --j)
//             std::printf(" %u ", inputTT[i][static_cast<unsigned>(j)]);
//         for (unsigned j = 0; j < numOutputs; ++j)
//             std::printf(outputTT[i][j] == DONTCARE ? " * " : " %u ", outputTT[i][j]);
//     }
//     std::printf("\n\n");
// }

/* -------------------------------------------------------
 * printTT
 * Displays the packed truth table on stdout. Values are
 * shown in hex since each one packs 32 individual rows.
 * ------------------------------------------------------- */
void CircuitData::printTT() const
{
    std::printf("      Number of inputs in the truth table: %u\n",  numInputs);
    std::printf("     Number of outputs in the truth table: %u\n",  numOutputs);
    std::printf("    Number of logical rows (2^numInputs): %u\n", numRowsTT);
    std::printf("   Number of packed instances (.p, 32 rows each): %u\n\n", numRowsTT);

    std::printf("Packed Truth Table (each value packs up to 32 rows)\n\n");

    for (unsigned i = 0; i < numInstances; ++i) {
        std::printf("Instance %2u: ", i);
        for (unsigned j = 0; j < numInputs; ++j) {
            std::printf("I%u=0x%08x ", j, inputTT[i][j]);
        }
        std::printf(" | ");
        for (unsigned j = 0; j < numOutputs; ++j) {
            std::printf("O%u=0x%08x ", j, outputTT[i][j]);
        }
        std::printf("   (valid bits: %u)\n", chunkValidBits[i]);
    }
    std::printf("\n");
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
    chunkValidBits.clear();
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
 * readToken (static, new)
 * Reads the next whitespace-delimited token verbatim --
 * including a leading '.', so ".i", ".o", ".p", ".e" come
 * back intact (unlike readNumber(), which treats a leading
 * '.' as the start of a decimal number and would mis-parse
 * these keywords). Also used for the plain numeric tokens
 * inside the .plu block, via readUnsigned() below.
 * ------------------------------------------------------- */
 std::string CircuitData::readToken(std::FILE* pf)
{
    char car = '\0';
    std::string str;
 
    /* Skip whitespace and comment lines */
    while (!std::feof(pf)) {
        car = static_cast<char>(std::fgetc(pf));
        if (car == ';') {
            while (car != EndLine && !std::feof(pf))
                car = static_cast<char>(std::fgetc(pf));
            continue;
        }
        if (car != ' ' && car != '\t' && car != EndLine) break;
    }
 
    while (!std::feof(pf) && car != ' ' && car != '\t' && car != EndLine && car != ';') {
        str += car;
        car = static_cast<char>(std::fgetc(pf));
    }
 
    if (car == ';')
        while ((car = static_cast<char>(std::fgetc(pf))) != EndLine && !std::feof(pf));
 
    return str;
}

/* -------------------------------------------------------
 * readUnsigned (static, new)
 * Reads the next token and parses it as an unsigned value
 * via strtoul rather than atoi/atof, since .plu packed
 * words can reach 4294967295 (0xFFFFFFFF) -- well past
 * what atoi (signed int) can hold without overflow.
 * ------------------------------------------------------- */
unsigned CircuitData::readUnsigned(std::FILE* pf)
{
    std::string tok = readToken(pf);
    return static_cast<unsigned>(std::strtoul(tok.c_str(), nullptr, 10));
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
