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
/* File: statistics.cpp                                 */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Library to generate statistics per      */
/* generation and per run.  Statistics include the best */
/* and worst particle, mean fitness, standard deviation,*/
/* variance, and the Boolean expression of the best     */
/* particle found in each generation.                   */
/********************************************************/

#include "statistics.h"
#include "psomatrixcircuit.h"   /* Full Swarm definition */
#include <cstdio>
#include <cmath>
#include <string>

namespace PSwarm {

/* -------------------------------------------------------
 * Helpers
 * ------------------------------------------------------- */

/* Returns a string representation for the representation type */
static const char* reprName(unsigned rep) {
    if (rep == BINARY)    return "Binary";
    if (rep == INTEGER_A) return "Integer A";
    return "Integer B";
}

/* -------------------------------------------------------
 * initStatistics
 * Zero-initialises one statistics block.
 * ------------------------------------------------------- */
void initStatistics(Statistics& st)
{
    st.squaredFitness = 0.0;
    st.meanFitness    = 0.0;
}

/* -------------------------------------------------------
 * dtoa
 * Converts an unsigned integer to a string.
 * Replaces the original fcvt-based dtoa with std::to_string.
 * ------------------------------------------------------- */
std::string dtoa(unsigned num)
{
    return std::to_string(num);
}

/* -------------------------------------------------------
 * chromToString
 * Converts a particle's chromosome into a printable string.
 * Returns std::string instead of writing into a char buffer.
 * ------------------------------------------------------- */
std::string chromToString(const Particle& par, const Swarm& swarm)
{
    std::string str = "'";
    for (unsigned i = 0; i < swarm.nAllele; ++i) {
        str += dtoa(par.chromX[i]);
        if (swarm.representation != BINARY) str += ' ';
    }
    return str;
}

/* -------------------------------------------------------
 * globalHeader
 * Writes the CSV header shared across all runs.
 * All Swarm configuration data is passed via const Swarm&
 * instead of being read from global variables.
 * ------------------------------------------------------- */
void globalHeader(const std::string& filename, const Swarm& swarm)
{
    std::FILE* pf = std::fopen(filename.c_str(), "w");
    if (!pf) return;

    std::fprintf(pf, "Centro de Investigación y de Estudios Avanzados del IPN\n");
    std::fprintf(pf, "Design of Combinational Logic Circuits using Particle Swarm Optimization\n");
    std::fprintf(pf, "Department of Electrical Engineering - Computing Section\n");
    std::fprintf(pf, "Erika Hernández Luna  eluna@computacion.cs.cinvestav.mx\n\n");
    std::fprintf(pf, "Statistics of the results obtained in all program runs\n\n");

    std::fprintf(pf, "     Processed File:,   %s\n", swarm.nfInput.c_str());
    std::fprintf(pf, "   Population Size:,   %u\n",  swarm.tPop);
    std::fprintf(pf, "Number of Generations:,   %u\n", swarm.nGen);
    std::fprintf(pf, "                 Phi 1:, %.3f\n", swarm.phi1);
    std::fprintf(pf, "                 Phi 2:, %.3f\n", swarm.phi2);
    std::fprintf(pf, "                  Vmax:, %.3f\n", swarm.vMax);
    std::fprintf(pf, "Mutation Percentage:, %.3f\n",    swarm.pMut * 100.0);
    std::fprintf(pf, "  Number of Gates:,   %u\n",     swarm.circuit.numGates);
    std::fprintf(pf, "        Cardinality:,   %u\n",   swarm.cardinality);
    std::fprintf(pf, "Rows in the Matrix:,   %u\n",    swarm.decoder.numRows);
    std::fprintf(pf, "Columns in the Matrix:,   %u\n", swarm.decoder.numCols);
    std::fprintf(pf, "       Neighborhood:,   %u\n",   swarm.tNeigh);
    std::fprintf(pf, "   Representation:, %s\n",       reprName(swarm.representation));

    std::fprintf(pf, "\n");
    std::fprintf(pf, "Run,Mean Fitness,Generation,");
    std::fprintf(pf, "Best Chromosome,Fitness,");
    std::fprintf(pf, "Worst Chromosome,Fitness,");
    std::fprintf(pf, "Standard Deviation,Variance, Expression\n");

    std::fclose(pf);
}

/* -------------------------------------------------------
 * runStatistics
 * Appends one run's summary row to the global CSV.
 * ------------------------------------------------------- */
void runStatistics(const std::string& filename, unsigned run, Swarm& swarm)
{
    std::FILE* pf = std::fopen(filename.c_str(), "a+");
    if (!pf) return;

    double variance = swarm.Run.squaredFitness
                    - swarm.Run.meanFitness * swarm.Run.meanFitness;

    std::string best  = chromToString(swarm.Run.best,  swarm);
    std::string worst = chromToString(swarm.Run.worst, swarm);

    std::fprintf(pf, "%u,%.8f,%u,%s,%.8f,%s,%.8f,%.8f,%.12f",
        run,
        swarm.Run.meanFitness,
        swarm.Run.generation,
        best.c_str(),
        swarm.Run.best.fitness,
        worst.c_str(),
        swarm.Run.worst.fitness,
        std::sqrt(variance),
        variance);

    for (unsigned i = 0; i < swarm.circuit.numOutputs; ++i) {
        std::string expr = swarm.decoder.expression(
            swarm.Run.best.chromX, i, swarm, swarm.circuit);
        std::fprintf(pf, ",%s", expr.c_str());
    }
    std::fprintf(pf, "\n");
    std::fclose(pf);
}

/* -------------------------------------------------------
 * runHeader
 * Writes the per-run CSV header.
 * ------------------------------------------------------- */
void runHeader(const std::string& filename, unsigned seed, const Swarm& swarm)
{
    std::FILE* pf = std::fopen(filename.c_str(), "w");
    if (!pf) return;

    std::fprintf(pf, "Centro de Investigación y de Estudios Avanzados del IPN\n");
    std::fprintf(pf, "Design of Combinational Logic Circuits using Particle Swarm Optimization\n");
    std::fprintf(pf, "Department of Electrical Engineering - Computing Section\n");
    std::fprintf(pf, "Erika Hernández Luna  eluna@computacion.cs.cinvestav.mx\n\n");
    std::fprintf(pf, "Results obtained in a program run\n\n");

    std::fprintf(pf, "     Processed File:,   %s\n", swarm.nfInput.c_str());
    std::fprintf(pf, "   Population Size:,   %u\n",  swarm.tPop);
    std::fprintf(pf, "Number of Generations:,   %u\n", swarm.nGen);
    std::fprintf(pf, "                 Phi 1:, %.3f\n", swarm.phi1);
    std::fprintf(pf, "                 Phi 2:, %.3f\n", swarm.phi2);
    std::fprintf(pf, "                  Vmax:, %.3f\n", swarm.vMax);
    std::fprintf(pf, "Mutation Percentage:, %.3f\n",    swarm.pMut * 100.0);
    std::fprintf(pf, "  Number of Gates:,   %u\n",     swarm.circuit.numGates);
    std::fprintf(pf, "        Cardinality:,   %u\n",   swarm.cardinality);
    std::fprintf(pf, "Rows in the Matrix:,   %u\n",    swarm.decoder.numRows);
    std::fprintf(pf, "Columns in the Matrix:,   %u\n", swarm.decoder.numCols);
    std::fprintf(pf, "       Neighborhood:,   %u\n",   swarm.tNeigh);
    std::fprintf(pf, "             Seed:,   %u\n",     seed);
    std::fprintf(pf, "   Representation:, %s\n",       reprName(swarm.representation));

    std::fprintf(pf, "\n");
    std::fprintf(pf, "Generation,Mean Fitness,Best Fitness,Worst Fitness,Best Chromosome\n");
    std::fclose(pf);
}

/* -------------------------------------------------------
 * generationStats
 * Appends one generation row to the per-run CSV.
 * ------------------------------------------------------- */
void generationStats(const std::string& filename, unsigned gen, const Swarm& swarm)
{
    std::FILE* pf = std::fopen(filename.c_str(), "a+");
    if (!pf) return;

    std::string best = chromToString(swarm.Gen.best, swarm);
    std::fprintf(pf, "%u,%.8f,%.8f,%.8f,%s,\n",
        gen,
        swarm.Gen.meanFitness,
        swarm.Gen.best.fitness,
        swarm.Gen.worst.fitness,
        best.c_str());

    std::fclose(pf);
}

/* -------------------------------------------------------
 * runFooter
 * Appends the final run summary (best individual, Boolean
 * expressions) to the per-run CSV.
 * ------------------------------------------------------- */
void runFooter(const std::string& filename, Swarm& swarm)
{
    std::FILE* pf = std::fopen(filename.c_str(), "a+");
    if (!pf) return;

    std::string best = chromToString(swarm.Run.best, swarm);

    std::fprintf(pf, "\n");
    std::fprintf(pf, "Mean Fitness:, %.6f\n\n", swarm.Run.meanFitness);
    std::fprintf(pf, "   Best Individual\n\n");
    std::fprintf(pf, "          Fitness:, %.8f\n", swarm.Run.best.fitness);
    std::fprintf(pf, "       Violations:,   %d\n",
        static_cast<int>(swarm.circuit.numTotalOutputs)
        - static_cast<int>(swarm.Run.best.numEqual));
    std::fprintf(pf, "           Gates:,   %u\n",  swarm.Run.best.numGates);
    std::fprintf(pf, "       Chromosome:,   %s\n", best.c_str());
    std::fprintf(pf, "   Boolean Expression:\n");

    for (unsigned i = 0; i < swarm.circuit.numOutputs; ++i) {
        std::string expr = swarm.decoder.expression(
            swarm.Run.best.chromX, i, swarm, swarm.circuit);
        std::fprintf(pf, "%s\n", expr.c_str());
    }

    std::fclose(pf);
}

} // namespace PSwarm
