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
/* File: psomatrixcircuit.cpp                           */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Main program that applies Particle      */
/* Swarm Optimization to the design of logic circuits.  */
/*                                                      */
/* EAPSO modification (April 2026):                     */
/*   Added computeEntropyRanks() and computeQuality()   */
/*   implementing the entropy-based diversity mechanism  */
/*   from Guo et al. (Mathematics 2019, 7, 414).        */
/*   Modified evaluatePopulation() to call these after  */
/*   fitness evaluation and to select neighbourhood     */
/*   bests by quality Q rather than raw fitness.        */
/*   All other methods are unchanged.                   */
/********************************************************/

#include "psomatrixcircuit.h"
#include "random.h"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <string>

/* -------------------------------------------------------
 * Entry point
 * ------------------------------------------------------- */
int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::printf("Usage: %s <input_file>\n", argv[0]);
        return 0;
    }

    PSwarm::Swarm swarm;

    if (swarm.loadParameters(argv[1])) {
        swarm.initVariables();
        swarm.reserveMemory();

        globalHeader(swarm.nfRun, swarm);

        for (unsigned run = 0; run < swarm.nRun; ++run) {
            std::printf("\n\nRun %02u Started\n", run);
            initStatistics(swarm.Run);
            swarm.pSwarm(run);
            runStatistics(swarm.nfRun, run, swarm);
            std::printf("\n\nRun %02u Finished\n", run);
        }

        swarm.freeMemory();
        return 1;
    }

    return 0;
}

/* -------------------------------------------------------
 * Swarm member implementations
 * ------------------------------------------------------- */

namespace PSwarm {

/* Load all parameters and the truth table from file */
bool Swarm::loadParameters(const std::string& filename)
{
    std::FILE* input = std::fopen(filename.c_str(), "r");
    if (!input) {
        std::printf("\nError: could not load data file '%s'\n", filename.c_str());
        return false;
    }

    std::printf("Currently loading parameters!\n");

    nfInput        = filename;
    circuit.loadTT(input);

    tPop           = static_cast<unsigned>(CircuitData::readNumber(input));
    nGen           = static_cast<unsigned>(CircuitData::readNumber(input));
    nRun           = static_cast<unsigned>(CircuitData::readNumber(input));
    tNeigh         = static_cast<unsigned>(CircuitData::readNumber(input));
    phi1           = CircuitData::readNumber(input);
    phi2           = CircuitData::readNumber(input);
    vMax           = CircuitData::readNumber(input);
    pMut           = CircuitData::readNumber(input) / 100.0;
    representation = static_cast<unsigned>(CircuitData::readNumber(input));
    circuit.numGates = static_cast<unsigned>(CircuitData::readNumber(input));
    cardinality    = static_cast<unsigned>(CircuitData::readNumber(input));
    decoder.numRows = static_cast<unsigned>(CircuitData::readNumber(input));
    decoder.numCols = static_cast<unsigned>(CircuitData::readNumber(input));

    nfGen = CircuitData::readString(input);
    nfRun = CircuitData::readString(input) + ".csv";

    std::fclose(input);

    circuit.printTT();
    printParameters();
    return true;
}

/* Print all loaded parameters to stdout */
void Swarm::printParameters() const
{
    std::printf("\n           Population size: %u",   tPop);
    std::printf("\n         Number of generations: %u", nGen);
    std::printf("\n                 Number of runs: %u", nRun);
    std::printf("\n             Neighbourhood size: %u",  tNeigh);
    std::printf("\n              PSO phi1 parameter: %.2f", phi1);
    std::printf("\n              PSO phi2 parameter: %.2f", phi2);
    std::printf("\n                 PSO vmax: %.2f",   vMax);
    std::printf("\n           Mutation percentage: %.2f", pMut * 100.0);
    std::printf("\n                 Representation: %s",
        representation == BINARY    ? "Binary" :
        representation == INTEGER_A ? "Integer A" : "Integer B");
    std::printf("\n      Number of available gates: %u", circuit.numGates);
    std::printf("\n                    Cardinality: %u", cardinality);
    std::printf("\n          Rows in the matrix: %u",   decoder.numRows);
    std::printf("\n       Columns in the matrix: %u",   decoder.numCols);
    std::printf("\n   Per-run results file: %sX.csv",   nfGen.c_str());
    std::printf("\nGlobal results file: %s\n",          nfRun.c_str());

    /* EAPSO parameters */
    std::printf("\n--- EAPSO Diversity Parameters ---");
    std::printf("\n       Fitness segments (m): %u",   EAPSO_NUM_SEGMENTS);
    std::printf("\n         Fitness weight (alpha): %.2f\n\n", EAPSO_ALPHA);
}

/* Initialise derived variables */
void Swarm::initVariables()
{
    if (tNeigh > tPop) tNeigh = tPop;

    decoder.tMat = decoder.numCols * decoder.numRows;

    initBounds();
}

/* Initialise variable bounds and compute nAllele */
void Swarm::initBounds()
{
    nVar = decoder.numRows * decoder.numCols * 3; /* input1, input2, gateType per cell */

    lInf.assign(nVar, 0u);
    lSup.assign(nVar, 0u);

    for (unsigned i = 0; i < nVar; ++i) {
        lInf[i] = 0;
        lSup[i] = (i + 1) % 3 ? decoder.numRows - 1 : circuit.numGates - 1;
    }

    if (representation == BINARY) {
        bitVariable.resize(nVar);
        nAllele = 0;
        for (unsigned i = 0; i < nVar; ++i) {
            bitVariable[i] = static_cast<unsigned>(
                std::ceil(std::log(std::fabs(static_cast<double>(lInf[i]))
                                 + std::fabs(static_cast<double>(lSup[i])) + 1.0)
                         / std::log(2.0)));
            nAllele += bitVariable[i];
        }
    } else {
        nAllele = nVar;
    }
}

/* Allocate all population and working memory */
void Swarm::reserveMemory()
{
    population.assign(tPop, Particle(nAllele));
    bestSocialExp.assign(tPop, Particle(nAllele));
    bestIndividualExp.assign(tPop, Particle(nAllele));

    reserveParticleMemory(Gen.best,  true);
    reserveParticleMemory(Gen.worst, true);
    reserveParticleMemory(Run.best,  true);
    reserveParticleMemory(Run.worst, true);

    decoder.reserve(circuit);
}

/* Initialise vi and chromX vectors of a statistics particle */
void Swarm::reserveParticleMemory(Particle& par, bool /*withChrom*/)
{
    par.vi.assign(nAllele, 0.0);
    par.chromX.assign(nAllele, 0u);
}

/* Release all memory */
void Swarm::freeMemory()
{
    circuit.freeMemory();
    decoder.free();

    population.clear();
    bestSocialExp.clear();
    bestIndividualExp.clear();

    lInf.clear();
    lSup.clear();
    bitVariable.clear();
}

/* Execute one complete PSO run */
void Swarm::pSwarm(unsigned run)
{
    unsigned seed = initRandom(0);
    std::string fileGen;
    runFileName(run, fileGen);

    initPopulation();
    runHeader(fileGen, seed, *this);

    for (unsigned gen = 0; gen < nGen; ++gen) {
        initStatistics(Gen);
        evaluatePopulation(gen);
        runInfo(fileGen, gen);
        PSOAlgorithm(gen);
        mutation();
    }

    runFooter(fileGen, *this);
}

/* Build per-run filename */
void Swarm::runFileName(unsigned run, std::string& outName) const
{
    outName = nfGen + (run == 0 ? "0" : std::to_string(run)) + ".csv";
}

/* Generate the initial random population */
void Swarm::initPopulation()
{
    for (unsigned i = 0; i < tPop; ++i)
        for (unsigned j = 0; j < nAllele; ++j)
            population[i].chromX[j] = (representation == BINARY)
                ? flip(0.5)
                : rndIR(lInf[j], (lSup[j] + 1) * cardinality) % (lSup[j] + 1);
}

/* -------------------------------------------------------
 * computeEntropyRanks  [EAPSO addition]
 *
 * Implements the entropy diversity measurement from
 * Guo et al. (2019), Section 3, Equations (3)-(4).
 *
 * Steps:
 *   1. Find the fitness range [fMin, fMax] across all
 *      bestIndividualExp particles.
 *   2. Divide that range into EAPSO_NUM_SEGMENTS equal
 *      bins and count how many particles fall in each.
 *   3. Assign fitnessRank: sort particles by descending
 *      fitness; rank 1 = best (highest) fitness.
 *   4. Assign entropyRank: sort segments by ascending
 *      occupancy count; the segment with fewest particles
 *      gets rank 1 (most novel / most sparse region).
 *      Particles in the same segment share a rank.
 *
 * Both ranks are written into bestIndividualExp[i] so
 * that computeQuality() can combine them immediately
 * after this call.
 * ------------------------------------------------------- */
void Swarm::computeEntropyRanks()
{
    /* ---- 1. Fitness range ---- */
    double fMin = bestIndividualExp[0].fitness;
    double fMax = bestIndividualExp[0].fitness;
    for (unsigned i = 1; i < tPop; ++i) {
        fMin = std::min(fMin, bestIndividualExp[i].fitness);
        fMax = std::max(fMax, bestIndividualExp[i].fitness);
    }

    /* ---- 2. Segment occupancy ---- */
    /* If all particles have identical fitness the range is
     * zero; assign every particle to segment 0 to avoid
     * a division-by-zero.  All ranks will be equal, so Q
     * degenerates to raw fitness rank — correct behaviour. */
    const double range = fMax - fMin;
    const unsigned m   = EAPSO_NUM_SEGMENTS;

    std::vector<unsigned> segCount(m, 0u);
    std::vector<unsigned> particleSeg(tPop, 0u);

    for (unsigned i = 0; i < tPop; ++i) {
        unsigned seg = 0u;
        if (range > 0.0) {
            /* Map fitness linearly into [0, m-1].
             * Clamp the top value so fMax lands in the
             * last segment rather than one past it. */
            double normalised = (bestIndividualExp[i].fitness - fMin) / range;
            seg = static_cast<unsigned>(normalised * static_cast<double>(m));
            if (seg >= m) seg = m - 1u;
        }
        particleSeg[i] = seg;
        ++segCount[seg];
    }

    /* ---- 3. fitnessRank ----
     * Build an index array sorted by descending fitness,
     * then assign rank 1 to the particle at position 0. */
    std::vector<unsigned> fitnessOrder(tPop);
    std::iota(fitnessOrder.begin(), fitnessOrder.end(), 0u);
    std::sort(fitnessOrder.begin(), fitnessOrder.end(),
        [this](unsigned a, unsigned b) {
            return bestIndividualExp[a].fitness > bestIndividualExp[b].fitness;
        });

    for (unsigned rank = 0; rank < tPop; ++rank)
        bestIndividualExp[fitnessOrder[rank]].fitnessRank = rank + 1u;

    /* ---- 4. entropyRank ----
     * Build a segment index array sorted by ascending
     * occupancy (fewest particles = best diversity rank).
     * All particles in the same segment receive the same
     * entropyRank equal to that segment's rank. */
    std::vector<unsigned> segOrder(m);
    std::iota(segOrder.begin(), segOrder.end(), 0u);
    std::sort(segOrder.begin(), segOrder.end(),
        [&segCount](unsigned a, unsigned b) {
            return segCount[a] < segCount[b];
        });

    /* Map segment index → its entropy rank */
    std::vector<unsigned> segRank(m, 0u);
    for (unsigned rank = 0; rank < m; ++rank)
        segRank[segOrder[rank]] = rank + 1u;

    for (unsigned i = 0; i < tPop; ++i)
        bestIndividualExp[i].entropyRank = segRank[particleSeg[i]];
}

/* -------------------------------------------------------
 * computeQuality  [EAPSO addition]
 *
 * Implements Equation (5) from Guo et al. (2019):
 *   Q_i = alpha * fitnessRank_i + entropyRank_i
 *
 * beta is fixed at 1 (as in the paper) so only alpha
 * needs tuning.  A lower Q is better — the particle is
 * both high-quality (low fitnessRank) and occupies a
 * sparsely populated region of the fitness landscape
 * (low entropyRank).
 *
 * Must be called after computeEntropyRanks().
 * Results are written into bestIndividualExp[i].quality.
 * ------------------------------------------------------- */
void Swarm::computeQuality()
{
    for (unsigned i = 0; i < tPop; ++i) {
        bestIndividualExp[i].quality =
            EAPSO_ALPHA * static_cast<double>(bestIndividualExp[i].fitnessRank)
            +             static_cast<double>(bestIndividualExp[i].entropyRank);
    }
}

/* -------------------------------------------------------
 * evaluatePopulation
 *
 * EAPSO modification: after Step 1 (fitness evaluation
 * and personal best update), call computeEntropyRanks()
 * and computeQuality() to assign Q scores.
 *
 * In Step 2 (neighbourhood scan), the comparison that
 * selects the neighbourhood best now uses .quality
 * (lower is better) instead of .fitness (higher is
 * better).  This causes the social exemplar to be the
 * particle that is both competitively fit AND occupies a
 * sparse region of the fitness landscape, which is the
 * core anti-premature-convergence mechanism of EAPSO.
 *
 * All downstream reporting (Gen.best, Gen.worst,
 * Run.best, Run.worst, meanFitness, squaredFitness)
 * continues to use raw .fitness so that statistics are
 * comparable to unmodified runs.
 * ------------------------------------------------------- */
void Swarm::evaluatePopulation(unsigned gen)
{
    /* Step 1 – evaluate fitness and update personal bests */
    for (unsigned i = 0; i < tPop; ++i) {
        evaluateParticle(population[i]);

        if (gen == 0 || population[i].fitness > bestIndividualExp[i].fitness)
            copyParticle(population[i], bestIndividualExp[i]);
    }

    /* EAPSO Step 1b – compute entropy ranks and quality scores
     * for all personal bests before neighbourhood selection.
     * This is the frequency-domain diversity measurement from
     * Section 3 of Guo et al. (2019), Equations (3)-(5). */
    computeEntropyRanks();
    computeQuality();

    /* Step 2 – find neighbourhood bests using Q score,
     * update social bests and generation statistics.
     *
     * The ring scan is identical to the original except
     * the comparison flips from "higher fitness is better"
     * to "lower quality is better", because Q encodes
     * both fitness rank and entropy rank in a single
     * minimisation objective. */
    for (unsigned i = 0; i < tPop; ++i) {
        int p  = static_cast<int>(i);
        int p0 = p;
        int p1 = p;

        for (unsigned j = 0; j < tNeigh; j += 2) {
            p0 = (p0 + 1 < static_cast<int>(tPop)) ? p0 + 1 : 0;
            p1 = (p1 - 1 >= 0) ? p1 - 1 : static_cast<int>(tPop) - 1;

            /* EAPSO: compare by quality (lower = better)
             * instead of fitness (higher = better) */
            if (bestIndividualExp[p0].quality < bestIndividualExp[p].quality) p = p0;
            if (bestIndividualExp[p1].quality < bestIndividualExp[p].quality) p = p1;
        }

        copyParticle(bestIndividualExp[p], bestSocialExp[i]);

        /* Statistics accumulation uses raw fitness, unchanged */
        Gen.meanFitness    += bestSocialExp[i].fitness / tPop;
        Gen.squaredFitness += std::pow(bestSocialExp[i].fitness, 2.0) / tPop;

        if (i == 0) {
            copyParticle(bestSocialExp[i], Gen.best);
            copyParticle(bestSocialExp[i], Gen.worst);
        } else {
            if (bestSocialExp[i].fitness > Gen.best.fitness)
                copyParticle(bestSocialExp[i], Gen.best);
            if (bestSocialExp[i].fitness < Gen.worst.fitness)
                copyParticle(bestSocialExp[i], Gen.worst);
        }
    }
}

/* Compute fitness of a single particle */
void Swarm::evaluateParticle(Particle& par)
{
    decoder.evaluate(par.chromX, par.numEqual, *this, circuit);

    par.numNoGates = 0;
    if (par.numEqual >= circuit.numTotalOutputs) {
        par.numGates   = decoder.countGates(circuit);
        par.numNoGates = decoder.tMat - par.numGates;
    }

    par.fitness = static_cast<double>(par.numEqual + par.numNoGates);
}

/* Deep-copy one particle into another */
void Swarm::copyParticle(const Particle& src, Particle& dst)
{
    dst.numGates    = src.numGates;
    dst.numEqual    = src.numEqual;
    dst.numNoGates  = src.numNoGates;
    dst.fitness     = src.fitness;
    dst.chromX      = src.chromX;
    dst.vi          = src.vi;
    /* EAPSO fields */
    dst.fitnessRank = src.fitnessRank;
    dst.entropyRank = src.entropyRank;
    dst.quality     = src.quality;
}

/* Apply mutation to the population */
void Swarm::mutation()
{
    for (unsigned i = 0; i < tPop; ++i) {
        for (unsigned j = 0; j < nAllele; ++j) {
            if (flip(pMut)) {
                switch (representation) {
                    case BINARY:
                        population[i].chromX[j] = population[i].chromX[j] ? 0u : 1u;
                        break;
                    case INTEGER_A:
                    case INTEGER_B:
                        population[i].chromX[j] =
                            rndIR(lInf[j], (lSup[j] + 1) * cardinality) % (lSup[j] + 1);
                        break;
                }
            }
        }
    }
}

/* Record generation statistics */
void Swarm::runInfo(const std::string& file, unsigned gen)
{
    Run.meanFitness    += Gen.meanFitness    / nGen;
    Run.squaredFitness += Gen.squaredFitness / nGen;

    if (gen == 0) {
        copyParticle(Gen.best,  Run.best);
        copyParticle(Gen.worst, Run.worst);
        generationStats(file, gen, *this);
        std::printf("\nGeneration: %04u -- Fitness: %f", gen, Gen.best.fitness);
        Run.generation = gen;
    } else {
        if (Gen.best.fitness > Run.best.fitness) {
            copyParticle(Gen.best, Run.best);
            Run.generation = gen;
        }
        generationStats(file, gen, *this);
        std::printf("\nGeneration: %04u -- Fitness: %f", gen, Gen.best.fitness);

        if (Gen.worst.fitness < Run.worst.fitness)
            copyParticle(Gen.worst, Run.worst);
    }
}

/* Core PSO velocity and position update — unchanged */
void Swarm::PSOAlgorithm(unsigned gen)
{
    for (unsigned i = 0; i < tPop; ++i) {
        for (unsigned d = 0; d < nAllele; ++d) {
            const double phi1p = rndF() * phi1;
            const double phi2p = rndF() * phi2;

            if (gen) {
                population[i].vi[d] +=
                    phi1p * (static_cast<double>(bestIndividualExp[i].chromX[d])
                           - static_cast<double>(population[i].chromX[d]));
                population[i].vi[d] +=
                    phi2p * (static_cast<double>(bestSocialExp[i].chromX[d])
                           - static_cast<double>(population[i].chromX[d]));
            } else {
                if      (!bestSocialExp[i].chromX[d] && !bestIndividualExp[i].chromX[d])
                    population[i].vi[d] = -vMax;
                else if ( bestSocialExp[i].chromX[d] &&  bestIndividualExp[i].chromX[d])
                    population[i].vi[d] =  vMax;
                else
                    population[i].vi[d] =  0.0;

                population[i].vi[d] +=
                    phi1p * static_cast<double>(bestIndividualExp[i].chromX[d]);
                population[i].vi[d] +=
                    phi2p * static_cast<double>(bestSocialExp[i].chromX[d]);
            }

            /* Clamp velocity */
            population[i].vi[d] = std::clamp(population[i].vi[d], -vMax, vMax);

            /* Update position via sigmoid-mapped velocity */
            const double vNorm = sigmoid(population[i].vi[d]);
            switch (representation) {
                case BINARY:
                    population[i].chromX[d] = flip(vNorm);
                    break;
                case INTEGER_A:
                    population[i].chromX[d] = flip(vNorm)
                        ? bestSocialExp[i].chromX[d]
                        : population[i].chromX[d];
                    break;
                case INTEGER_B:
                    population[i].chromX[d] = flip(vNorm)
                        ? bestSocialExp[i].chromX[d]
                        : flip(1.0 - vNorm)
                            ? bestIndividualExp[i].chromX[d]
                            : population[i].chromX[d];
                    break;
            }
        }
    }
}

} // namespace PSwarm