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
/********************************************************/

#include "psomatrixcircuit.h"
#include "random.h"
#include <cstdio>
#include <cmath>
#include <algorithm>
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
    std::printf("\nGlobal results file: %s\n\n",        nfRun.c_str());
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

    initPopulationEntropy();   /* IPSO: diversity-preserving initialisation */
    runHeader(fileGen, seed, *this);

    for (unsigned gen = 0; gen < nGen; ++gen) {
        initStatistics(Gen);
        evaluatePopulation(gen);
        updateAccelerationCoefficients();  /* IPSO: rank-based ci per particle */
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

/* Generate the initial random population (kept for reference) */
void Swarm::initPopulation()
{
    for (unsigned i = 0; i < tPop; ++i)
        for (unsigned j = 0; j < nAllele; ++j)
            population[i].chromX[j] = (representation == BINARY)
                ? flip(0.5)
                : rndIR(lInf[j], (lSup[j] + 1) * cardinality) % (lSup[j] + 1);
}

/* -------------------------------------------------------
 * IPSO §II-A  —  Entropy-based population initialisation
 * (Yang et al. 2020)
 *
 * Idea: for every allele dimension d, compute how "similar"
 * a candidate value is to all already-accepted particles.
 * The similarity probability G_mN between particle m and
 * candidate N in dimension d is:
 *
 *   G_mN^d = |y_m^d - y_N^d| / (lSup[d] - lInf[d])
 *
 * The per-dimension entropy of candidate N is:
 *
 *   S_d = - sum_m ( G_mN^d * ln(G_mN^d) )
 *
 * A candidate allele is accepted only when S_d > entropyThreshold,
 * guaranteeing the initial swarm is spread across the search
 * space rather than clustered in one region.
 * ------------------------------------------------------- */
void Swarm::initPopulationEntropy()
{
    /* Safety fallback: if the range of any allele is zero we
     * cannot compute entropy and fall back to uniform random. */
    bool canUseEntropy = true;
    for (unsigned d = 0; d < nAllele; ++d)
        if (lSup[d] == lInf[d]) { canUseEntropy = false; break; }

    if (!canUseEntropy) {
        std::printf("[IPSO] Warning: zero-range allele detected – "
                    "falling back to uniform random initialisation.\n");
        initPopulation();
        return;
    }

    const unsigned maxRetries = 200;   /* per-allele retry cap */

    for (unsigned i = 0; i < tPop; ++i) {
        for (unsigned d = 0; d < nAllele; ++d) {

            unsigned candidate = 0;
            bool     accepted  = false;

            for (unsigned attempt = 0; attempt < maxRetries; ++attempt) {

                /* Draw a random candidate value for allele d */
                candidate = (representation == BINARY)
                    ? flip(0.5)
                    : rndIR(lInf[d], (lSup[d] + 1) * cardinality) % (lSup[d] + 1);

                /* For the very first particle there are no prior particles
                 * to compare against, so accept immediately.             */
                if (i == 0) { accepted = true; break; }

                /* Compute the per-dimension entropy for this candidate
                 * relative to all previously accepted particles.         */
                double range = static_cast<double>(lSup[d] - lInf[d]);
                double entropy = 0.0;

                for (unsigned m = 0; m < i; ++m) {
                    double diff = std::fabs(
                        static_cast<double>(population[m].chromX[d])
                        - static_cast<double>(candidate));

                    double G = diff / range;   /* similarity probability */

                    /* G == 0 means identical values: entropy contribution
                     * is 0*ln(0) = 0 by convention (L'Hôpital), so skip. */
                    if (G > 1e-12)
                        entropy -= G * std::log(G);
                }

                if (entropy > entropyThreshold) { accepted = true; break; }
            }

            /* If no candidate passed the entropy test within the retry
             * budget, use the last drawn value to avoid an infinite loop. */
            if (!accepted)
                std::printf("[IPSO] Entropy threshold not met for particle %u "
                            "allele %u after %u retries; using last candidate.\n",
                            i, d, maxRetries);

            population[i].chromX[d] = candidate;
        }
    }
}

/* -------------------------------------------------------
 * IPSO §II-B  —  Fitness-based per-particle acceleration
 * (Yang et al. 2020)
 *
 * After each generation's evaluation:
 *   Fi = (fit_i - fit_min) / (fit_max - fit_min)
 *   ci = c1 + (c2 - c1) * Fi
 *
 * Fi == 0 -> best particle  -> ci = c1 (slow, local search)
 * Fi == 1 -> worst particle -> ci = c2 (fast, global search)
 *
 * This replaces the fixed phi1/phi2 used by standard PSO so
 * that every particle has a learning rate that is adapted to
 * its current quality, balancing exploration and exploitation.
 * ------------------------------------------------------- */
void Swarm::updateAccelerationCoefficients()
{
    /* Collect fitness values from personal-best experiences */
    double fitMin =  bestIndividualExp[0].fitness;
    double fitMax =  bestIndividualExp[0].fitness;
    for (unsigned i = 1; i < tPop; ++i) {
        if (bestIndividualExp[i].fitness < fitMin) fitMin = bestIndividualExp[i].fitness;
        if (bestIndividualExp[i].fitness > fitMax) fitMax = bestIndividualExp[i].fitness;
    }

    double range = fitMax - fitMin;

    for (unsigned i = 0; i < tPop; ++i) {
        double Fi = (range > 1e-12)
            ? (bestIndividualExp[i].fitness - fitMin) / range
            : 0.5;   /* all particles identical: give everyone the midpoint */

        /* Note: higher fitness is *better* in this codebase, so a particle
         * with Fi == 1.0 is the BEST.  We invert so that the best particle
         * gets the smallest coefficient (local search) and the worst gets
         * the largest (global search), matching the paper's intent.        */
        population[i].ci = c1 + (c2 - c1) * (1.0 - Fi);
    }
}

/* Evaluate all particles and update personal/social bests */
void Swarm::evaluatePopulation(unsigned gen)
{
    /* Step 1 – evaluate fitness and update personal bests */
    for (unsigned i = 0; i < tPop; ++i) {
        evaluateParticle(population[i]);

        if (gen == 0 || population[i].fitness > bestIndividualExp[i].fitness)
            copyParticle(population[i], bestIndividualExp[i]);
    }

    /* Step 2 – find neighbourhood bests, update social bests and gen stats */
    for (unsigned i = 0; i < tPop; ++i) {
        int p  = static_cast<int>(i);
        int p0 = p;
        int p1 = p;

        for (unsigned j = 0; j < tNeigh; j += 2) {
            p0 = (p0 + 1 < static_cast<int>(tPop)) ? p0 + 1 : 0;
            p1 = (p1 - 1 >= 0) ? p1 - 1 : static_cast<int>(tPop) - 1;

            if (bestIndividualExp[p0].fitness > bestIndividualExp[p].fitness) p = p0;
            if (bestIndividualExp[p1].fitness > bestIndividualExp[p].fitness) p = p1;
        }

        copyParticle(bestIndividualExp[p], bestSocialExp[i]);

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
    dst.numGates   = src.numGates;
    dst.numEqual   = src.numEqual;
    dst.numNoGates = src.numNoGates;
    dst.fitness    = src.fitness;
    dst.ci         = src.ci;       /* IPSO: carry per-particle coefficient */
    dst.chromX     = src.chromX;
    dst.vi         = src.vi;
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

/* Core PSO velocity and position update
 * IPSO modification: phi1 and phi2 are both replaced by the
 * per-particle acceleration coefficient ci (Yang et al. 2020).
 * ci was computed by updateAccelerationCoefficients() just
 * before this function is called each generation.             */
void Swarm::PSOAlgorithm(unsigned gen)
{
    for (unsigned i = 0; i < tPop; ++i) {

        /* IPSO: use this particle's individual learning rate.
         * Random scaling [0, ci] is applied to each attraction
         * term, matching equation (7) in the paper.           */
        const double ci = population[i].ci;

        for (unsigned d = 0; d < nAllele; ++d) {
            const double phi1p = rndF() * ci;   /* cognitive: scaled by ci */
            const double phi2p = rndF() * ci;   /* social:    scaled by ci */

            if (gen) {
                population[i].vi[d] +=
                    phi1p * (static_cast<double>(bestIndividualExp[i].chromX[d])
                           - static_cast<double>(population[i].chromX[d]));
                population[i].vi[d] +=
                    phi2p * (static_cast<double>(bestSocialExp[i].chromX[d])
                           - static_cast<double>(population[i].chromX[d]));
            } else {
                /* Generation 0: initialise velocity with a bias toward bests */
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