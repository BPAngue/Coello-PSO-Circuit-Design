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

#ifdef _OPENMP
#  include <omp.h>
#endif

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

    /* ---- PSO-PSO parameters (new) ---- */
    numSubswarms   = static_cast<unsigned>(CircuitData::readNumber(input));
    subswarmIters  = static_cast<unsigned>(CircuitData::readNumber(input));
    phi3           = CircuitData::readNumber(input);

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
    if (numSubswarms > 1) {
        std::printf("\n    PSO-PSO sub-swarms (K): %u", numSubswarms);
        std::printf("\n    Multi-phase iters (N): %u",  subswarmIters);
        std::printf("\n    gbest acceleration (phi3): %.2f", phi3);
    } else {
        std::printf("\n    PSO-PSO: disabled (single swarm)");
    }
    std::printf("\n   Per-run results file: %sX.csv",   nfGen.c_str());
    std::printf("\nGlobal results file: %s\n\n",        nfRun.c_str());
}

/* Initialise derived variables */
void Swarm::initVariables()
{
    if (tNeigh > tPop) tNeigh = tPop;

    /* ---- Normalise PSO-PSO parameters ---- */
    if (numSubswarms == 0)        numSubswarms = 1;
    if (numSubswarms > tPop)      numSubswarms = tPop;
    if (numSubswarms > 1 && subswarmIters == 0) subswarmIters = 1;

    /* Build the contiguous sub-swarm partition (K+1 boundary indices). */
    subswarmStart.assign(numSubswarms + 1, 0u);
    for (unsigned k = 0; k <= numSubswarms; ++k)
        subswarmStart[k] = (tPop * k) / numSubswarms;

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

    /* PSO-PSO structures */
    subswarmBest.assign(numSubswarms, Particle(nAllele));
    reserveParticleMemory(globalBest, true);

    reserveParticleMemory(Gen.best,  true);
    reserveParticleMemory(Gen.worst, true);
    reserveParticleMemory(Run.best,  true);
    reserveParticleMemory(Run.worst, true);

    decoder.reserve(circuit);

    /* Per-sub-swarm decoders: one MatrixDecoder per sub-swarm gives
       every OpenMP thread its own private scratch buffers, so the
       parallel multi-evolutionary phase has zero shared mutable state
       in the inner loop. */
    if (numSubswarms > 1) {
        subDecoders.clear();
        subDecoders.resize(numSubswarms);
        for (unsigned k = 0; k < numSubswarms; ++k) {
            subDecoders[k].numRows = decoder.numRows;
            subDecoders[k].numCols = decoder.numCols;
            subDecoders[k].tMat    = decoder.tMat;
            subDecoders[k].reserve(circuit);
        }

        /* Pre-allocate the per-sub-swarm, per-inner-iter statistics
           buffers used for reduction after each multi-phase barrier. */
        innerStats.assign(numSubswarms, std::vector<InnerStat>(subswarmIters));
        for (unsigned k = 0; k < numSubswarms; ++k) {
            for (unsigned n = 0; n < subswarmIters; ++n) {
                reserveParticleMemory(innerStats[k][n].best,  true);
                reserveParticleMemory(innerStats[k][n].worst, true);
            }
        }
    }
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

    for (auto& d : subDecoders) d.free();
    subDecoders.clear();
    innerStats.clear();

    population.clear();
    bestSocialExp.clear();
    bestIndividualExp.clear();
    subswarmBest.clear();
    subswarmStart.clear();

    lInf.clear();
    lSup.clear();
    bitVariable.clear();
}

/* Execute one complete PSO run.
 *
 * When numSubswarms == 1 this behaves exactly like the original
 * neighbourhood PSO.
 *
 * When numSubswarms > 1 control is delegated to pSwarmPSOPSO, which
 * implements the island-model PSO-PSO algorithm of Gonsalves &
 * Egashira (2013) with OpenMP thread-level parallelism across
 * sub-swarms.
 */
void Swarm::pSwarm(unsigned run)
{
    if (numSubswarms > 1) {
        pSwarmPSOPSO(run);
        return;
    }

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

/* ------------------------------------------------------------------
 * pSwarmPSOPSO -- island-model PSO-PSO with thread-level parallelism.
 *
 *   for each outer cycle m = 0 .. M-1:
 *       #pragma omp parallel for              (one thread per sub-swarm)
 *       for k = 0 .. K-1:
 *           for n = 0 .. N-1:                 (multi-evolutionary phase)
 *               evaluate sub-swarm k, update pbest & sbest[k]
 *               2-term velocity/position update (pbest + sbest)
 *               mutate sub-swarm k
 *               record innerStats[k][n]
 *       // implicit barrier
 *
 *       // Serial: emit the N multi-phase CSV rows, then single phase.
 *       for n = 0 .. N-1:
 *           reduceInnerStats(n)               (into Gen)
 *           runInfo(fileGen, currentGen++)
 *
 *       evaluatePopulation()                  (serial -- updates gbest)
 *       runInfo(fileGen, currentGen++)
 *       PSOSinglePhase()                      (3-term update)
 *       mutation()
 *
 * Every sub-swarm is algorithmically independent during the multi
 * phase (no pbest / sbest sharing), and the single phase is where
 * gbest information is propagated back to every particle.
 *
 * Thread safety is guaranteed by:
 *   - each thread owning an exclusive contiguous range of particles
 *     and the matching slices of bestIndividualExp / bestSocialExp;
 *   - each thread using its own MatrixDecoder instance from
 *     subDecoders[k];
 *   - the RNG being thread_local (see random.cpp).
 * ------------------------------------------------------------------ */
void Swarm::pSwarmPSOPSO(unsigned run)
{
    unsigned seed = initRandom(0);
    std::string fileGen;
    runFileName(run, fileGen);

    initPopulation();
    runHeader(fileGen, seed, *this);

    const unsigned K = numSubswarms;
    const unsigned N = subswarmIters;

#   ifdef _OPENMP
    if (run == 0) {
        int maxT = 0, nProc = 0;
#       pragma omp parallel
        {
#           pragma omp single
            {
                maxT  = omp_get_num_threads();
                nProc = omp_get_num_procs();
            }
        }
        std::printf("[OpenMP] procs=%d default_threads=%d requested K=%u\n",
                    nProc, maxT, K);
    }
#   else
    if (run == 0) std::printf("[OpenMP] NOT ENABLED -- compiled without -fopenmp\n");
#   endif

    /* We will advance a logical "generation counter" that matches the
       original serial PSO for CSV output.  Each outer cycle produces
       N (multi) + 1 (single) rows.  If nGen is not a multiple of
       (N+1) the remainder is played out as a trailing multi-phase. */
    unsigned currentGen = 0;

    while (currentGen < nGen) {
        const bool fullCycle    = (currentGen + N + 1u) <= nGen;
        const unsigned nThisMul = fullCycle
            ? N
            : std::min<unsigned>(N, nGen - currentGen);
        const bool runSingle    = fullCycle;

        /* ---------------- MULTI-EVOLUTIONARY PHASE ---------------- */
        /* Reset per-sub-swarm stat scratch so leftover particle data
           from the previous cycle never leaks into this cycle. */
        for (unsigned k = 0; k < K; ++k) {
            for (unsigned n = 0; n < nThisMul; ++n) {
                innerStats[k][n].meanSum    = 0.0;
                innerStats[k][n].squaredSum = 0.0;
            }
        }

#       ifdef _OPENMP
#       pragma omp parallel for schedule(static) num_threads(static_cast<int>(K))
#       endif
        for (int ki = 0; ki < static_cast<int>(K); ++ki) {
            const unsigned k = static_cast<unsigned>(ki);

            /* Seed the per-thread RNG deterministically from the run
               seed + sub-swarm id + outer cycle counter.  This keeps
               the stream reproducible across threads for a given
               (run, k, cycle) triple. */
            initRandom(seed ^ (0x9E3779B9u * (k + 1u)) ^ (currentGen * 2654435761u));

            for (unsigned n = 0; n < nThisMul; ++n) {
                const bool firstGen = (currentGen == 0u) && (n == 0u);
                evaluateSubswarm(k, n, firstGen, subDecoders[k]);
                PSOAlgorithmSubswarm(k, firstGen);
                mutationSubswarm(k);
            }
        }
        /* Implicit OpenMP barrier here -- all K sub-swarms have now
           completed their N multi-phase iterations independently. */

        /* ------- Serial: emit N CSV rows, one per multi-phase gen ------- */
        for (unsigned n = 0; n < nThisMul; ++n) {
            initStatistics(Gen);
            reduceInnerStats(n);
            runInfo(fileGen, currentGen);
            ++currentGen;
        }

        if (!runSingle || currentGen >= nGen) break;

        /* ---------------- SINGLE-EVOLUTIONARY PHASE ---------------- */
        /* Serial: evaluate the whole population (updates pbest, each
           sub-swarm's sbest, and gbest).  evaluatePopulation uses the
           main `decoder` instance which is fine because we're outside
           the parallel region. */
        initStatistics(Gen);
        evaluatePopulation(currentGen);
        runInfo(fileGen, currentGen);
        ++currentGen;

        PSOSinglePhase(currentGen - 1);
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

/* Evaluate all particles and update personal / social / global bests.
 *
 * Step 1 - evaluate fitness and update personal bests (pbest).
 * Step 2 - update the social best (sbest) for each particle:
 *            numSubswarms == 1 : original index-neighbourhood logic
 *            numSubswarms  > 1 : PSO-PSO sub-swarm best (best pbest
 *                                among all particles in the same
 *                                sub-swarm, kept across generations).
 *          Also updates the global best (gbest) across sub-swarms,
 *          and accumulates generation statistics.
 */
void Swarm::evaluatePopulation(unsigned gen)
{
    /* ---- Step 1 : per-particle fitness + pbest ---- */
    for (unsigned i = 0; i < tPop; ++i) {
        evaluateParticle(population[i], decoder);

        if (gen == 0 || population[i].fitness > bestIndividualExp[i].fitness)
            copyParticle(population[i], bestIndividualExp[i]);
    }

    /* ---- Step 2 : social best ---- */
    if (numSubswarms <= 1) {
        /* Original neighbourhood PSO behaviour. */
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
            Gen.squaredFitness += bestSocialExp[i].fitness
                                * bestSocialExp[i].fitness / tPop;

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
        return;
    }

    /* ---- PSO-PSO path (numSubswarms > 1) ---- */

    /* 2a.  For each sub-swarm find the best pbest and update sbest (kept
            across generations, so we only overwrite on improvement). */
    for (unsigned k = 0; k < numSubswarms; ++k) {
        const unsigned lo = subswarmStart[k];
        const unsigned hi = subswarmStart[k + 1];
        if (lo >= hi) continue;

        unsigned bestIdx = lo;
        for (unsigned i = lo + 1; i < hi; ++i)
            if (bestIndividualExp[i].fitness > bestIndividualExp[bestIdx].fitness)
                bestIdx = i;

        if (gen == 0 || bestIndividualExp[bestIdx].fitness > subswarmBest[k].fitness)
            copyParticle(bestIndividualExp[bestIdx], subswarmBest[k]);
    }

    /* 2b.  Propagate each sub-swarm's sbest into bestSocialExp[i]
            for every particle that belongs to that sub-swarm. */
    for (unsigned k = 0; k < numSubswarms; ++k) {
        const unsigned lo = subswarmStart[k];
        const unsigned hi = subswarmStart[k + 1];
        for (unsigned i = lo; i < hi; ++i)
            copyParticle(subswarmBest[k], bestSocialExp[i]);
    }

    /* 2c.  Update the global best: max over all sbests. */
    unsigned bestSub = 0;
    for (unsigned k = 1; k < numSubswarms; ++k)
        if (subswarmBest[k].fitness > subswarmBest[bestSub].fitness)
            bestSub = k;

    if (gen == 0 || subswarmBest[bestSub].fitness > globalBest.fitness)
        copyParticle(subswarmBest[bestSub], globalBest);

    /* 2d.  Generation statistics (same semantics as the single-swarm
            path: aggregate over each particle's social reference). */
    for (unsigned i = 0; i < tPop; ++i) {
        Gen.meanFitness    += bestSocialExp[i].fitness / tPop;
        Gen.squaredFitness += bestSocialExp[i].fitness
                            * bestSocialExp[i].fitness / tPop;

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

/* Compute fitness of a single particle using the supplied decoder.
   In single-threaded code callers pass this->decoder.  In the parallel
   multi-evolutionary phase each thread passes its own private
   subDecoders[k] so the internal gate / wire scratch buffers do not
   race. */
void Swarm::evaluateParticle(Particle& par, MatrixDecoder& dec)
{
    dec.evaluate(par.chromX, par.numEqual, *this, circuit);

    par.numNoGates = 0;
    if (par.numEqual >= circuit.numTotalOutputs) {
        par.numGates   = dec.countGates(circuit);
        par.numNoGates = dec.tMat - par.numGates;
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

/* Mutate only the particles in sub-swarm k -- used by the parallel
   multi-evolutionary phase where each thread mutates its own disjoint
   range of the population. */
void Swarm::mutationSubswarm(unsigned k)
{
    const unsigned lo = subswarmStart[k];
    const unsigned hi = subswarmStart[k + 1];
    for (unsigned i = lo; i < hi; ++i) {
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

/* ------------------------------------------------------------------
 * evaluateSubswarm  --  one multi-phase evaluation of a single island.
 *
 * Steps:
 *   1. Compute fitness of every particle in sub-swarm k and update its
 *      personal best (pbest).
 *   2. Find the best pbest inside the sub-swarm and update
 *      subswarmBest[k].
 *   3. Propagate subswarmBest[k] to bestSocialExp[i] for every i in
 *      sub-swarm k.  This is what makes the subsequent 2-term velocity
 *      update attract particles toward their sub-swarm best.
 *   4. Record scalar statistics and best/worst chromosome snapshots
 *      into innerStats[k][multiIter % N] so they can be aggregated
 *      after the barrier.
 *
 * Thread safety:
 *   Every index (population[i], bestIndividualExp[i], bestSocialExp[i]
 *   with i in [lo, hi), and subswarmBest[k], innerStats[k][n]) is
 *   owned exclusively by thread k, so there are no races.
 * ------------------------------------------------------------------ */
void Swarm::evaluateSubswarm(unsigned k, unsigned statSlot, bool firstGen,
                             MatrixDecoder& dec)
{
    const unsigned lo = subswarmStart[k];
    const unsigned hi = subswarmStart[k + 1];
    if (lo >= hi) return;

    /* 1. Evaluate + update pbest. */
    for (unsigned i = lo; i < hi; ++i) {
        evaluateParticle(population[i], dec);
        if (firstGen || population[i].fitness > bestIndividualExp[i].fitness)
            copyParticle(population[i], bestIndividualExp[i]);
    }

    /* 2. Best pbest inside the sub-swarm -> subswarmBest[k]. */
    unsigned bestIdx  = lo;
    unsigned worstIdx = lo;
    for (unsigned i = lo + 1; i < hi; ++i) {
        if (bestIndividualExp[i].fitness > bestIndividualExp[bestIdx].fitness)
            bestIdx = i;
        if (bestIndividualExp[i].fitness < bestIndividualExp[worstIdx].fitness)
            worstIdx = i;
    }
    if (firstGen ||
        bestIndividualExp[bestIdx].fitness > subswarmBest[k].fitness)
        copyParticle(bestIndividualExp[bestIdx], subswarmBest[k]);

    /* 3. Propagate sbest -> bestSocialExp[i] for this sub-swarm. */
    for (unsigned i = lo; i < hi; ++i)
        copyParticle(subswarmBest[k], bestSocialExp[i]);

    /* 4. Scalar stats + best/worst snapshot for the reducer. */
    InnerStat& st = innerStats[k][statSlot];
    st.meanSum    = 0.0;
    st.squaredSum = 0.0;
    for (unsigned i = lo; i < hi; ++i) {
        const double f = bestIndividualExp[i].fitness;
        st.meanSum    += f;
        st.squaredSum += f * f;
    }
    copyParticle(bestIndividualExp[bestIdx],  st.best);
    copyParticle(bestIndividualExp[worstIdx], st.worst);
}

/* ------------------------------------------------------------------
 * PSOAlgorithmSubswarm  --  2-term velocity + position update,
 * restricted to sub-swarm k.  Same math as PSOAlgorithm; only the
 * particle index range differs.
 * ------------------------------------------------------------------ */
void Swarm::PSOAlgorithmSubswarm(unsigned k, bool firstGen)
{
    const unsigned lo = subswarmStart[k];
    const unsigned hi = subswarmStart[k + 1];

    for (unsigned i = lo; i < hi; ++i) {
        for (unsigned d = 0; d < nAllele; ++d) {
            const double phi1p = rndF() * phi1;
            const double phi2p = rndF() * phi2;

            if (!firstGen) {
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

            population[i].vi[d] = std::clamp(population[i].vi[d], -vMax, vMax);

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

/* ------------------------------------------------------------------
 * reduceInnerStats  --  combine the K sub-swarm stat records at
 * multi-iteration n into this->Gen so the existing runInfo /
 * generationStats pipeline can emit a CSV row that matches the
 * format of the original serial PSO.
 *
 * Gen.meanFitness / squaredFitness are weighted averages across all
 * tPop particles.  Gen.best / Gen.worst are the global max / min
 * chromosomes taken across all sub-swarms at that iteration.
 * ------------------------------------------------------------------ */
void Swarm::reduceInnerStats(unsigned n)
{
    Gen.meanFitness    = 0.0;
    Gen.squaredFitness = 0.0;

    double totalMean   = 0.0;
    double totalSqr    = 0.0;

    unsigned bestK  = 0;
    unsigned worstK = 0;
    bool     haveSeed = false;

    for (unsigned k = 0; k < numSubswarms; ++k) {
        totalMean += innerStats[k][n].meanSum;
        totalSqr  += innerStats[k][n].squaredSum;

        if (!haveSeed) {
            bestK  = k;
            worstK = k;
            haveSeed = true;
        } else {
            if (innerStats[k][n].best.fitness  > innerStats[bestK ][n].best.fitness)
                bestK  = k;
            if (innerStats[k][n].worst.fitness < innerStats[worstK][n].worst.fitness)
                worstK = k;
        }
    }

    Gen.meanFitness    = totalMean / static_cast<double>(tPop);
    Gen.squaredFitness = totalSqr  / static_cast<double>(tPop);

    copyParticle(innerStats[bestK ][n].best,  Gen.best);
    copyParticle(innerStats[worstK][n].worst, Gen.worst);
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

/* ------------------------------------------------------------------
 * PSOSinglePhase  -  PSO-PSO single-evolutionary phase.
 *
 * Fires once at the end of every (N+1)-generation cycle.  Each
 * particle is attracted toward three references:
 *     - pbest                (its own historical best)
 *     - sbest (sub-swarm)    (stored in bestSocialExp[i])
 *     - gbest (global)       (stored in this->globalBest)
 *
 * Velocity is clamped to [-vMax, vMax] and the discrete position
 * update is the same as PSOAlgorithm so that representation-specific
 * semantics (binary / integer-A / integer-B) are preserved.
 * ------------------------------------------------------------------ */
void Swarm::PSOSinglePhase(unsigned /*gen*/)
{
    for (unsigned i = 0; i < tPop; ++i) {
        for (unsigned d = 0; d < nAllele; ++d) {
            const double phi1p = rndF() * phi1;
            const double phi2p = rndF() * phi2;
            const double phi3p = rndF() * phi3;

            population[i].vi[d] +=
                phi1p * (static_cast<double>(bestIndividualExp[i].chromX[d])
                       - static_cast<double>(population[i].chromX[d]));
            population[i].vi[d] +=
                phi2p * (static_cast<double>(bestSocialExp[i].chromX[d])
                       - static_cast<double>(population[i].chromX[d]));
            population[i].vi[d] +=
                phi3p * (static_cast<double>(globalBest.chromX[d])
                       - static_cast<double>(population[i].chromX[d]));

            population[i].vi[d] = std::clamp(population[i].vi[d], -vMax, vMax);

            const double vNorm = sigmoid(population[i].vi[d]);
            switch (representation) {
                case BINARY:
                    population[i].chromX[d] = flip(vNorm);
                    break;
                case INTEGER_A:
                    /* Attract toward sbest or gbest (whichever wins a coin flip). */
                    if (flip(vNorm)) {
                        population[i].chromX[d] = flip(0.5)
                            ? bestSocialExp[i].chromX[d]
                            : globalBest.chromX[d];
                    }
                    break;
                case INTEGER_B:
                    if (flip(vNorm)) {
                        population[i].chromX[d] = flip(0.5)
                            ? bestSocialExp[i].chromX[d]
                            : globalBest.chromX[d];
                    } else if (flip(1.0 - vNorm)) {
                        population[i].chromX[d] = bestIndividualExp[i].chromX[d];
                    }
                    break;
            }
        }
    }
}

/* Core PSO velocity and position update */
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
