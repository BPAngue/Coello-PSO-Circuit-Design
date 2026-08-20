/********************************************************/
/* File: epso_island.cpp                                 */
/* See epso_island.h for design notes and references.    */
/********************************************************/

#include "epso_island.h"
#include "random.h"

#include <algorithm>
#include <numeric>
#include <cmath>
#include <cstdio>

namespace PSwarm {

namespace {

constexpr double kTwoPi = 6.283185307179586;

/* Standard Gaussian noise N(0,1) via Box-Muller, built on
   top of the existing Random::rndF() uniform generator so
   no new global RNG state is introduced. */
double gaussianNoise()
{
    double u1 = std::max(rndF(), 1e-12); /* avoid log(0) */
    double u2 = rndF();
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(kTwoPi * u2);
}

} // namespace

/* -------------------------------------------------------
 * Construction
 * ------------------------------------------------------- */
IslandPDEPSO::IslandPDEPSO(Swarm& swarm, const IslandParams& params)
    : swarm_(swarm), params_(params)
{
}

/* -------------------------------------------------------
 * evaluate
 * Delegates fitness computation to the Swarm's existing
 * MatrixDecoder, exactly like Swarm::evaluatePopulation
 * does for the non-island pipeline.
 * ------------------------------------------------------- */
void IslandPDEPSO::evaluate(Particle& p) const
{
    swarm_.evaluateParticle(p);
}

/* -------------------------------------------------------
 * decodeCandidate
 * Builds a new Particle from a velocity vector using the
 * sigmoid + Bernoulli-flip rule of eqns (14)-(15)/(20)-(21).
 * Binary representation only -- see class-level notes.
 * ------------------------------------------------------- */
Particle IslandPDEPSO::decodeCandidate(const std::vector<double>& v) const
{
    Particle p(swarm_.nAllele);
    p.vi = v;
    for (unsigned d = 0; d < swarm_.nAllele; ++d) {
        const double s = sigmoid(v[d]);
        p.chromX[d] = flip(s);
    }
    return p;
}

/* -------------------------------------------------------
 * partitionPopulation (Step 1 / Step 2)
 * Splits swarm.tPop particles as evenly as possible across
 * nIslands sub-populations and initialises each.
 * ------------------------------------------------------- */
void IslandPDEPSO::partitionPopulation()
{
    islands_.assign(params_.nIslands, Island{});

    const unsigned base      = swarm_.tPop / params_.nIslands;
    const unsigned remainder = swarm_.tPop % params_.nIslands;

    for (unsigned isl = 0; isl < params_.nIslands; ++isl) {
        const unsigned size = base + (isl < remainder ? 1u : 0u);
        initIsland(islands_[isl], size);
    }
}

/* -------------------------------------------------------
 * initIsland (Step 2-1)
 * Random binary initialisation, then evaluate and seed
 * pbest / gbest for the island.
 * ------------------------------------------------------- */
void IslandPDEPSO::initIsland(Island& isl, unsigned size)
{
    isl.population.assign(size, Particle(swarm_.nAllele));
    isl.pbest.assign(size, Particle(swarm_.nAllele));
    isl.w1.assign(size, swarm_.phi1);
    isl.w2.assign(size, swarm_.phi2);
    isl.hasGbest = false;

    for (unsigned i = 0; i < size; ++i) {
        for (unsigned d = 0; d < swarm_.nAllele; ++d)
            isl.population[i].chromX[d] = flip(0.5);

        evaluate(isl.population[i]);
        isl.pbest[i] = isl.population[i];

        if (!isl.hasGbest || isl.population[i].fitness > isl.gbest.fitness) {
            isl.gbest    = isl.population[i];
            isl.hasGbest = true;
        }
    }
}

/* -------------------------------------------------------
 * stepIsland (Steps 2-2..2-6, one generation)
 * For every particle: replicate `replicationRate` mutated
 * children (EPSO weight + gbest-noise mutation, eqns 16-19),
 * update velocity/position with the particle speed limit
 * (eqn 22), evaluate, and keep the best of {parent, children}.
 * ------------------------------------------------------- */
void IslandPDEPSO::stepIsland(Island& isl)
{
    const unsigned n = static_cast<unsigned>(isl.population.size());

    for (unsigned i = 0; i < n; ++i) {

        /* (mu + lambda) baseline: the parent is always a
           candidate, so fitness can never regress. */
        Particle bestChild = isl.population[i];
        double    bestW1    = isl.w1[i];
        double    bestW2    = isl.w2[i];

        for (unsigned r = 0; r < params_.replicationRate; ++r) {
            /* Eqn (18): mutate acceleration weights */
            const double w1r = std::max(0.0, isl.w1[i] + params_.tau * gaussianNoise());
            const double w2r = std::max(0.0, isl.w2[i] + params_.tau * gaussianNoise());

            /* Eqn (19): gbest mutation, realised as extra
               velocity noise (see header design notes) */
            const double gbestNoise = params_.tauPrime * gaussianNoise();

            std::vector<double> vNew(swarm_.nAllele);
            for (unsigned d = 0; d < swarm_.nAllele; ++d) {
                double v = isl.population[i].vi[d]
                    + w1r * (static_cast<double>(isl.pbest[i].chromX[d])
                           - static_cast<double>(isl.population[i].chromX[d]))
                    + w2r * (static_cast<double>(isl.gbest.chromX[d])
                           - static_cast<double>(isl.population[i].chromX[d]))
                    + gbestNoise;

                /* Eqn (22): particle speed limit */
                v = std::clamp(v, -swarm_.vMax, swarm_.vMax);
                vNew[d] = v;
            }

            Particle child = decodeCandidate(vNew);
            evaluate(child);

            if (child.fitness > bestChild.fitness) {
                bestChild = child;
                bestW1    = w1r;
                bestW2    = w2r;
            }
        }

        /* Step 2-5: replace particle with the selected candidate */
        isl.population[i] = bestChild;
        isl.w1[i] = bestW1;
        isl.w2[i] = bestW2;

        /* Step 2-6: update pbest / gbest */
        if (bestChild.fitness > isl.pbest[i].fitness)
            isl.pbest[i] = bestChild;

        if (!isl.hasGbest || bestChild.fitness > isl.gbest.fitness) {
            isl.gbest    = bestChild;
            isl.hasGbest = true;
        }
    }
}

/* -------------------------------------------------------
 * migrate (Step 3 / Fig. 3 - ring migration)
 * Each island sends its top (migrationRate * size) particles
 * to the next island in the ring, replacing that island's
 * worst particles of the same count. Index-aligned so that
 * w1/w2/pbest stay consistent with the replaced population
 * slots.
 * ------------------------------------------------------- */
void IslandPDEPSO::migrate()
{
    if (params_.nIslands < 2) return;

    std::vector<std::vector<Particle>> outgoing(params_.nIslands);

    /* Collect each island's best migrants (by value, no index
       alignment needed since these are just copies to send). */
    for (unsigned isl = 0; isl < params_.nIslands; ++isl) {
        auto sorted = islands_[isl].population;
        std::sort(sorted.begin(), sorted.end(),
                   [](const Particle& a, const Particle& b) {
                       return a.fitness > b.fitness;
                   });

        unsigned nMigrants = static_cast<unsigned>(
            std::lround(params_.migrationRate * sorted.size()));
        nMigrants = std::min(nMigrants, static_cast<unsigned>(sorted.size()));

        outgoing[isl].assign(sorted.begin(), sorted.begin() + nMigrants);
    }

    /* Insert incoming migrants into the worst slots of the next
       island, keeping w1/w2/pbest arrays index-aligned with
       the (unsorted) population array. */
    for (unsigned isl = 0; isl < params_.nIslands; ++isl) {
        const unsigned src = (isl == 0) ? params_.nIslands - 1 : isl - 1;
        auto& incoming = outgoing[src];
        if (incoming.empty()) continue;

        auto& pop = islands_[isl].population;

        std::vector<unsigned> idx(pop.size());
        std::iota(idx.begin(), idx.end(), 0u);
        const unsigned nReplace = std::min(static_cast<unsigned>(incoming.size()),
                                            static_cast<unsigned>(pop.size()));
        std::partial_sort(idx.begin(), idx.begin() + nReplace, idx.end(),
                            [&](unsigned a, unsigned b) {
                                return pop[a].fitness < pop[b].fitness;
                            });

        for (unsigned k = 0; k < nReplace; ++k) {
            const unsigned target = idx[k];
            pop[target]                     = incoming[k];
            islands_[isl].pbest[target]     = incoming[k];
            islands_[isl].w1[target]        = swarm_.phi1;
            islands_[isl].w2[target]        = swarm_.phi2;

            if (incoming[k].fitness > islands_[isl].gbest.fitness)
                islands_[isl].gbest = incoming[k];
        }
    }
}

/* -------------------------------------------------------
 * globalBest / globalWorst
 * Best/worst particle across all islands combined.
 * ------------------------------------------------------- */
Particle IslandPDEPSO::globalBest() const
{
    Particle best = islands_[0].gbest;
    for (unsigned isl = 1; isl < params_.nIslands; ++isl)
        if (islands_[isl].gbest.fitness > best.fitness)
            best = islands_[isl].gbest;
    return best;
}

Particle IslandPDEPSO::globalWorst() const
{
    Particle worst = islands_[0].population[0];
    for (const auto& isl : islands_)
        for (const auto& p : isl.population)
            if (p.fitness < worst.fitness)
                worst = p;
    return worst;
}

/* -------------------------------------------------------
 * collectGenerationStats
 * Populates swarm.Gen exactly like Swarm::evaluatePopulation
 * would, so the existing generationStats()/runInfo() free
 * functions can be reused unmodified for CSV output.
 * ------------------------------------------------------- */
void IslandPDEPSO::collectGenerationStats()
{
    initStatistics(swarm_.Gen);

    unsigned total = 0;
    for (const auto& isl : islands_) {
        for (const auto& p : isl.population) {
            swarm_.Gen.meanFitness    += p.fitness;
            swarm_.Gen.squaredFitness += p.fitness * p.fitness;
            ++total;
        }
    }
    if (total) {
        swarm_.Gen.meanFitness    /= total;
        swarm_.Gen.squaredFitness /= total;
    }

    swarm_.Gen.best  = globalBest();
    swarm_.Gen.worst = globalWorst();
}

/* -------------------------------------------------------
 * execute
 * Full algorithm of section 4.4: for each of swarm.nRun
 * independent runs, partition into islands, iterate
 * generations (mutate/update/evaluate/select + periodic
 * migration), and write output via the existing
 * statistics.h free functions.
 * ------------------------------------------------------- */
void IslandPDEPSO::execute()
{
    if (swarm_.representation != BINARY) {
        std::printf(
            "\nError: IslandPDEPSO requires swarm.representation == BINARY "
            "(the paper's coding is a bit-string decoded via sigmoid).\n");
        return;
    }

    globalHeader(swarm_.nfRun, swarm_);

    for (unsigned run = 0; run < swarm_.nRun; ++run) {
        std::printf("\n\nIsland-PDEPSO Run %02u Started\n", run);

        initStatistics(swarm_.Run);
        const unsigned seed = initRandom(0);

        std::string fileGen;
        swarm_.runFileName(run, fileGen);

        partitionPopulation();
        runHeader(fileGen, seed, swarm_);

        for (unsigned gen = 0; gen < swarm_.nGen; ++gen) {
            for (auto& isl : islands_)
                stepIsland(isl);

            collectGenerationStats();
            swarm_.runInfo(fileGen, gen);   /* also writes generationStats() row */

            if (params_.nIslands > 1 &&
                (gen + 1) % params_.migrationInterval == 0)
                migrate();
        }

        runFooter(fileGen, swarm_);
        runStatistics(swarm_.nfRun, run, swarm_);

        std::printf("\n\nIsland-PDEPSO Run %02u Finished\n", run);
    }
}

} // namespace PSwarm
