/********************************************************/
/*  Parallel Discrete EPSO with Island Model              */
/*  and Particle Speed Limit                              */
/*                                                        */
/*  Implements the pipeline of:                            */
/*    Ikegami, H. and Mori, H. (2018). "Development of     */
/*    DEPSO Island Model with Particle Speed Limit for     */
/*    Distribution Network Reconfigurations."              */
/*    IFAC PapersOnLine 51-28, pp. 552-557.                */
/*                                                        */
/*  File: epso_island.h                                    */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*  Design notes                                            */
/*  -----------                                             */
/*  * This is an ADDITIVE extension: it does not modify    */
/*    Swarm, Particle, MatrixDecoder, or CircuitData. It   */
/*    reuses swarm.decoder / swarm.circuit for fitness      */
/*    evaluation and reuses the existing statistics free-   */
/*    functions (globalHeader/runHeader/generationStats/    */
/*    runFooter/runStatistics) so the CSV output format     */
/*    matches the original pipeline.                        */
/*                                                        */
/*  * The paper's chromosome coding is a flat bit-string    */
/*    decoded through a sigmoid (eqns 14-15 / 20-21), so    */
/*    this class requires swarm.representation == BINARY.  */
/*                                                        */
/*  * EPSO (eqns 8-11 / 16-19) mutates the acceleration     */
/*    weights w1, w2 with Gaussian noise and mutates gbest  */
/*    with Gaussian noise scaled by tauPrime. Because our   */
/*    gbest is a discrete chromosome (not a continuous      */
/*    vector), the gbest mutation is realised as an extra   */
/*    additive noise term on the velocity update rather     */
/*    than a perturbation of gbest's bits directly. This    */
/*    preserves the "gbest randomisation" role the term      */
/*    plays in eqn (9)/(17) without corrupting the           */
/*    discrete chromosome itself.                            */
/*                                                        */
/*  * Selection uses a (parent + replicationRate children)  */
/*    scheme: the parent's current position is always one  */
/*    of the candidates, so a generation can never regress. */
/*    This is a standard (mu+lambda) ES safeguard and is     */
/*    consistent with "select the best solution... and      */
/*    regard it as the initial solution at the next step"    */
/*    described in section 3 for EPSO.                       */
/*                                                        */
/*  * Migration (Step 3 / Fig. 3-4) uses a ring topology:   */
/*    each island sends its top (migrationRate * size)       */
/*    particles to the next island, replacing that          */
/*    island's worst particles of the same count.            */
/********************************************************/

#pragma once

#include <vector>
#include <string>
#include "psomatrixcircuit.h"

namespace PSwarm {

/*------------------------------------------------------*/
/* IslandParams: tunable parameters for the PDEPSO        */
/* pipeline. Defaults follow Table 1 of the paper          */
/* (Methods C / F, i.e. PDEPSO and PDEPSO with speed       */
/* limit).                                                 */
/*------------------------------------------------------*/
struct IslandParams {
    unsigned nIslands          = 4;     /* No. of islands                        */
    unsigned migrationInterval = 400;   /* Generations between migrations        */
    double   migrationRate     = 0.1;   /* Fraction of each island's pop migrated */
    unsigned replicationRate   = 2;     /* EPSO children produced per particle    */
    double   tau               = 0.01;  /* Learning rate for w1/w2 mutation        */
    double   tauPrime          = 0.02;  /* Learning rate for gbest-noise mutation  */
};

/*------------------------------------------------------*/
/* Island: one independent sub-population -- corresponds  */
/* to one "processor" in the paper's Fig. 4. Each island   */
/* owns its own EPSO weights, personal bests, and          */
/* island-level best (gbest), and evolves independently    */
/* between migrations.                                     */
/*------------------------------------------------------*/
struct Island {
    std::vector<Particle> population;  /* Xi  - current particle positions */
    std::vector<double>   w1;          /* per-particle mutated weight w*_i1 */
    std::vector<double>   w2;          /* per-particle mutated weight w*_i2 */
    std::vector<Particle> pbest;       /* personal best per particle        */
    Particle               gbest;      /* island-level best particle        */
    bool                    hasGbest = false;
};

/*------------------------------------------------------*/
/* IslandPDEPSO: drives the full island-model generation  */
/* loop (Steps 1-3 of section 4.4) for one Swarm instance. */
/* The Swarm supplies circuit/decoder/bounds/tPop/nGen/    */
/* vMax/file-name configuration; this class supplies the   */
/* EPSO + island + speed-limit algorithm itself.            */
/*------------------------------------------------------*/
class IslandPDEPSO {
public:
    IslandPDEPSO(Swarm& swarm, const IslandParams& params = IslandParams{});

    /* Run swarm.nRun independent island-model runs, writing */
    /* the same global + per-run CSV files Swarm::pSwarm      */
    /* would produce.                                          */
    void execute();

private:
    Swarm&        swarm_;
    IslandParams  params_;
    std::vector<Island> islands_;

    /* ---- Setup (Step 1 / Step 2-1) ---- */
    void partitionPopulation();
    void initIsland(Island& isl, unsigned size);

    /* ---- Per-generation island update (Steps 2-2..2-6) ---- */
    void stepIsland(Island& isl);

    /* ---- Migration (Step 3 / Fig. 3-4) ---- */
    void migrate();

    /* ---- Bookkeeping reused by statistics.h free functions ---- */
    Particle globalBest() const;
    Particle globalWorst() const;
    void collectGenerationStats();

    /* ---- Helpers ---- */
    void evaluate(Particle& p) const;
    Particle decodeCandidate(const std::vector<double>& v) const;
};

} // namespace PSwarm
