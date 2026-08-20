/********************************************************/
/*                  CINVESTAV - IPN                     */
/*        Department of Electrical Engineering          */
/*                 Computing Section                    */
/*                                                      */
/*               Evolutionary Computation               */
/*                                                      */
/*                    PSO-PSO orchestrator               */
/*                     April 2026                       */
/*                                                      */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* File: psopso.h                                       */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Implements the Parallel Swarms Oriented  */
/* PSO (PSO-PSO) pipeline of Gonsalves & Egashira (2013), */
/* "Parallel Swarms Oriented Particle Swarm Optimization" */
/* (Applied Computational Intelligence and Soft Computing,*/
/* vol. 2013, article 756719), on top of the existing     */
/* PSwarm::Swarm circuit-design solver.                   */
/*                                                        */
/* Algorithm summary (paper Section 4):                   */
/*                                                        */
/*   Multi-evolutionary phase (K sub-swarms, N iterations):*/
/*     Step 1: randomly initialise K independent swarms.  */
/*     Step 2: evaluate fitness of every particle.        */
/*     Step 3: determine pbest and sbest of each swarm.    */
/*     Step 4: update velocity/position, eq. (1)/(2).      */
/*     Step 5: repeat steps 2-4 for N iterations.          */
/*                                                        */
/*   Single-evolutionary phase (M iterations):             */
/*     Step 6: gbest = best sbest across all K swarms.     */
/*     Step 7: update velocity using pbest, sbest, AND      */
/*             gbest, eq. (11).                            */
/*     Step 8: update position, eq. (12).                  */
/*     Step 9: repeat steps 2-8 (with the gbest term) for   */
/*             M iterations.                               */
/*                                                        */
/*   The two phases alternate for a configured number of   */
/*   cycles, after which the best-ever particle (gbest) is */
/*   reported.                                              */
/*                                                        */
/* Parallelism: each of the K sub-swarms is a completely   */
/* independent PSwarm::Swarm (own population, decoder       */
/* buffers, thread-local RNG stream, own output files) --   */
/* exactly like psomatrixcircuit_mt.cpp, but here the K      */
/* swarms are *not* independent runs: after every N          */
/* iterations they synchronise (via a thread join) to        */
/* exchange their swarm-bests and compute a shared global-   */
/* best, then run M more iterations, now steered in part by  */
/* that shared value, before resynchronising again.          */
/********************************************************/

#pragma once

#include <string>

namespace PSwarm {

struct PSOPSOConfig {
    std::string inputFile;

    /* Structural parameters of the PSO-PSO topology. These are
       not "new algorithm parameters" in the paper's sense (they
       don't tune convergence like an inertia weight would) --
       they describe the swarm-of-swarms structure itself, the
       same way population size does for a single swarm. */
    unsigned numSwarms        = 4;   /* K: number of sub-swarms       */
    unsigned multiPhaseIters  = 20;  /* N: multi-evolutionary iters   */
    unsigned singlePhaseIters = 10;  /* M: single-evolutionary iters  */
    unsigned numCycles        = 10;  /* number of (multi,single) alternations */

    /* phi3 / c3: acceleration coefficient toward gbest, used only
       during the single-evolutionary phase (eq. 11). phi1/phi2
       (c1/c2, toward pbest/sbest) are read from the input file
       like every other PSO parameter. */
    double phi3 = 2.0;

    /* 0 => seed from std::random_device once, then derive a
       distinct per-swarm/per-cycle/per-phase seed from it, so
       repeated runs of the program don't correlate. */
    unsigned baseSeed = 0;
};

/* Runs the full PSO-PSO pipeline described above and prints a
   summary (best fitness + Boolean expression of every output)
   to stdout. Per-generation CSVs are written per sub-swarm, per
   cycle, per phase, following this codebase's existing
   statistics.cpp format; a top-level "<input>_psopso_cycles.csv"
   log records how the global-best evolves cycle by cycle. */
void runPSOPSO(const PSOPSOConfig& cfg);

} // namespace PSwarm
