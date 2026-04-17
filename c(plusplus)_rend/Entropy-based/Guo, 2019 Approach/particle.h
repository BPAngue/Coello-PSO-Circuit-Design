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
/* File: particle.h                                     */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Standalone Particle class definition.   */
/* Separated from psomatrixcircuit.h to break the       */
/* circular dependency between statistics.h and         */
/* psomatrixcircuit.h.                                  */
/*                                                      */
/* EAPSO modification (April 2026):                     */
/*   Added fitnessRank, entropyRank, and quality fields */
/*   to support entropy-assisted exemplar selection.    */
/*   These are computed each generation in              */
/*   Swarm::computeEntropyRanks() and                   */
/*   Swarm::computeQuality() and used during            */
/*   neighbourhood best selection in                    */
/*   Swarm::evaluatePopulation().                       */
/********************************************************/

#pragma once

#include <vector>

namespace PSwarm {

/*------------------------------------------------------*/
/* Particle: one individual in the swarm population.    */
/*------------------------------------------------------*/
class Particle {
public:
    /* PSO data */
    std::vector<double>   vi;       /* Velocity per allele */
    std::vector<unsigned> chromX;   /* Chromosome that evolves */

    /* Evaluation results */
    unsigned numGates   = 0;    /* Gates used in solution */
    unsigned numNoGates = 0;    /* Wire cells in solution */
    unsigned numEqual   = 0;    /* Outputs matching the target truth table */
    double   fitness    = 0.0;  /* Particle fitness */

    /* EAPSO diversity fields --------------------------------
     * fitnessRank  : rank of this particle by raw fitness
     *                (1 = best, tPop = worst).  Assigned in
     *                Swarm::computeEntropyRanks().
     * entropyRank  : rank of the fitness-domain segment this
     *                particle occupies, ordered by segment
     *                occupancy (1 = least crowded segment,
     *                i.e. most novel).  Assigned in
     *                Swarm::computeEntropyRanks().
     * quality      : combined score Q = alpha*fitnessRank
     *                + entropyRank.  Lower is better.
     *                Assigned in Swarm::computeQuality().
     *                Used instead of raw fitness when
     *                selecting the neighbourhood best in
     *                Swarm::evaluatePopulation().
     * ------------------------------------------------------ */
    unsigned fitnessRank = 0;
    unsigned entropyRank = 0;
    double   quality     = 0.0;

    Particle() = default;
    explicit Particle(unsigned nAllele)
        : vi(nAllele, 0.0), chromX(nAllele, 0u) {}
};

} // namespace PSwarm