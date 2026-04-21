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

    Particle() = default;
    explicit Particle(unsigned nAllele)
        : vi(nAllele, 0.0), chromX(nAllele, 0u) {}
};

} // namespace PSwarm
