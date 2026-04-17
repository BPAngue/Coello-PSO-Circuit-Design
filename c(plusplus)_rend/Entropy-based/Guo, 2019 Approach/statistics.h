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
/* File: statistics.h                                   */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Header for the per-generation and       */
/* per-run statistics library.                          */
/********************************************************/

#pragma once

#include <string>
#include "particle.h"   /* Particle - separated to break circular deps */

namespace PSwarm {

/* Forward-declare Swarm to avoid pulling in the full
   psomatrixcircuit.h from every statistics translation unit. */
class Swarm;

/*------------------------------------------------------*/
/* Statistics: data accumulated for one generation or   */
/* one complete run.  Replaces the original C struct    */
/* with default member initialisers.                    */
/*------------------------------------------------------*/
struct Statistics {
    double   meanFitness    = 0.0;  /* Mean fitness across the population */
    double   squaredFitness = 0.0;  /* E[f²] — used to compute variance */
    unsigned generation     = 0;    /* Generation when the optimum was found */
    Particle best;                  /* Best particle seen so far */
    Particle worst;                 /* Worst particle seen so far */
};

/* -------------------------------------------------------
 * Free functions
 *
 * Functions that only touch a Statistics struct take it by
 * reference.  Functions that also read Swarm configuration
 * (file names, population size, phi values, etc.) take a
 * const Swarm& as well.
 * ------------------------------------------------------- */

/* Zero-initialise a statistics block */
void initStatistics(Statistics& st);

/* Write the CSV header shared by all runs to filename */
void globalHeader(const std::string& filename, const Swarm& swarm);

/* Append one run's summary row to the global CSV */
void runStatistics(const std::string& filename, unsigned run, Swarm& swarm);

/* Write the per-run CSV header */
void runHeader(const std::string& filename, unsigned seed, const Swarm& swarm);

/* Append one generation row to the per-run CSV */
void generationStats(const std::string& filename, unsigned gen, const Swarm& swarm);

/* Append the final summary to the per-run CSV */
void runFooter(const std::string& filename, Swarm& swarm);

/* Convert a particle's chromosome to a printable string */
std::string chromToString(const Particle& par, const Swarm& swarm);

/* Convert an integer to a string (retained for compatibility;
   internally delegates to std::to_string) */
std::string dtoa(unsigned num);

} // namespace PSwarm
