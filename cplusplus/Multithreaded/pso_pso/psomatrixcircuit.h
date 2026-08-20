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
/* File: psomatrixcircuit.h                             */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Header for the PSO circuit design       */
/* library.  The Swarm class owns all PSO state,        */
/* circuit data, matrix decoder, and statistics.        */
/********************************************************/

#pragma once

#include <string>
#include <vector>
#include <cstdlib>

#include "particle.h"    /* Particle class */
#include "circuits.h"    /* CircuitData struct + gate constants */
#include "matrixpso.h"   /* MatrixDecoder struct */
#include "statistics.h"  /* Statistics struct + free functions */

namespace PSwarm {

/* -------------------------------------------------------
 * Representation type constants
 * ------------------------------------------------------- */
constexpr unsigned BINARY    = 2;
constexpr unsigned INTEGER_A = 1;
constexpr unsigned INTEGER_B = 0;

/*------------------------------------------------------*/
/* Swarm: the central class that owns all program state  */
/* that was previously spread across global variables.   */
/*------------------------------------------------------*/
class Swarm {
public:

    /* ---- Owned sub-systems ---- */
    CircuitData   circuit;   /* Truth table + circuit meta-data */
    MatrixDecoder decoder;   /* Matrix decode / evaluate engine */
    Statistics    Gen;       /* Current generation statistics */
    Statistics    Run;       /* Current run statistics */

    /* ---- File names ---- */
    std::string nfInput;     /* Input data file */
    std::string nfRun;       /* Global run statistics CSV */
    std::string nfGen;       /* Per-run generation statistics prefix */

    /* ---- User-specified PSO parameters ---- */
    unsigned tPop           = 0;    /* Population size */
    unsigned nGen           = 0;    /* Max generations per run */
    unsigned nRun           = 0;    /* Number of independent runs */
    unsigned tNeigh         = 0;    /* Neighbourhood size */
    double   phi1           = 0.0;  /* Cognitive acceleration */
    double   phi2           = 0.0;  /* Social acceleration */
    double   vMax           = 0.0;  /* Maximum velocity */
    double   pMut           = 0.0;  /* Mutation probability */
    unsigned representation = 0;    /* BINARY=2, INTEGER_A=1, INTEGER_B=0 */
    unsigned cardinality    = 0;    /* Cardinality for integer representations */

    /* ---- Problem variables ---- */
    unsigned nVar    = 0;   /* Number of chromosome variables */
    unsigned nAllele = 0;   /* Total alleles (bits or integers) */
    unsigned instanceSeed = 0; /* Random seed for this run */

    std::vector<unsigned> lInf;         /* Lower bounds per variable */
    std::vector<unsigned> lSup;         /* Upper bounds per variable */
    std::vector<unsigned> bitVariable;  /* Bits per variable (binary only) */

    /* ---- Populations ---- */
    std::vector<Particle> population;           /* Working population */
    std::vector<Particle> bestSocialExp;        /* Best social experiences */
    std::vector<Particle> bestIndividualExp;    /* Best individual experiences */

    /* ---- Public interface ---- */

    /* Load all parameters + truth table from file */
    bool loadParameters(const std::string& filename);

    /* Print loaded parameters to stdout */
    void printParameters() const;

    /* Initialise derived variables (tMat, bounds, allele count) */
    void initVariables();

    /* Allocate working memory for all populations */
    void reserveMemory();

    /* Release all working memory */
    void freeMemory();

    /* Execute one complete PSO run */
    void pSwarm(unsigned runIndex);

    /* Build the per-run CSV filename into outName */
    void runFileName(unsigned runIndex, std::string& outName) const;

    /* Generate the initial random population */
    void initPopulation();

    /* Evaluate all particles and update experiences */
    void evaluatePopulation(unsigned generation);

    /* Compute fitness for a single particle */
    void evaluateParticle(Particle& p);

    /* Deep-copy src into dst */
    void copyParticle(const Particle& src, Particle& dst);

    /* Apply mutation to the working population */
    void mutation();

    /* Record generation statistics and update run best/worst */
    void runInfo(const std::string& filename, unsigned generation);

    /* Core PSO velocity + position update */
    void PSOAlgorithm(unsigned generation);

    /* -----------------------------------------------------------
     * PSO-PSO extensions (Gonsalves & Egashira, 2013 -- "Parallel
     * Swarms Oriented Particle Swarm Optimization").
     *
     * These add a *swarm-best* (sBest) on top of the existing
     * per-particle pbest (bestIndividualExp): sBest is the best
     * particle this sub-swarm has ever produced, and stands in
     * for "pg" in the paper's eq. (1) during the multi-evolutionary
     * phase. During the single-evolutionary phase, a *global-best*
     * (gBest) computed across every sub-swarm is additionally
     * used, matching eq. (11).
     *
     * These methods are independent of PSOAlgorithm()/
     * evaluatePopulation() above (which implement this codebase's
     * original ring-neighbourhood social topology) so that both
     * the single-swarm program and the PSO-PSO orchestrator
     * (psopso.h/.cpp) can share the same Swarm class.
     * ----------------------------------------------------------- */

    /* Swarm-best: best particle found by this sub-swarm across
       its whole PSO-PSO run (persists across phases/cycles,
       unlike Gen.best/Run.best which are scoped to one pSwarm() run). */
    Particle sBest;

    /* Size sBest's vectors and set its fitness to -infinity so the
       first evaluated particle always becomes the swarm-best. */
    void initSwarmBest();

    /* Evaluate the population, update each particle's pbest
       (bestIndividualExp), this swarm's sBest, and Gen statistics.
       Does NOT use the ring-neighbourhood social topology. */
    void evaluatePopulationPSOPSO(unsigned generation);

    /* Velocity + position update using pbest and sBest (eq. 1/2,
       multi-evolutionary phase) and, when singlePhase is true,
       additionally gBest (eq. 11/12, single-evolutionary phase).
       isFirstGeneration selects the same special-case velocity
       initialisation PSOAlgorithm() uses at generation 0. */
    void PSOAlgorithmPSOPSO(unsigned generation, bool isFirstGeneration,
                             bool singlePhase, const Particle& gBest, double phi3);

private:
    /* Initialise variable bounds and calculate nAllele */
    void initBounds();

    /* Initialise the vi/chromX vectors of a statistics particle */
    void reserveParticleMemory(Particle& par, bool withChrom);
};

} // namespace PSwarm
