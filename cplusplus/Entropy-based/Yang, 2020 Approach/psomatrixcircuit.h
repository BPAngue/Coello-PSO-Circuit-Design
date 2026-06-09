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
    double   phi1           = 0.0;  /* Cognitive acceleration (kept for BINARY gen-0 init) */
    double   phi2           = 0.0;  /* Social acceleration   (kept for BINARY gen-0 init) */
    double   vMax           = 0.0;  /* Maximum velocity */
    double   pMut           = 0.0;  /* Mutation probability */
    unsigned representation = 0;    /* BINARY=2, INTEGER_A=1, INTEGER_B=0 */
    unsigned cardinality    = 0;    /* Cardinality for integer representations */

    /* ---- IPSO parameters (Yang et al. 2020) ---- */
    /* Lower and upper bounds for the per-particle acceleration coefficient.
     * Each particle's ci = c1 + (c2-c1)*Fi  where Fi is its normalised
     * fitness rank (0 = best, 1 = worst).
     * Worst particle -> ci -> c2 (aggressive global search).
     * Best  particle -> ci -> c1 (fine-grained local search).         */
    double   c1             = 1.5;  /* Lower bound of acceleration range */
    double   c2             = 2.5;  /* Upper bound of acceleration range */

    /* Entropy threshold for population initialisation.
     * A candidate allele is only accepted when its contribution raises
     * the running per-dimension entropy above this value, ensuring the
     * initial swarm spans the search space evenly.                    */
    double   entropyThreshold = 0.5;

    /* ---- Problem variables ---- */
    unsigned nVar    = 0;   /* Number of chromosome variables */
    unsigned nAllele = 0;   /* Total alleles (bits or integers) */

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

    /* IPSO: entropy-based population initialisation (Yang et al. 2020, §II-A)
     * Replaces pure-random seeding; ensures alleles are spread evenly
     * across the search space so the swarm starts with high diversity. */
    void initPopulationEntropy();

    /* IPSO: compute and assign per-particle acceleration coefficients
     * based on individual fitness rank (Yang et al. 2020, §II-B).
     * Must be called after evaluatePopulation() each generation.      */
    void updateAccelerationCoefficients();

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

private:
    /* Initialise variable bounds and calculate nAllele */
    void initBounds();

    /* Initialise the vi/chromX vectors of a statistics particle */
    void reserveParticleMemory(Particle& par, bool withChrom);
};

} // namespace PSwarm