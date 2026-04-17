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
/*                                                      */
/* EAPSO modification (April 2026):                     */
/*   Added numSegments and alpha constants for the       */
/*   entropy-assisted diversity mechanism from:         */
/*     Guo et al., "An Entropy-Assisted Particle Swarm  */
/*     Optimizer for Large-Scale Optimization Problem"  */
/*     Mathematics 2019, 7, 414.                        */
/*   Added computeEntropyRanks() and computeQuality()   */
/*   private methods.  evaluatePopulation() now uses    */
/*   particle quality Q (not raw fitness) when          */
/*   selecting neighbourhood bests, while all           */
/*   reporting metrics continue to use raw fitness.     */
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

/* -------------------------------------------------------
 * EAPSO constants (hardcoded as per paper Section 4)
 *
 * EAPSO_NUM_SEGMENTS: number of fitness-domain bins used
 *   when computing Shannon entropy H over the population.
 *   Paper uses m = 30.
 *
 * EAPSO_ALPHA: weight on fitnessRank in the quality score
 *   Q = alpha * fitnessRank + entropyRank.
 *   Paper tests alpha in {0.1, 0.2, 0.3, 0.4}; best
 *   results were at 0.1 and 0.2, so 0.1 is used here.
 *   A smaller alpha means diversity (entropyRank) has
 *   more influence on exemplar selection, which is the
 *   key lever against premature convergence.
 * ------------------------------------------------------- */
constexpr unsigned EAPSO_NUM_SEGMENTS = 30u;
constexpr double   EAPSO_ALPHA        = 0.05;

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

    /* Evaluate all particles and update experiences.
     * EAPSO modification: after raw fitness evaluation,
     * calls computeEntropyRanks() and computeQuality()
     * so that neighbourhood best selection uses Q scores. */
    void evaluatePopulation(unsigned generation);

    /* Compute fitness for a single particle */
    void evaluateParticle(Particle& p);

    /* Deep-copy src into dst */
    void copyParticle(const Particle& src, Particle& dst);

    /* Apply mutation to the working population */
    void mutation();

    /* Record generation statistics and update run best/worst */
    void runInfo(const std::string& filename, unsigned generation);

    /* Core PSO velocity + position update — unchanged from
     * original; EAPSO diversity operates only on exemplar
     * selection, not on the discrete update mechanics. */
    void PSOAlgorithm(unsigned generation);

private:
    /* Initialise variable bounds and calculate nAllele */
    void initBounds();

    /* Initialise the vi/chromX vectors of a statistics particle */
    void reserveParticleMemory(Particle& par, bool withChrom);

    /* ---- EAPSO private methods ---- */

    /* Bin all bestIndividualExp fitness values into
     * EAPSO_NUM_SEGMENTS segments, assign fitnessRank and
     * entropyRank to every particle in bestIndividualExp.
     *
     * fitnessRank: 1-based rank by descending fitness
     *   (rank 1 = highest fitness = best solution quality).
     *
     * entropyRank: 1-based rank of each particle's segment
     *   ordered by ascending occupancy count
     *   (rank 1 = least occupied segment = most novel/sparse
     *   region of the fitness landscape, i.e. best for
     *   diversity).  Particles in the same segment share
     *   the same entropyRank.
     *
     * Both ranks are written directly into the
     * bestIndividualExp[i] Particle fields. */
    void computeEntropyRanks();

    /* Compute Q = EAPSO_ALPHA * fitnessRank + entropyRank
     * for every particle in bestIndividualExp and write it
     * into the quality field.  Lower Q is better.
     * Must be called after computeEntropyRanks(). */
    void computeQuality();
};

} // namespace PSwarm