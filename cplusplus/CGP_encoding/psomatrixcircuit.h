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
    double   phi1           = 0.0;  /* Cognitive acceleration (pbest) */
    double   phi2           = 0.0;  /* Social acceleration (sbest) */
    double   vMax           = 0.0;  /* Maximum velocity */
    double   pMut           = 0.0;  /* Mutation probability */
    unsigned representation = 0;    /* BINARY=2, INTEGER_A=1, INTEGER_B=0 */
    unsigned cardinality    = 0;    /* Cardinality for integer representations */

    /* ---- Parallel-Swarms-Oriented PSO (PSO-PSO) parameters ----
       Gonsalves & Egashira (2013).  Set numSubswarms = 1 to
       fall back to the original neighbourhood PSO.            */
    unsigned numSubswarms   = 1;    /* K : number of sub-swarms          */
    unsigned subswarmIters  = 0;    /* N : multi-evolutionary iters      */
    double   phi3           = 0.0;  /* Global-best acceleration (gbest)  */

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

    /* ---- PSO-PSO bookkeeping (only used when numSubswarms > 1) ---- */
    std::vector<Particle> subswarmBest;         /* sbest per sub-swarm  */
    Particle              globalBest;           /* gbest across all sbests */
    std::vector<unsigned> subswarmStart;        /* Size K+1: particle ranges */

    /* Per-sub-swarm fitness evaluator.  Each thread in the parallel
       multi-evolutionary phase uses its own MatrixDecoder instance so
       that internal scratch buffers do not race.  subDecoders has
       exactly numSubswarms entries (empty when numSubswarms == 1). */
    std::vector<MatrixDecoder> subDecoders;

    /* ---- Per-sub-swarm per-inner-iter statistics (PSO-PSO) ----
       Populated by each thread during the parallel multi-phase and
       reduced/emitted serially after the barrier.  Outer index is the
       sub-swarm k (size K), inner index is the multi-phase iteration n
       (size N).  These buffers are allocated once in reserveMemory. */
    struct InnerStat {
        double   meanSum      = 0.0;  /* sum of pbest fitness over the sub-swarm */
        double   squaredSum   = 0.0;  /* sum of pbest fitness^2 over the sub-swarm */
        Particle best;                /* best pbest inside the sub-swarm   */
        Particle worst;               /* worst pbest inside the sub-swarm  */
    };
    std::vector<std::vector<InnerStat>> innerStats;

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

    /* Compute fitness for a single particle using the given decoder.
       Passing the decoder explicitly makes this thread-safe when each
       sub-swarm has its own MatrixDecoder instance. */
    void evaluateParticle(Particle& p, MatrixDecoder& dec);

    /* Deep-copy src into dst */
    void copyParticle(const Particle& src, Particle& dst);

    /* Apply mutation to the working population */
    void mutation();

    /* Apply mutation to only the particles in sub-swarm k
       (range [subswarmStart[k], subswarmStart[k+1]) ). */
    void mutationSubswarm(unsigned k);

    /* ---- PSO-PSO island-model helpers ----
       These act on a single sub-swarm only (particle range
       [subswarmStart[k], subswarmStart[k+1]) ) and must be callable
       concurrently by different threads as long as each thread uses
       its own sub-swarm index. */

    /* Evaluate all particles in sub-swarm k, update their pbest, and
       refresh subswarmBest[k].  `statSlot` is the index into
       innerStats[k] where per-iter statistics are written.  `firstGen`
       is true only for the very first multi-phase iteration of the
       whole run and triggers unconditional pbest / sbest seeding. */
    void evaluateSubswarm(unsigned k, unsigned statSlot, bool firstGen,
                          MatrixDecoder& dec);

    /* Standard 2-term velocity/position update, restricted to the
       particles of sub-swarm k.  `firstGen` selects the velocity
       initialisation branch used on the first ever iteration. */
    void PSOAlgorithmSubswarm(unsigned k, bool firstGen);

    /* Aggregate the K sub-swarm InnerStat records at multi-iteration n
       into this->Gen (so CSV output stays identical to the original
       serial PSO). */
    void reduceInnerStats(unsigned n);

    /* Execute one complete PSO-PSO run (island model + OpenMP). */
    void pSwarmPSOPSO(unsigned runIndex);

    /* Record generation statistics and update run best/worst */
    void runInfo(const std::string& filename, unsigned generation);

    /* Core PSO velocity + position update (multi / standard phase).
       Uses bestIndividualExp[i] (pbest) and bestSocialExp[i] (sbest). */
    void PSOAlgorithm(unsigned generation);

    /* PSO-PSO single-evolutionary phase: 3-term velocity update that
       adds an attraction toward the global best (gbest) in addition
       to the sub-swarm best (sbest) and personal best (pbest).      */
    void PSOSinglePhase(unsigned generation);

private:
    /* Initialise variable bounds and calculate nAllele */
    void initBounds();

    /* Initialise the vi/chromX vectors of a statistics particle */
    void reserveParticleMemory(Particle& par, bool withChrom);
};

} // namespace PSwarm
