/********************************************************/
/*                  CINVESTAV - IPN                     */
/*        Department of Electrical Engineering          */
/*                 Computing Section                    */
/*                                                      */
/*               Evolutionary Computation               */
/*                                                      */
/*                Erika Hernandez Luna                  */
/*         eluna@computacion.cs.cinvestav.mx            */
/*                   January 21, 2003                   */
/*         Converted to C++ - April 2026                */
/*                                                      */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* File: random.h                                       */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Header for the random utility library.  */
/********************************************************/

#pragma once

#include <ctime>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>

namespace Random {

/* -------------------------------------------------------
 * Thread-local engine.
 *
 * std::rand()/std::srand() are PROCESS-GLOBAL state, so
 * calling them from multiple threads concurrently is a
 * data race (and even race-free, every thread would draw
 * from the same shared stream). Each thread that runs an
 * independent Swarm instance must get its own RNG stream,
 * seeded independently, so this is declared thread_local:
 * each OS thread gets its own private instance the first
 * time it's touched.
 * ------------------------------------------------------- */
extern thread_local std::mt19937 engine;

/*------------------------------------------------------*/
/* SortElement: used by the merge sort algorithm.       */
/*------------------------------------------------------*/
struct SortElement {
    double   value = 0.0;
    unsigned pos   = 0;
};

/* Seed initialisation */
unsigned initRandom(unsigned seed = 0);

/* Basic random number generators */
double   rndF();                            /* Uniform [0, 1] */
unsigned rndI(unsigned rng);               /* Integer in [0, rng-1] */
double   rndFR(double lower, double upper); /* Real in [lower, upper] */
unsigned rndIR(unsigned lower, unsigned upper); /* Unsigned in [lower, upper) */
unsigned flip(double prob);                /* Bernoulli trial */
int      rndInt(int limit);                /* Integer in [0, limit) */
int      rndIntRange(int lower, int upper);/* Integer in [lower, upper) */

/* Utility */
void     swapInts(int& a, int& b);
void     swapUnsigned(unsigned& a, unsigned& b);
void     shuffle(std::vector<unsigned>& arr);
unsigned roundCustom(double a);

/* Math */
double sigmoid(double val);

/* Sorting */
void algMergeSort(std::vector<SortElement>& data);

} // namespace Random

/* -------------------------------------------------------
 * Convenience aliases – bring the most-used symbols into
 * global scope so existing call sites need no changes.
 * ------------------------------------------------------- */
using Random::flip;
using Random::rndF;
using Random::rndFR;
using Random::rndI;
using Random::rndIR;
using Random::rndInt;
using Random::rndIntRange;
using Random::sigmoid;
using Random::initRandom;
using Random::shuffle;
