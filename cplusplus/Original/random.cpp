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
/* File: random.cpp                                     */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Library for generating random numbers   */
/* with specific properties, and for random shuffling.  */
/********************************************************/

#include "random.h"
#include <algorithm>
#include <random>
#include <chrono>
 
namespace Random {

/* -------------------------------------------------------
 * Thread local Random Number Generator (RNG) instance
 * ------------------------------------------------------- */
thread_local std::mt19937 rng;

/* -------------------------------------------------------
 * Auxiliary array for merge sort.
 * thread_local prevents multiple PSO threads from modifying
 * the same temporary array.
 * ------------------------------------------------------- */
thread_local std::vector<SortElement> auxSortArray;

/* -------------------------------------------------------
 * Seed initialisation
 * ------------------------------------------------------- */
unsigned initRandom(unsigned seed)
{
    // unsigned realSeed = seed ? seed : static_cast<unsigned>(std::time(nullptr));
    // std::srand(realSeed);
    // return realSeed;

    unsigned realSeed;

    if (seed != 0) {
        realSeed = seed;
    } else {
        /*
         * Generate a seed using the current high-resolution clock.
         * Each thread will seed it's own RNG.
         */
        realSeed = static_cast<unsigned>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    }

    rng.seed(realSeed);

    return realSeed;
}

/* -------------------------------------------------------
 * Basic random number generators
 * ------------------------------------------------------- */

double rndF()
{
    // return static_cast<double>(std::rand()) / RAND_MAX;
    /*
     * Generate a floating-point number in [0, 1)
     */
    return std::generate_canonical<double, 53>(rng);
}

unsigned rndI(unsigned rngLimit)
{
    // double val = rndF() * static_cast<double>(rng) - 1.0;
    // return val < 0.0 ? 0u : static_cast<unsigned>(val);
    
    double val = rndF() * static_cast<double>(rngLimit) - 1.0;

    return val < 0.0 ? 0u : static_cast<unsigned>(val);
}

double rndFR(double lower, double upper)
{
    return lower + rndF() * (upper - lower);
}

unsigned rndIR(unsigned lower, unsigned upper)
{
    return static_cast<unsigned>(lower + rndF() * static_cast<double>(upper - lower));
}

unsigned flip(double prob)
{
    return rndF() <= prob ? 1u : 0u;
}

int rndInt(int limit)
{
    // return static_cast<int>(static_cast<double>(std::rand()) / RAND_MAX * limit);
    
    return static_cast<int>(rndF() * static_cast<double>(limit));
}

int rndIntRange(int lower, int upper)
{
    if (lower > upper) std::swap(lower, upper);
    return rndInt(upper - lower) + lower;
}

/* -------------------------------------------------------
 * Utility
 * ------------------------------------------------------- */

void swapInts(int& a, int& b)
{
    std::swap(a, b);
}

void swapUnsigned(unsigned& a, unsigned& b)
{
    std::swap(a, b);
}

void shuffle(std::vector<unsigned>& arr)
{
    for (int i = static_cast<int>(arr.size()) - 1; i > 0; --i)
        swapUnsigned(arr[static_cast<unsigned>(i)], arr[rndI(static_cast<unsigned>(i))]);
}

unsigned roundCustom(double a)
{
    double w = std::fabs(a) - std::floor(std::fabs(a));
    w = (w < 0.5) ? std::floor(std::fabs(a)) : std::ceil(std::fabs(a));
    w = (a < 0.0) ? w * -1.0 : w;
    return static_cast<unsigned>(static_cast<int>(w));
}

/* -------------------------------------------------------
 * Math
 * ------------------------------------------------------- */

double sigmoid(double val)
{
    return 1.0 / (1.0 + std::exp(-val));
}

/* -------------------------------------------------------
 * Merge sort (internal helpers are file-local)
 * ------------------------------------------------------- */

static void merge(std::vector<SortElement>& data,
                  unsigned start, unsigned mid, unsigned end)
{
    unsigned j = start, k = mid + 1;

    for (unsigned i = start; i <= end; ++i) {
        if ((data[j].value >= data[k].value && j <= mid) || k > end) {
            auxSortArray[i] = data[j++];
        } else {
            auxSortArray[i] = data[k++];
        }
    }
    for (unsigned i = start; i <= end; ++i)
        data[i] = auxSortArray[i];
}

static void mergeSort(std::vector<SortElement>& data,
                      unsigned start, unsigned end)
{
    if (start < end) {
        unsigned mid = (start + end) / 2;
        mergeSort(data, start, mid);
        mergeSort(data, mid + 1, end);
        merge(data, start, mid, end);
    }
}

void algMergeSort(std::vector<SortElement>& data)
{
    if (data.empty()) return;
    auxSortArray.resize(data.size());
    mergeSort(data, 0, static_cast<unsigned>(data.size() - 1));
    auxSortArray.clear();
}

} // namespace Random
