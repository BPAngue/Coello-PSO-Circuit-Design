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

namespace Random {

/* Definition of the thread-local engine declared in random.h.
   Each thread gets its own instance, default-constructed with
   a fixed seed until initRandom() reseeds it for that thread. */
thread_local std::mt19937 engine{ std::random_device{}() };

/* -------------------------------------------------------
 * Auxiliary array for merge sort.
 *
 * Was a single shared `static` (file-local, but still one
 * instance for the whole process) — two threads sorting at
 * the same time would corrupt each other's data. Made
 * thread_local so each thread owns its own scratch buffer.
 * ------------------------------------------------------- */
static thread_local std::vector<SortElement> auxSortArray;

/* -------------------------------------------------------
 * Seed initialisation
 *
 * Reseeds THIS THREAD's engine only. Callers running one
 * Swarm per thread should call initRandom() once at the
 * start of that thread (e.g. with a per-thread/per-run seed
 * such as std::random_device{}() ^ threadIndex, or 0 to let
 * it fall back to a time-derived seed) so different threads
 * don't end up correlated.
 * ------------------------------------------------------- */
unsigned initRandom(unsigned seed)
{
    unsigned realSeed = seed ? seed : static_cast<unsigned>(std::time(nullptr));
    engine.seed(realSeed);
    return realSeed;
}

/* -------------------------------------------------------
 * Basic random number generators
 * ------------------------------------------------------- */

double rndF()
{
    static thread_local std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(engine);
}

unsigned rndI(unsigned rng)
{
    double val = rndF() * static_cast<double>(rng) - 1.0;
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
    if (limit <= 0) return 0;
    std::uniform_int_distribution<int> dist(0, limit - 1);
    return dist(engine);
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
