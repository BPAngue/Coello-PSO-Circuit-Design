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

/* -------------------------------------------------------
 * Auxiliary array for merge sort – kept file-local so it
 * replaces the old global SortElement *auxSortArray.
 * ------------------------------------------------------- */
static std::vector<SortElement> auxSortArray;

/* -------------------------------------------------------
 * Seed initialisation
 * ------------------------------------------------------- */
unsigned initRandom(unsigned seed)
{
    unsigned realSeed = seed ? seed : static_cast<unsigned>(std::time(nullptr));
    std::srand(realSeed);
    return realSeed;
}

/* -------------------------------------------------------
 * Basic random number generators
 * ------------------------------------------------------- */

double rndF()
{
    return static_cast<double>(std::rand()) / RAND_MAX;
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
    return static_cast<int>(static_cast<double>(std::rand()) / RAND_MAX * limit);
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
