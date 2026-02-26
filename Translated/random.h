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
/*                                                      */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* File: random.h                                       */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Header file of the random.c library     */
/********************************************************/

#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

#ifndef RANDOM_LIB
#define RANDOM_LIB

typedef struct eo {
    double value;
    unsigned pos;
} SortElement;

SortElement *auxSortArray;     /* Auxiliary array for the merge sort algorithm */

/* Functions of the random library */
unsigned initRandom(unsigned);                 /* Initialization of the seed for random number generation */
double rndF(void);                             /* Generates random numbers between 0 and 1 */
unsigned rndI(unsigned);                       /* Generates a random integer between 0 and rng */
double rndFR(double, double);                  /* Generates a random real number between lower and upper bound */
unsigned rndIR(unsigned, unsigned);            /* Generates a positive random integer between lower and upper bound */
unsigned flip(double);                         /* Returns true with probability prob */
void swap(int*, int*);                         /* Swaps the values of two integers */
void swapUnsigned(unsigned*, unsigned*);       /* Swaps the values of two unsigned integers */
void shuffle(unsigned*, unsigned);             /* Shuffling of numbers */
unsigned round_custom(double);                 /* Rounds a number */
int rndInt(int);                               /* Generates an integer not greater than limit-1 */
int rndIntRange(int, int);                     /* Generates an integer between lowerLimit and upperLimit-1 */
double sigmoid(double);                        /* Sigmoid function that returns a value between 0 and 1 based on the parameter */
void algMergeSort(SortElement*, unsigned);     /* Main merge sort algorithm, reserves memory for the auxiliary array */
void mergeSort(SortElement*, unsigned, unsigned); /* Merge sort algorithm for data ordering */
void merge(SortElement*, unsigned, unsigned, unsigned); /* Function used by merge sort that “merges” data */

#endif
