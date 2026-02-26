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
/* File: random.c                                       */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Library for generating random numbers   */
/* with specific properties, and for random shuffling   */
/* of numbers.                                          */
/********************************************************/

#include "random.h"

/* Initialization of the seed for random number generation */
unsigned initRandom(unsigned seed)
{    
    unsigned realSeed = !seed ? (unsigned)(time(NULL)) : seed;
	srand(realSeed);
	return realSeed;
}

/* Generates random numbers between 0 and 1 */
double rndF(void)
{
	return ((double)rand()) / RAND_MAX;
}

/* Generates a random integer between 0 and rng */
unsigned rndI(unsigned rng)
{
	double val;
	if((val = rndF() * (rng) - 1) < 0) return 0;
	
    return (unsigned)(val);
}

/* Generates a random real number between lower and upper bound */
double rndFR(double lower, double upper)
{
    return lower + rndF() * (upper - lower);
}

/* Generates a positive random integer between lower and upper bound */
unsigned rndIR(unsigned lower, unsigned upper)
{
    return (unsigned)(lower + rndF() * (upper - lower));
}


/* Returns true with probability prob */
unsigned flip(double prob)
{
    if(rndF() <= prob) return(1);
    return(0);
}

/* Swaps the values of two integers */
void swap(int *a, int *b)
{
	int temp = *a;

	*a = *b;
	*b = temp;
}

/* Swaps the values of two unsigned integers */
void swapUnsigned(unsigned *a, unsigned *b)
{
	int temp = *a;

	*a = *b;
	*b = temp;
}

/* Shuffling of numbers */
void shuffle(unsigned *arr, unsigned size)
{
	int i;
	
    for(i = size - 1; i > 0; i--)
        swapUnsigned(&(arr[i]), &(arr[rndI(i)]));
}

/* Rounding of numbers */
unsigned round_custom(double a)
{
	double w = fabs(a) - floor(fabs(a));
	w = w<0.5 ? floor(fabs(a)) : ceil(fabs(a));
	w = a<0 ? w*(-1) : w;
	return((int)w);
}

/* Generates an integer not greater than limit-1 */
int rndInt(int limit)
{	
	double div = ((double)rand() / RAND_MAX);
	return (int)(div * limit);
}

/* Generates an integer between lowerLimit and upperLimit-1 */
int rndIntRange(int lowerLimit, int upperLimit)
{
    if(lowerLimit > upperLimit) swap(&lowerLimit, &upperLimit);
    return rndInt(upperLimit - lowerLimit) + lowerLimit;
}


/* Sigmoid function that returns a value between 0 and 1 based on the parameter */
double sigmoid(double val)
{
    return 1 / (1 + exp(-val));
}

/* Main merge sort algorithm, reserves memory for the auxiliary array */
void algMergeSort(SortElement *data, unsigned size)
{
    auxSortArray = (SortElement *)malloc(sizeof(SortElement) * size);
    mergeSort(data, 0, size - 1);
    free(auxSortArray);
}

/* Merge sort algorithm for data ordering */
void mergeSort(SortElement *data, unsigned start, unsigned end)
{
    unsigned mid;

    if(start < end){
        mid = (end + start) / 2;
        mergeSort(data, start, mid);
        mergeSort(data, mid + 1, end);
        merge(data, start, mid, end);
	}
}

/* Function used by merge sort that “merges” data */
void merge(SortElement *data, unsigned start, unsigned mid, unsigned end)
{
    unsigned i, j = start, k = mid + 1;

    for(i = start; i <= end; i++){
		if((data[j].value >= data[k].value && j <= mid) || k > end){
            auxSortArray[i].value = data[j].value;
            auxSortArray[i].pos   = data[j++].pos;
        }
        else{
            auxSortArray[i].value = data[k].value;
            auxSortArray[i].pos   = data[k++].pos;
        }
	}
			
	for(i = start; i <= end; i++){
        data[i].value = auxSortArray[i].value;
        data[i].pos   = auxSortArray[i].pos;
    }
}