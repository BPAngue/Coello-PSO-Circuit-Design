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
/*                                                      */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* File: statistics.h                                   */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Header file for the library             */
/* statistics.c.                                        */
/********************************************************/

#include "psomatrixcircuit.h"

#ifndef STATISTICS
#define STATISTICS

#define sizeChrom 10000        /* Max characters in chromosome string */
#define sizeNumDouble 100      /* Max digits for a double */    

/* Structure to store statistics of a generation or run */
typedef struct{                
    double meanFitness;        /* Mean fitness of particles */
    double squaredFitness;     /* Sum of squared fitness values (for variance) */
    unsigned generation;       /* Generation when optimum found */
    particle best;             /* Best particle */
    particle worst;            /* Worst particle */
} statistics;

statistics Gen;                /* Generation statistics */
statistics Run;                /* Run statistics */

/* Function prototypes */
void initStatistics(statistics*);           /* Initializes the information of the run or generation statistics */
void globalHeader(char*);                   /* Generates the header for the global statistics of all runs */             
void runStatistics(char*,unsigned);         /* Generates the statistics of each run in the global runs file */  
void runHeader(char*,unsigned);             /* Generates the header for the run statistics */        
void generationStats(char*,unsigned);       /* Generates the statistics of one generation of the run */  
void runFooter(char*);                      /* Generates the final results of the run statistics */             
void chromToString(particle,char*);         /* Converts the chromosome of real numbers or binary to a character string */    
void dtoa(double,unsigned,char*);           /* Converts a real number into a character string */      

#endif
