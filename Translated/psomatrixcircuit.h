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
/* File: psomatrixcircuit.h                             */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Header file for the library             */
/* psomatrixcircuit.c.                                  */ 
/********************************************************/

#include <malloc.h>
#include <stdlib.h>

#ifndef PSWARM
#define PSWARM

/* Macros for memory allocation */
#define memPopulation   (particle *)malloc(sizeof(particle)*tPop)      /* Population memory */
#define memVarX         (unsigned *)malloc(sizeof(unsigned)*nVar)      /* Variable memory */
#define memChromX       (unsigned *)malloc(sizeof(unsigned)*nAllele)   /* Chromosome memory */
#define memVi           (double *)malloc(sizeof(double)*nAllele)       /* Velocity memory */

#define sizeFileName 200   /* Max characters in filenames */
#define BINARY   2         /* Binary representation */
#define INTEGERA 1         /* Integer representation A */
#define INTEGERB 0         /* Integer representation B */

/* Structure of a population particle */
typedef struct{                    
    /*** PSO Information ***/
    double *vi;           /* Velocities for each allele */
    unsigned *chromX;     /* Chromosome that evolves in PSO */

    /*** Individual Information ***/
    unsigned numGates;    /* Number of gates used in solution */
    unsigned numNoGates;  /* Number of wires in solution */
    unsigned numEqual;    /* Number of matching outputs between circuit truth table and target */
    double fitness;       /* Particle fitness */
} particle;

/* File name variables */
char nfInput[sizeFileName];   /* Input data file name */
char nfRun[sizeFileName];     /* Run statistics file name */
char nfGen[sizeFileName];     /* Generation statistics file name */

/* User-specified variables */
unsigned tPop;                /* Population size */
unsigned nGen;                /* Max generations */
unsigned nRun;                /* Max runs */
unsigned tNeigh;              /* Neighborhood size */
double phi1;                  /* Acceleration constant 1 */
double phi2;                  /* Acceleration constant 2 */
double vMax;                  /* Max velocity */
double pMut;                  /* Mutation probability */
unsigned representation;      /* Representation Binary=2, IntA=1, IntB=0 */
unsigned cardinality;

/* Program variables */
unsigned nVar;                /* Number of problem variables */
unsigned nAllele;             /* Number of alleles in chromosome */    
unsigned *lInf;               /* Lower bounds */
unsigned *lSup;               /* Upper bounds */
unsigned *bitVariable;        /* Bits per variable (binary encoding) */

/* Population variables */
particle *population;              /* Structure that stores a population */
particle *bestSocialExp;           /* Structure that stores the population of the best social experiences */      
particle *bestIndividualExp;       /* Structure that stores the population of the best individual experiences */ 

/* Function prototypes */
unsigned loadParameters(char*);                     /* Stores in memory the data obtained from the input file */             
void printParameters(void);                         /* Prints on screen the program input parameters received from the file */               
void initVariables(void);                           /* Initializes the variables required for the program */                  
void initBounds(void);                              /* Initializes the ranges of the variables of the particle’s chromosome */                     
void reserveMemory(void);                           /* Reserves the memory required for the program variables */                  
void reserveParticleMemory(particle*,unsigned);     /* Reserves memory for the elements of the particle */ 
void freeMemory(void);                              /* Frees all memory reserved throughout the program once the runs are completed */                      
void freeParticleMemory(particle*,unsigned);        /* Frees the memory reserved for each particle */
void pSwarm(unsigned);                              /* Main program of the Particle Swarm */                    
void runFileName(unsigned,char*);                   /* Generates the filename for the statistics of the current run */           
void initPopulation(void);                          /* Generates the initial random population of individuals */                  
void evaluatePopulation(unsigned);                  /* Calculates the characteristics of each particle */         
void evaluateParticle(particle*);                   /* Obtains the fitness of the particle */           
void copyParticle(particle,particle*);              /* Copies the data from particle 'in' to 'out' */      
void mutation(void);                        
void runInfo(char*,unsigned);                      /* Computes the data for the run statistics */
void PSOAlgorithm(unsigned);                        /* Particle Swarm Optimization algorithm */               

#endif
