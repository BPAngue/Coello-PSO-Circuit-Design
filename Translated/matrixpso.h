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
/* File: matrixpso.h                                    */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Header file for the library matrixpso.c */
/********************************************************/

#include <malloc.h>

#ifndef MATRIXPSO
#define MATRIXPSO

/* Define the matrix data type */
#define matrix unsigned*    

/* Reserve memory for decoding values of each cell in the matrix */
#define memValCell (unsigned)malloc(sizeof(unsigned)*numRows*numCols)    

/* Maximum size of Boolean expression and its type */
#define maxBoolExprSize 4000            
typedef char boolExpr[maxBoolExprSize];

unsigned numRows;       /* Number of rows in the matrix */
unsigned numCols;       /* Number of columns in the matrix */
unsigned *input1;       /* Decoded values of gate input 1 (temporary) */
unsigned *input2;       /* Decoded values of gate input 2 (temporary) */
unsigned *gateType;     /* Decoded values of gate type (temporary) */
unsigned **inTT;        /* Stores input/output of each evaluation stage (temporary) */
unsigned *gateCount;    /* Counts the gates in a matrix (temporary) */
unsigned *output;       /* Stores the output of an evaluation stage (temporary) */
unsigned tMat;          /* Size of matrix = numRows * numCols */

/* Matrix handling functions */
void reserveMatrixMemory(void);               /* Reserve memory for decoding variables */
void freeMatrixMemory(void);                  /* Free memory for decoding variables */
void decode(matrix);                          /* Convert matrix cells into input1, input2, and gateType */
void evaluate(matrix,unsigned*);              /* Evaluate the matrix, return violations and WIRE count */
unsigned countGates(void);                    /* Count total gates in the solution */
void countGate(unsigned,int,unsigned*);       /* Count gates involved in a solution */
void expression(matrix,boolExpr,unsigned);    /* Get Boolean expression of the matrix */
void booleanString(unsigned,int,boolExpr);    /* Recursive builder of Boolean expression string */
void printMatrix(matrix M);                   /* Print the matrix of gates */

#endif
