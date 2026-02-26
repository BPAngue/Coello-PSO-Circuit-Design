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
/* File: circuits.h                                     */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Header file of the library circuits.c   */
/********************************************************/

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#ifndef CIRCUITS
#define CIRCUITS

#define EndLine   '\n'
#define EndString '\x0'

/* Functions that may be contained in the tree */
#define AND   0
#define OR    1
#define NOT   2
#define WIRE  3
#define XOR   4
#define NOT1  5
#define WIRE1 6
#define XOR1  7

#define DONTCARE 2   /* "Don't care" condition in the output truth table */

unsigned numGates;       /* Number of available gates */
unsigned numInputs;      /* Number of inputs in the truth table */
unsigned numOutputs;     /* Number of outputs defined in the truth table */
unsigned numRowsTT;      /* Number of rows in the truth table */
unsigned numTotalOutputs;/* Total number of outputs the circuit must satisfy, numOutputs * numRowsTT */
unsigned **inputTT;      /* Inputs in the truth table */
unsigned **outputTT;     /* Outputs in the truth table */

void generateTT(void);                /* Generates the inputs of the circuit's truth table */
void loadTT(FILE *);                  /* Loads the truth table data from the input file into memory */
void printTT(void);                   /* Displays the truth table data obtained from the input file */
void freeMemoryTT(void);              /* Frees the memory reserved for the circuit’s truth table */
double readNumber(FILE *);            /* Reads a number from the file pointed to by pf */
void readString(FILE *pf,char*);      /* Reads a string from the file pointed to by pf */
unsigned isDigit(char);               /* Returns 1 if the character is 0–9, 0 otherwise */

#endif
