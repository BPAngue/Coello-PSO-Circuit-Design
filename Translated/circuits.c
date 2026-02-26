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
/* File: circuits.c                                     */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Library to load into memory a file with */
/* the design of the logic circuit to optimize          */
/********************************************************/

#include "circuits.h"

/* Generates the inputs of the circuit's truth table */
void generateTT(void)
{
	unsigned i, j;

	for(i=0; i<numRowsTT; i++){
        inputTT[i] = (unsigned *)malloc(sizeof(unsigned)*numInputs);
        for(j=0; j<numInputs; j++)
            inputTT[i][j] = (i>>j)&1;
    }
}

/* Loads the truth table data from the input file into memory */
void loadTT(FILE *pf)
{
	int j;
	unsigned i;

	numRowsTT = 1;
    numInputs = (unsigned)readNumber(pf);
    numRowsTT <<= numInputs;
    numOutputs = (unsigned)readNumber(pf);
    numTotalOutputs = numOutputs * numRowsTT;

	inputTT  = (unsigned **)malloc(sizeof(unsigned *)*numRowsTT);
    outputTT = (unsigned **)malloc(sizeof(unsigned *)*numRowsTT);
	
	generateTT();

	for(i=0; i<numRowsTT; i++){
        outputTT[i] = (unsigned *)malloc(sizeof(unsigned)*numOutputs);
        for(j=0; (unsigned)j<numOutputs; j++)
            outputTT[i][j] = (unsigned)readNumber(pf);
    }
}

/* Displays the truth table data obtained from the input file */
void printTT(void)
{
	int j;
	unsigned i;

	printf("      Number of inputs in the truth table: %d\n",numInputs);
    printf("     Number of outputs in the truth table: %d\n",numOutputs);
    printf("       Number of rows in the truth table: %d\n",numRowsTT);


	printf("\nTruth Table\n\n");
    for(i=0; i<numInputs; i++)  printf("E%d ",i);
    for(i=0; i<numOutputs; i++) printf("S%d ",i);
	
	for(i=0; i<numRowsTT; i++){
        printf("\n");
        for(j=numInputs-1; j>=0; j--)
            printf(" %d ",inputTT[i][j]);
        for(j=0; (unsigned)j<numOutputs; j++)
            if(outputTT[i][j] == DONTCARE)
                printf(" * ");
            else    
                printf(" %d ",outputTT[i][j]);
    }
	printf("\n\n");
}

/* Frees the memory reserved for the circuit’s truth table */
void freeMemoryTT(void)
{
	unsigned i;

	for(i=0; i<numRowsTT; i++){
        free(inputTT[i]);
        free(outputTT[i]);
    }
    free(inputTT);
    free(outputTT);
}

/* Reads a number from the file pointed to by pf */
double readNumber(FILE *pf)
{
	char car, num[100];
	unsigned i = 0;

    while(!feof(pf) && !isDigit(car = fgetc(pf)) && car != '.' && car != ';');
	if(car == ';'){ 
		while(car != EndLine && !feof(pf)) car = fgetc(pf);
		car = fgetc(pf);
	}

	while(!feof(pf) && (car != EndLine && car!=' ' && car!=';')){
		if(isDigit(car) || car=='.')num[i++] = car;
		car = fgetc(pf);
	}
	if(car == ';') while((car = fgetc(pf)) != EndLine && !feof(pf));
	
	num[i] = EndString;
	return atof(num);
}

/* Reads a string from the file pointed to by pf */
void readString(FILE *pf, char *str)
{
	unsigned i = 0;
	char car;

	while((car = fgetc(pf)) != ';' && car!=EndLine && !isalpha(car));
	if(car == ';'){
		while(car != EndLine && !feof(pf)) car = fgetc(pf);
		car = fgetc(pf);
	}
	
	while(!feof(pf) && (car != EndLine && car!=' ' && car!=';' && car!='\t')){
		str[i++] = car;
		car = fgetc(pf);
	}
	if(car == ';') while((car = fgetc(pf)) != EndLine && !feof(pf));
	str[i] = EndString;
}

/* Returns 1 if the character is between 0–9, 0 otherwise */
unsigned isDigit(char a)
{
	return a=='0' || a=='1' || a=='2' || a=='3' || a=='4' || a=='5' || a=='6' || a== '7' || a=='8' || a=='9' ? 1 : 0;
}
