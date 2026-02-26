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
/* File: matrixpso.c                                    */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Library for handling the matrix that     */
/* represents a circuit.                                */
/********************************************************/

#include <math.h>
#include <string.h>
#include "matrixpso.h"
#include "circuits.h"
#include "psomatrixcircuit.h"

/* Reserve memory for global variables used in decoding */
void reserveMatrixMemory(void)
{
	unsigned i;

    input1   = memVarX;
    input2   = memVarX;
    gateType = memVarX;
    output   = (unsigned *)malloc(sizeof(unsigned) * numRows);
    gateCount= (unsigned *)malloc(sizeof(unsigned) * tMat);
    inTT     = (unsigned **)malloc(sizeof(unsigned *) * numRowsTT);

    for(i = 0; i < numRowsTT; i++)
        inTT[i] = (unsigned *)malloc(sizeof(unsigned) * numRows);
}

/* Free memory used by global variables in decoding */
void freeMatrixMemory(void)
{
	unsigned i;

	free(input1);
	free(input2);
	free(gateType);
	free(output);
	free(gateCount);
	
    for(i = 0; i < numRowsTT; i++)
        inTT[i] = (unsigned *)malloc(sizeof(unsigned) * numRows); // <- Bug? Should free, not malloc
    free(inTT);
}

/* Decode the matrix into the values of the variables */
void decode(matrix M)
{
	int k;
	unsigned i,j,num,in;

	for(i=0,j=0,in=0; i<nVar; i++){

        if(representation != BINARY)
			num = M[i];
		else {
            num = 0;
            for (k = bitVariable[i]-1; k >= 0; k--, j++)
                if(M[j]) num += (unsigned)pow(2, k);
        }

		switch(i % 3) {
            case 0: /* Input 1 */
                input1[in] = (num % numRows);
                break;
            case 1: /* Input 2 */
                input2[in] = (num % numRows);
                break;
            case 2: /* Gate type */
                gateType[in++] = (num % numGates);
        }
	}
}

/* Evaluate the matrix: obtain number of matches with target circuit and count WIREs */
void evaluate(matrix M, unsigned *numEqual)
{
	unsigned i,j,k,in;

    *numEqual = 0;
    decode(M);

    /* Initialize truth table inputs */
	for(i = 0; i < numRowsTT; i++)
        for(j = 0; j < numRows; j++)
            inTT[i][j] = j < numInputs ? inputTT[i][j] : inTT[i][j - numInputs];	

    /* Evaluate for each row in truth table */
	for(i = 0; i < numRowsTT; i++) {
        for(j = 0, in = 0; j < tMat; j++) {
            switch(gateType[j]) {
                case AND : 
                    output[in] = inTT[i][input1[j]] & inTT[i][input2[j]];
                    break;
                case OR : 
                    output[in] = inTT[i][input1[j]] | inTT[i][input2[j]];
                    break;
                case NOT : 
                case NOT1:
                    output[in] = inTT[i][input1[j]] ? 0 : 1;
                    break;
                case WIRE :
                case WIRE1:
                    output[in] = inTT[i][input1[j]];
                    break;
                case XOR : 
                case XOR1 : 
                    output[in] = inTT[i][input1[j]] ^ inTT[i][input2[j]];
                    break;
            }

            if(!((j + 1) % numRows)) {
                in = 0;
                for (k = 0; k < numRows; k++) 
                    inTT[i][k] = output[k];
            }
            else in++;
        }

		/* Compare outputs with target truth table */
        for(j = 0; j < numOutputs; j++)
            if(output[j] == outputTT[i][j]) (*numEqual)++; 
    }
}

/* Count the total number of gates used in the solution */
unsigned countGates(void)
{
    unsigned i, gates = 0; 
    unsigned startOutput = numRows * (numCols - 1);

    for(i = 0; i < tMat; i++) gateCount[i] = 0;
    for(i = startOutput; i < startOutput + numOutputs; i++)
        countGate(i, numCols - 1, gateCount);
    for(i = 0; i < tMat; i++) if(gateCount[i]) gates++;

    return gates;
}

/* Count the gates involved in a solution (recursive) */
void countGate(unsigned out, int col, unsigned *gateCount)
{
    if(col >= 0 && out >= 0) {
        if(gateType[out] != WIRE && gateType[out] != WIRE1)
            if(gateType[out] != AND && gateType[out] != OR)
                gateCount[out] = 1;
            else
                if(input1[out] != input2[out])
                    gateCount[out] = 1;

        countGate(numRows * (col - 1) + input1[out], col - 1, gateCount);
        countGate(numRows * (col - 1) + input2[out], col - 1, gateCount);  
    }
}

/* Get the Boolean expression that represents the matrix */
void expression(matrix M, boolExpr expBool, unsigned out)
{
    strcpy(expBool, "");
    decode(M);
    booleanString(numRows * (numCols - 1) + out, numCols - 1, expBool);
}

/* Recursive builder of the Boolean expression string */
void booleanString(unsigned cell, int col, boolExpr str)
{
    unsigned new1, new2;
    char gateStr[7], varStr[2];

    varStr[1] = EndString;
    switch(gateType[cell]) {
			case AND:
				strcpy(gateStr,"(AND1 ");
				break;
			case OR:
				strcpy(gateStr,"(OR1 ");
				break;
			case NOT:
			case NOT1:
				strcpy(gateStr,"(NOT1 ");
				break;
			case XOR:
			case XOR1:
				strcpy(gateStr,"(XOR1 ");
				break;
			case WIRE:
			case WIRE1:
				strcpy(gateStr,"(WIRE ");
				break;
	}

    strcat(str, gateStr);
	if(col<=0){
        varStr[0] = input1[cell] % numInputs + (int)'A';
        strcat(str, varStr);
		if(gateType[cell] != WIRE && gateType[cell] != WIRE1 && gateType[cell] != NOT  && gateType[cell] != NOT){ // change last NOT1 to NOT
            varStr[0] = input2[cell] % numInputs + (int)'A';
            strcat(str, " ");
            strcat(str, varStr);
		}
	}
	else{
        new1 = numRows * (col - 1) + input1[cell];
        new2 = numRows * (col - 1) + input2[cell];

        booleanString(new1, col - 1, str);
		if(gateType[cell] != WIRE && gateType[cell] != WIRE1 && gateType[cell] != NOT  && gateType[cell] != NOT){ // change last NOT1 to NOT
            strcat(str, " ");
            booleanString(new2, col - 1, str);
		}
	}
    strcat(str, ")");
}

/* Print the matrix of gates */
void printMatrix(matrix M)
{
	int dec,sign;
	unsigned i, j;
    char str[5000];
	
	printf("\n");
    decode(M);
    for(i = 0; i < numRows; i++) {
        strcpy(str, "");
        for(j = 0; j < numCols; j++) {
            switch(gateType[j * numRows + i]) {
					case AND : 
						strcat(str,"AND(");
						break;
					case OR: 
						strcat(str,"OR(");
						break;
					case NOT : 
						strcat(str,"NOT(");
						break;
					case WIRE : 
						strcat(str,"WIRE(");
						break;
					case XOR : 
						strcat(str,"XOR(");
						break;
					case NOT1:
						strcat(str,"NOT(");
						break;
					case WIRE1:
						strcat(str,"WIRE(");
						break;
					case XOR1 : 
						strcat(str,"XOR(");
						break;

			}

            strcat(str, input1[j*numRows + i] == 0 ? "0" : fcvt(input1[j*numRows + i],0,&dec,&sign));
			strcat(str, " ");
            strcat(str, input2[j*numRows + i] == 0 ? "0" : fcvt(input2[j*numRows + i],0,&dec,&sign));
			strcat(str,")");
            if(j != numCols-1) strcat(str,",");
		}
		printf("%s\n",str);
	}
}
