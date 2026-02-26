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
/* File: statistics.c                                   */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Library to generate statistics per       */
/* generation and per run of the Particle Swarm          */
/* Optimization. Statistics include the best and worst   */
/* particle, mean fitness, best and worst fitness,       */
/* standard deviation, variance, and the chromosome      */
/* string of the best particle found in each generation. */
/********************************************************/

#include "statistics.h"
#include "psomatrixcircuit.h"
#include "circuits.h"
#include "matrixpso.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* Initializes the information of the run or generation statistics */
void initStatistics(statistics *st)
{
    st->squaredFitness = 0;
    st->meanFitness = 0;
}

/* Generates the header for the global statistics of all runs */
void globalHeader(char *filename)
{
    FILE *pf = fopen(filename,"w");

	fprintf(pf,"Centro de Investigación y de Estudios Avanzados del IPN\n");
    fprintf(pf,"Design of Combinational Logic Circuits using Particle Swarm Optimization\n");
    fprintf(pf,"Department of Electrical Engineering - Computing Section\n");
    fprintf(pf,"Erika Hernández Luna  eluna@computacion.cs.cinvestav.mx\n\n");
    fprintf(pf,"Statistics of the results obtained in all program runs\n\n");

    fprintf(pf,"     Processed File:,   %s\n",nfInput);
    fprintf(pf,"   Population Size:,   %u\n",tPop);
    fprintf(pf,"Number of Generations:,   %u\n",nGen);
	fprintf(pf,"                 Phi 1:, %.3f\n",phi1);
	fprintf(pf,"                 Phi 2:, %.3f\n",phi2);
	fprintf(pf,"                  Vmax:, %.3f\n",vMax);
    fprintf(pf,"Mutation Percentage:, %.3f\n",pMut*100);
    fprintf(pf,"  Number of Gates:,   %u\n",numGates);
    fprintf(pf,"        Cardinality:,   %u\n",cardinality);
    fprintf(pf,"Rows in the Matrix:,   %u\n",numRows);
    fprintf(pf,"Columns in the Matrix:,   %u\n",numCols);    
    fprintf(pf,"       Neighborhood:,   %u\n",tNeigh);
	switch(representation){
		case BINARY:
            fprintf(pf,"   Representation:, Binary\n");
			break;
		case INTEGERA:
            fprintf(pf,"   Representation:, Integer A\n");
			break;
        case INTEGERB:
            fprintf(pf,"   Representation:, Integer B\n");
	}

	fprintf(pf,"\n");
    fprintf(pf,"Run,Mean Fitness,Generation,");
    fprintf(pf,"Best Chromosome,Fitness,");
    fprintf(pf,"Worst Chromosome,Fitness,");
    fprintf(pf,"Standard Deviation,Variance, Expression\n");
    fclose(pf);
}

/* Generates the statistics of each run in the global runs file */
void runStatistics(char *filename, unsigned run)
{
    unsigned i;
    double variance;
    char best[sizeChrom], worst[sizeChrom], boolExpr[maxBoolExprSize];
    FILE *pf = fopen(filename,"a+");

    variance = Run.squaredFitness - Run.meanFitness * Run.meanFitness;
    chromToString(Run.best, best);
    chromToString(Run.worst, worst);

	fprintf(pf,"%u,%.8f,%u,%s,%.8f,%s,%.8f,%.8f,%.12f",
		run,
        Run.meanFitness,
        Run.generation,
		best,
        Run.best.fitness,
		worst,
        Run.worst.fitness,
        sqrt(variance),
		variance);

    for(i=0; i<numOutputs; i++){
        expression(Run.best.chromX, boolExpr, i);
        fprintf(pf,",%s",boolExpr);
	}
	fprintf(pf,"\n");
	fclose(pf);

}

/* Generates the header for the run statistics */
void runHeader(char *filename, unsigned seed)
{
    FILE *pf = fopen(filename,"w");

	fprintf(pf,"Centro de Investigación y de Estudios Avanzados del IPN\n");
    fprintf(pf,"Design of Combinational Logic Circuits using Particle Swarm Optimization\n");
    fprintf(pf,"Department of Electrical Engineering - Computing Section\n");
    fprintf(pf,"Erika Hernández Luna  eluna@computacion.cs.cinvestav.mx\n\n");
    fprintf(pf,"Results obtained in a program run\n\n");
	
    fprintf(pf,"     Processed File:,   %s\n",nfInput);
    fprintf(pf,"   Population Size:,   %u\n",tPop);
    fprintf(pf,"Number of Generations:,   %u\n",nGen);
	fprintf(pf,"                 Phi 1:, %.3f\n",phi1);
	fprintf(pf,"                 Phi 2:, %.3f\n",phi2);
	fprintf(pf,"                  Vmax:, %.3f\n",vMax);
    fprintf(pf,"Mutation Percentage:, %.3f\n",pMut*100);
    fprintf(pf,"  Number of Gates:,   %u\n",numGates);
    fprintf(pf,"        Cardinality:,   %u\n",cardinality);
    fprintf(pf,"Rows in the Matrix:,   %u\n",numRows);
    fprintf(pf,"Columns in the Matrix:,   %u\n",numCols);    
    fprintf(pf,"       Neighborhood:,   %u\n",tNeigh);
    fprintf(pf,"             Seed:,   %u\n",seed);
	switch(representation){
		case BINARY:
            fprintf(pf,"   Representation:, Binary\n");
			break;
        case INTEGERA:
            fprintf(pf,"   Representation:, Integer A\n");
			break;
		case INTEGERB:
            fprintf(pf,"   Representation:, Integer B\n");
	}

	fprintf(pf,"\n");
    fprintf(pf,"Generation,Mean Fitness,Best Fitness,Worst Fitness,");
    fprintf(pf,"Best Chromosome");
    fprintf(pf,"\n");
	fclose(pf);
}

/* Generates the statistics of one generation of the run */
void generationStats(char *filename, unsigned gen)
{
    char best[sizeChrom];
    FILE *pf = fopen(filename,"a+");

    chromToString(Gen.best,best);
	fprintf(pf,"%u,%.8f,%.8f,%.8f,%s,",
		gen,
        Gen.meanFitness,
        Gen.best.fitness,
        Gen.worst.fitness,
		best);
	fprintf(pf,"\n");
	fclose(pf);
}

/* Generates the final results of the run statistics */
void runFooter(char *filename)
{
	unsigned i;
    char best[sizeChrom];
    boolExpr boolExpression;
    FILE *pf = fopen(filename,"a+");

    chromToString(Run.best, best);

	fprintf(pf,"\n");
    fprintf(pf,"Mean Fitness:, %.6f",Run.meanFitness);
	fprintf(pf,"\n\n");
	
    fprintf(pf,"   Best Individual\n\n");
    fprintf(pf,"          Fitness:, %.8f\n",Run.best.fitness);
    fprintf(pf,"       Violations:,   %d\n",numTotalOutputs - Run.best.numEqual);
    fprintf(pf,"           Gates:,   %d\n",Run.best.numGates);
    fprintf(pf,"       Chromosome:,   %s\n",best);
    fprintf(pf,"   Boolean Expression:\n");
    for(i=0; i<numOutputs; i++){
        expression(Run.best.chromX, boolExpression, i);
        fprintf(pf,"%s\n",boolExpression);
	}
	fclose(pf);
}



/* Converts the chromosome (real numbers or binary) into a character string */
void chromToString(particle par, char *str)
{
	unsigned i;
	char numDouble[sizeNumDouble];

    strcpy(str,"'");
    for(i=0; i<nAllele; i++){
        dtoa(par.chromX[i], 0 , numDouble);
        strcat(str,numDouble);
        if(representation != BINARY) strcat(str," ");
	}
}

/* Converts a real number into a character string */
void dtoa(double num, unsigned prec, char *str)
{
    int dec,sign, i=0;
	double num1 = num;
	char *aux = fcvt(num,prec,&dec,&sign);

	sign = sign != 0 ? 1 : 0;
    if(sign == 1) strcpy(str,"-"); else strcpy(str,"");

	while(dec <= 0){
        str[i++] = '0';
		dec++;
	}
	for(; i<dec; i++)
        str[i] = *(aux++);
	if(prec > 0){
        str[i++] = '.';
        str[i] = EndString;
        strcat(str,aux);
		return;
	}
	
    str[i] = EndString;
}
