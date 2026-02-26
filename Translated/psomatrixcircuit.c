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
/* File: psomatrixcircuit.c                             */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Description: Main program that applies Particle Swarm */
/* Optimization to the design of logic circuits.         */
/********************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "psomatrixcircuit.h"
#include "matrixpso.h"
#include "circuits.h"
#include "statistics.h"
#include "random.h"

/* Load input parameters and execute runs */
int main(int argc, char *argv[])
{
    unsigned run;

	/*unsigned testFile, numTestFile = 3, iniFile, finFile, testRep;
	char nameTestFile[12][sizeNameFile],nameFileAux[sizeNameFile];*/

	/*strcpy(nameTestFile[0],"D:\\PSOMatrixCircuit2003A\\circuitosDTA\\circuit1Bin1S.dta");
	strcpy(nameTestFile[1],"D:\\PSOMatrixCircuit2003A\\circuitosDTA\\circuit4VEGA1S.dta");
	strcpy(nameTestFile[2],"D:\\PSOMatrixCircuit2003A\\circuitosDTA\\circuit8VEGA1S.dta");
	strcpy(nameTestFile[3],"D:\\PSOMatrixCircuit2003A\\circuitosDTA\\circuitSasaoBin1S.dta");
	strcpy(nameTestFile[4],"D:\\PSOMatrixCircuit2003A\\circuitosDTA\\circuitCatherineVEGA1S.dta");
	strcpy(nameTestFile[5],"D:\\PSOMatrixCircuit2003A\\circuitosDTA\\circuitAdderBin3S.dta");
	strcpy(nameTestFile[6],"D:\\PSOMatrixCircuit2003A\\circuitosDTA\\circuitFddiVEGA5S.dta");
	strcpy(nameTestFile[7],"D:\\PSOMatrixCircuit2003A\\circuitosDTA\\circuitMultiplierBin4S.dta");*/
	
	/*strcpy(nameTestFile[ 0],"ieee_adder_2ip");
	strcpy(nameTestFile[ 1],"miller_2ip");
	strcpy(nameTestFile[ 2],"rnd4-0");
	strcpy(nameTestFile[ 3],"rnd4-1");
	strcpy(nameTestFile[ 4],"rnd4-2");
	strcpy(nameTestFile[ 5],"rnd4-3");
	strcpy(nameTestFile[ 6],"rnd4-4");
	strcpy(nameTestFile[ 7],"rnd4-5");
	strcpy(nameTestFile[ 8],"rnd4-6");
	strcpy(nameTestFile[ 9],"rnd4-7");
	strcpy(nameTestFile[10],"rnd4-8");
	strcpy(nameTestFile[11],"rnd4-9");
	switch(atoi(argv[1])){
		//Admin2003
		case 1:
			iniFile = 0;
			finFile = 3;
			break;
		//PCSandra
		case 2:
			iniFile = 3;
			finFile = 6;
			break;
		//PCGelus
		case 3:
			iniFile = 6;
			finFile = 8;
			break;
		//PCIsra
		case 4:
			iniFile = 8;
			finFile = 10;
			break;
		//PCJoan
		case 5:
			iniFile = 10;
			finFile = 12;
	}*/

	/*strcpy(nameTestFile[ 0],"planar_11");
	strcpy(nameTestFile[ 1],"rd53_2ip");
	strcpy(nameTestFile[ 2],"roth_2ip");
	iniFile = atoi(argv[1]);
	finFile = iniFile + 1;*/
		
	//for(testFile = iniFile; testFile < finFile/*numTestFile*/; testFile++){
	/*	for(testRep=0; testRep<3; testRep++){

			strcpy(nameFileAux, nameTestFile[testFile]);
			if(testRep==0) strcat(nameFileAux,"_B.dta");
			else if(testRep==1) strcat(nameFileAux,"_EA.dta");
			else strcat(nameFileAux,"_EB.dta");*/
	//if(cargaParametros(nameFileAux)){
        if(loadParameters(argv[1])) {        
                initVariables();                     /* Initialize program variables */
                reserveMemory();                     /* Reserve required memory */
                globalHeader(nfRun);                 /* Generate header for global statistics */
                for(run = 0; run < nRun; run++) {
                    printf("\n\nRun %02d Started\n", run);
                    initStatistics(&Run);            /* Initialize run statistics */
                    pSwarm(run);                     /* Execute PSO algorithm for circuit design */
                    runStatistics(nfRun, run);       /* Generate statistics of this run */
                    printf("\n\nRun %02d Finished\n", run);
                }
                freeMemory();                        /* Free all reserved memory */
                return 1;
        }
		//}
	//}
	return 0;
}

/* Load input data from file into memory */
unsigned loadParameters(char *finput)
{
    FILE *input = fopen(finput,"r");

	if(!input){
        printf("\nError: could not load data file '%s'\n", finput);
		return 0; 
	}

    strcpy(nfInput, finput);                 /* Input filename */
    loadTT(input);                           /* Load truth table */
    tPop         = (unsigned)readNumber(input); /* Population size */
    nGen         = (unsigned)readNumber(input); /* Max generations */
    nRun         = (unsigned)readNumber(input); /* Max runs */
    tNeigh       = (unsigned)readNumber(input); /* Neighborhood size */
    phi1         = readNumber(input);           /* Acceleration constant phi1 */
    phi2         = readNumber(input);           /* Acceleration constant phi2 */
    vMax         = readNumber(input);           /* Max velocity */
    pMut         = readNumber(input)/100;       /* Mutation probability */
    representation= (unsigned)readNumber(input);/* Representation type */
	//representacion = 2;
    numGates     = (unsigned)readNumber(input); /* Number of gates */
    cardinality  = (unsigned)readNumber(input); /* Cardinality for random numbers */
    numRows      = (unsigned)readNumber(input); /* Number of rows in circuit matrix */
    numCols      = (unsigned)readNumber(input); /* Number of columns in circuit matrix */
    readString(input, nfGen);                  /* Generation statistics file */
    readString(input, nfRun);                  /* Run statistics file */
    strcat(nfRun,".csv");
	fclose(input);

    printTT();             /* Print truth table */
    printParameters();     /* Print loaded parameters */
	return 1;
}

/* Print program parameters read from input file */
void printParameters(void)
{
    printf("\n           Population size: %d", tPop);
    printf("\n         Number of generations: %d", nGen);
    printf("\n                 Number of runs: %d", nRun);
    printf("\n             Neighborhood size: %d", tNeigh);
    printf("\n              PSO phi1 parameter: %.2f", phi1);
    printf("\n              PSO phi2 parameter: %.2f", phi2);
    printf("\n                 PSO vmax: %.2f", vMax);
    printf("\n           Mutation percentage: %.2f", pMut*100);
    printf("\n                 Representation: %s", representation == BINARY ? "Binary" : representation == INTEGERA ? "Integer A" : "Integer B");	
    printf("\n      Number of available gates: %d", numGates);
    printf("\n                    Cardinality: %d", cardinality);
    printf("\n          Rows in the matrix: %d", numRows);
    printf("\n       Columns in the matrix: %d", numCols);
    printf("\n   Per-run results file: %sX.csv", nfGen);
    printf("\nGlobal results file: %s\n\n", nfRun);
}

/* Initialize program variables */
void initVariables(void)
{
    if(tNeigh > tPop) tNeigh = tPop;   /* Neighborhood size <= population size */
	
    tMat = numCols * numRows;          /* Circuit matrix size */

    initBounds();                      /* Initialize bounds for particle chromosome */
}

/* Initialize bounds for particle chromosome variables */
void initBounds(void)
{
	unsigned i;

    nVar = numRows * numCols * 3;  /* 3 vars per gate cell: input1, input2, gate type */
    lInf = memVarX;
    lSup = memVarX;

    /* Input1 + Input2 + Gate type */
    for(i=0; i<nVar; i++) {
        lInf[i] = 0;
        lSup[i] = (i+1)%3 ? numRows-1 : numGates-1;
	}

    /* Calculate the number of bits required for each variable in the binary representation */
    if(representation == BINARY) {
        bitVariable = memVarX;
        nAllele = 0;
		for(i=0; i<nVar; i++){
			bitVariable[i] = (unsigned)ceil(log((fabs(lInf[i]) + fabs(lSup[i])) + 1) / log(2));
            nAllele += bitVariable[i];
		}
	}
    else 
        nAllele = nVar;
}

/* Reserve memory required for program variables */
void reserveMemory(void)
{
	unsigned i;

    population       = memPopulation;       /* Population */
    bestSocialExp    = memPopulation;       /* Best social experiences */
    bestIndividualExp= memPopulation;       /* Best individual experiences */

    /*Reserve memory for each particle of the population*/
	for(i=0; i<tPop; i++){
        reserveParticleMemory(&population[i],0);
        reserveParticleMemory(&bestIndividualExp[i],1);
        reserveParticleMemory(&bestSocialExp[i],1);
    }
	
    /*Reserve memory to store the best and worst particle of the generation*/
	reserveParticleMemory(&(Gen.best),1);
    reserveParticleMemory(&(Gen.worst),1);

    /*Reserve memory to store the best and worst particle of the run*/
	reserveParticleMemory(&(Run.best),1);
    reserveParticleMemory(&(Run.worst),1);

    /*Memory needed for handling the matrices*/
    reserveMatrixMemory();
}

/* Reserve memory for elements of a particle */
void reserveParticleMemory(particle *par, unsigned rec)
{
	par->vi = memVi;
    par->chromX = memChromX;
}


/* Free all reserved memory after runs */
void freeMemory(void)
{
	unsigned i;

    /*Frees the memory reserved for the truth table*/
    freeMemoryTT();

    /* Frees the memory of each particle in the populations */
	for(i=0; i<tPop; i++){
        freeParticleMemory(&population[i],0);
        freeParticleMemory(&bestSocialExp[i],1);
        freeParticleMemory(&bestIndividualExp[i],1);
    }

    /* Frees the memory of the populations */
	free(population);
    free(bestSocialExp);
    free(bestIndividualExp);

    /* Frees the memory of the best and worst particle of the generation */
	freeParticleMemory(&(Gen.best),1);
    freeParticleMemory(&(Gen.worst),1);

    /* Reserves memory to store the best and worst particle of the run */
	freeParticleMemory(&(Run.best),1);
    freeParticleMemory(&(Run.worst),1);

    /* Frees the memory associated with the variables */
	free(lInf);
	free(lSup);
	free(bitVariable);

    /* Frees the memory used for matrix handling */
	reserveMatrixMemory();  // reserveMatrixMemory()
}

/* Free memory of a particle */
void freeParticleMemory(particle *par, unsigned rec)
{
	free(par->vi);
    free(par->chromX);
}

/* Main Particle Swarm program */
void pSwarm(unsigned run)
{
    unsigned gen, seed = initRandom(0);
    char fileGen[sizeFileName];

    runFileName(run, fileGen);
    initPopulation();
    runHeader(fileGen, seed);
    for(gen=0; gen<nGen; gen++){		
        initStatistics(&Gen);           /* Initializes the information for the generation statistics */
		evaluatePopulation(gen);		/* Calculates the data of each particle */
		//mutacion();
		//iniEstadistica(&Gen);			/* Initializes the information for the generation statistics */
		//aptitudPoblacion(gen);		/* Calculates the data of each particle */
		runInfo(fileGen,gen);		    /* Data obtained in the generation */
		PSOAlgorithm(gen);				/* Main algorithm */
		mutation();
	}
	runFooter(fileGen);                 /* Final results of the run */
}

/* Generate filename for current run statistics */
void runFileName(unsigned run, char *name)
{
    char num[10];

    strcpy(name,nfGen);
    dtoa(run,0,num);
    strcat(name, run == 0 ? "0" : num);
    strcat(name,".csv");
}

/* Generate initial random population */
void initPopulation(void)
{
	unsigned i,j; 

	for(i=0; i<tPop; i++)
        for(j=0; j<nAllele; j++)
            population[i].chromX[j] = representation == BINARY ? flip(0.5) : rndIR(lInf[j],(lSup[j]+1)*cardinality)%(lSup[j]+1);
}

/* Evaluate characteristics of each particle */
void evaluatePopulation(unsigned gen)
{
	unsigned i,j;
    int p, p0, p1;

    for(i=0; i<tPop; i++){
		/* Obtains the fitness of the particle */
        evaluateParticle(&population[i]);

		/* Searches for the best local particle that has been found for each particle */
        if(gen == 0 || population[i].fitness > bestIndividualExp[i].fitness)
            copyParticle(population[i],&(bestIndividualExp[i]));
	}

    for(i=0; i<tPop; i++){
		/* Searches for the best particle that has been found in the neighborhood of the particle */
        p = p0 = p1 = i;
        for(j=0; j<tNeigh; j+=2){
            p0 = p0 + 1 < (int)tPop ? p0 + 1 : 0;
            p1 = p1 - 1 < 0 ? (int)tPop - 1 : p1 - 1;
            if(bestIndividualExp[p0].fitness > bestIndividualExp[p].fitness) p = p0;
            if(bestIndividualExp[p1].fitness > bestIndividualExp[p].fitness) p = p1;
		}
        copyParticle(bestIndividualExp[p],&bestSocialExp[i]);

		/* The data for the generation statistics are calculated */
        Gen.meanFitness += bestSocialExp[i].fitness/tPop;
        Gen.squaredFitness += pow(bestSocialExp[i].fitness,2)/tPop;

        /* The best and worst particle of the generation are stored */
		if(i == 0){
            copyParticle(bestSocialExp[i],&(Gen.best));
            copyParticle(bestSocialExp[i],&(Gen.worst));
		}
		else{
            if(bestSocialExp[i].fitness > Gen.best.fitness)
                copyParticle(bestSocialExp[i],&(Gen.best));
            if(bestSocialExp[i].fitness < Gen.worst.fitness)
                copyParticle(bestSocialExp[i],&(Gen.worst));
		}
	}
}

/* Evaluate fitness of a particle */
void evaluateParticle(particle *par)
{
	/* The matrix is evaluated and the number of equal outputs and the number of WIREs are obtained */	
    evaluate(par->chromX, &(par->numEqual));

	/* Sets to zero the variable that counts the number of WIREs */
    par->numNoGates = 0;

    /* If violations exist */
    if(par->numEqual >= numTotalOutputs){
		/* The total number of gates in the solution must be counted
            and subtracted from the total number of gates in the matrix,
            that is, all those not involved are taken as WIREs */
        par->numGates = countGates();
        par->numNoGates = tMat - par->numGates;
	}
    par->fitness = (double)(par->numEqual + par->numNoGates);
}

/* Copy data from particle in to out */
void copyParticle(particle in, particle *out)
{
	unsigned i;
	
    out->numGates = in.numGates;
    out->numEqual = in.numEqual;
    out->numNoGates = in.numNoGates;
    out->fitness = in.fitness;
	
    /* Copies the chromosome and the velocities of each allele in the chromosome */
    for(i=0; i<nAllele; i++){
        out->chromX[i] = in.chromX[i];
        out->vi[i] = in.vi[i];  
	}
}

/* Apply mutation to population */
void mutation(void)
{
	unsigned i, j, flag = 0;

    for(i=0; i<tPop; i++)
		//if(poblacion[i].aptitud <= bestSocialExp[i].aptitud )
            for(j=0; j<nAllele; j++){
                if(flip(pMut)){
                    switch(representation){
                        case BINARY:
                            population[i].chromX[j] = population[i].chromX[j] ? 0 : 1;
							break;
                        case INTEGERA:
                        case INTEGERB:
                            population[i].chromX[j] = rndIR(lInf[j],(lSup[j]+1)*cardinality)%(lSup[j]+1);
							break;
					}
				}
			}
				
}

/* Compute run statistics */
void runInfo(char *file, unsigned gen)
{	
    /* Calculates the mean and the variance data */
    Run.meanFitness += Gen.meanFitness/nGen;
    Run.squaredFitness += Gen.squaredFitness/nGen;

    /* Stores the best and worst individual of the run */
    if(gen == 0){
        copyParticle(Gen.best,&(Run.best));
        copyParticle(Gen.worst,&(Run.worst));
        /* If a better result was found, it is written to the file and displayed */
        generationStats(file,gen);                
        printf("\nGeneration: %04d -- Fitness: %f",gen,Gen.best.fitness);
        Run.generation = gen;
	}
	else{
        if(Gen.best.fitness > Run.best.fitness){
            copyParticle(Gen.best, &(Run.best));
            /* If a better result was found, it is written to the file and displayed */
            Run.generation = gen;
		}
        generationStats(file,gen);                
        printf("\nGeneration: %04d -- Fitness: %f",gen,Gen.best.fitness);

        if(Gen.worst.fitness < Run.worst.fitness)
            copyParticle(Gen.worst, &(Run.worst));
	}
}

/* Particle Swarm Optimization algorithm */
void PSOAlgorithm(unsigned gen)
{
	unsigned i,d;
	double phi1p, phi2p, vMaxNorm;

    for(i=0; i<tPop; i++){
		/* Calculates the individual predisposition */
        for(d = 0; d < nAllele; d++){
			phi1p = rndF()*phi1;
			phi2p = rndF()*phi2;

			if(gen){
                population[i].vi[d] += phi1p * ((double)bestIndividualExp[i].chromX[d] - (double)population[i].chromX[d]); 
                population[i].vi[d] += phi2p * ((double)bestSocialExp[i].chromX[d] - (double)population[i].chromX[d]);
			}
			else{
				//if(representacion == BINARIA){
                    if(!bestSocialExp[i].chromX[d] && !bestIndividualExp[i].chromX[d])
                        population[i].vi[d] = -vMax;
					else
                        if(bestSocialExp[i].chromX[d] && bestIndividualExp[i].chromX[d])							
                            population[i].vi[d] = vMax;
						else
                            population[i].vi[d] = 0;
				//}
				//else
				//	poblacion[i].vi[d] = rndFR(-vMax,vMax);

				population[i].vi[d] += phi1p * ((double)bestIndividualExp[i].chromX[d]); /*- (double)Gen[gen].mejorInd.cromosoma[d]*/
				population[i].vi[d] += phi2p * ((double)bestSocialExp[i].chromX[d]); /*- (double)Gen[gen].mejorInd.cromosoma[d]*/

			}

            population[i].vi[d] = 
                population[i].vi[d] > vMax ? vMax : 
                population[i].vi[d] < -vMax ? -vMax : population[i].vi[d];

            vMaxNorm = sigmoid(population[i].vi[d]);
            switch(representation){
                case BINARY:
                    population[i].chromX[d] = flip(vMaxNorm);
					break;
                case INTEGERA:
                    population[i].chromX[d] = flip(vMaxNorm) ? bestSocialExp[i].chromX[d] : population[i].chromX[d];
					break;
                case INTEGERB:
                    population[i].chromX[d] = flip(vMaxNorm) ? bestSocialExp[i].chromX[d] : flip(1 - vMaxNorm) ? bestIndividualExp[i].chromX[d] : population[i].chromX[d];
					break;

			}
		}
	}

	/*for(i=0; i<tPob; i++)
		for(d = 0; d < nAlelo; d++){
			vMaxNorm = sigmoide(poblacion[i].vi[d]);
			switch(representacion){
				case BINARIA:
					poblacion[i].cromX[d] = flip(vMaxNorm);
					break;
				case ENTERAA:
					poblacion[i].cromX[d] = flip(vMaxNorm) ? bestSocialExp[i].cromX[d] : poblacion[i].cromX[d];
					break;
			}
		}*/
}