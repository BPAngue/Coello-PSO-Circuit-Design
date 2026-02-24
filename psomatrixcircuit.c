/********************************************************/
/*					CINVESTAV - IPN						*/
/*			Departamento de Ingenería Eléctrica			*/
/*					Sección Computación					*/
/*														*/	
/*					Computacion Evolutiva				*/
/*														*/
/*					Erika Hernandez Luna				*/
/*			eluna@computacion.cs.cinvestav.mx			*/
/*					 2 / agosto / 2003					*/
/*														*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Archivo: psomatrixcircuit.c							*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Descripción:											*/
/********************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "psomatrixcircuit.h"
#include "matrixpso.h"
#include "circuitos.h"
#include "estadisticas.h"
#include "random.h"

/*Se carga los parametros de entrada y se realizan las corridas*/
int main(int argc, char *argv[])
{
	unsigned corr;

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
		if(cargaParametros(argv[1])){		
				iniVariables();						/*Se inicializan las variables necesarias para el programa*/
				reservaMemoria();					/*Se reserva la memoria necesaria para el programa*/
				estEncGlobal(nfCorr);				/*Se genera el encabezado de la estadistica global*/
				for(corr=0; corr<nCorr; corr++){
					printf("\n\nCorrida %02d Iniciada\n",corr);
					iniEstadistica(&Corr);			/*Inicializa las variables de la estadística global*/
					pSwarm(corr);					/*Algoritmo de PSO para diseño de circuitos*/
					estCorr(nfCorr,corr);			/*Genera las estadísticas de la corrida*/
					printf("\n\nCorrida %02d Terminada\n",corr);
				}
				liberaMemoria();					/*Se libera toda la memoria reservada*/
				return 1;
		}
		//}
	//}
	return 0;
}

/*Almacenan en la memoria los datos obtenidos del archivo de entrada*/
unsigned cargaParametros(char *finput)
{
	FILE *input = fopen(finput,"r");

	if(!input){
		printf("\nHa ocurrido un error al tratar de cargar el archivo de datos '%s'\n", finput);
		return 0; 
	}

	strcpy(nfEntrada,finput);							/*Nombre del archivo de entrada*/
	cargaTT(input);										/*Coloca los datos de la tabla de verdad en variables*/
	tPob =				(unsigned)leeNumero(input);		/*Tamaño de la poblacion*/ 
	nGen =				(unsigned)leeNumero(input);		/*Numero maximo de generaciones*/
	nCorr =				(unsigned)leeNumero(input);		/*Numero maximo de corridas*/
	tVec =				(unsigned)leeNumero(input);		/*Tamaño del vecinario*/
	phi1 =				leeNumero(input);				/*Constante de aceleracion 1, phi1*/
	phi2 =				leeNumero(input);				/*Constante de aceleracion 2, phi2*/
	vMax =				leeNumero(input);				/*Constante del sistema Vmax*/			
	pMut =				leeNumero(input)/100;			/*Porcentaje de mutación*/
	representacion =	(unsigned)leeNumero(input);		/*Representación Binaria=2, Entera A = 1, Entera B = 0*/
	//representacion = 2;
	numGates =			(unsigned)leeNumero(input);		/*Numero de compuertas disponibles en el programa*/
	cardinalidad =		(unsigned)leeNumero(input);		/*Cardinalidad para los números aleatorios de cada variables*/
	numReng =			(unsigned)leeNumero(input);		/*Numero de renglones en la matriz del circuito*/
	numCols =			(unsigned)leeNumero(input);		/*Numero de columnas en la matriz del circuito*/
	leeCadena(input,nfGen);								/*Estadisticas detalladas por Generación*/
	leeCadena(input,nfCorr);							/*Estadisticas detalladas por Corrida*/
	strcat(nfCorr,".csv");
	fclose(input);

	/*Imprime los datos referentes a la tabla de verdad*/
	imprimeTT();
	/*Imprime los datos recibidos*/
	imprimeParametros();			
	return 1;
}

/*Imprime en pantalla los parametros de entrada al programa recibidos en el archivo*/
void imprimeParametros(void)
{
	printf("\n           Tamaño de la poblacion: %d",      tPob);
	printf("\n           Numero de generaciones: %d",      nGen);
	printf("\n               Numero de corridas: %d",      nCorr);
	printf("\n            Tamaño del vecindario: %d",      tVec);
	printf("\n           Parametro phi1 del PSO: %.2f",    phi1);
	printf("\n           Parametro phi2 del PSO: %.2f",    phi2);
	printf("\n           Parametro vmax del PSO: %.2f",    vMax);
	printf("\n           Porcentaje de mutacion: %.2f",    pMut*100);
	printf("\n                   Representacion: %s",      representacion == BINARIA ? "Binaria" : representacion == ENTERAA ? "Entera A" : "Entero B");
	printf("\n Numero de compuertas disponibles: %d",	   numGates);
	printf("\n                     Cardinalidad: %d",	   cardinalidad);
	printf("\n Numero de Renglones de la matriz: %d",	   numReng);
	printf("\n  Numero de Columnas de la matriz: %d",	   numCols);
	printf("\nArchivo de resultados por corrida: %sX.csv", nfGen);
	printf("\n   Archivo de resultados globales: %s\n\n",  nfCorr);
}

/*Inicializa las variables necesarias para el programa*/
void iniVariables(void)
{
	/*El tam del vecindario debe ser menor o igual al tam de la poblacion*/
	if(tVec > tPob) tVec = tPob;			
	
	/*El tamaño de la matriz que esta representando al circuito*/
	tMat = numCols*numReng;

	/*Inicializa el rangos de las variables del cromosoma de la particula*/
	iniLimites();
}

/*Inicializa el rangos de las variables del cromosoma de la particula*/
void iniLimites(void)
{
	unsigned i;

	/*Hay tantas variables como elementos tenga la matriz de compuertas y cada una
	  de esas variables, tiene 3 variables a su vez, para codificar la entrada1 a la
	  compuerta, la entrada2 a la compuerta y por ultimo el tipo de compuerta.*/
	nVar = numReng * numCols * 3; 
	lInf = memVarX;
	lSup = memVarX;

	/*Entrada1 + Entrada2 + Tipo*/
	for (i=0; i<nVar; i++){
		lInf[i]=0;
		lSup[i] = (i+1)%3 ? numReng-1 : numGates-1;
	}

	/*Calcula el numero de bits que ocupa cada variable en la representacion binaria*/
	if(representacion == BINARIA){
		bitVariable = memVarX;
		nAlelo = 0;
		for(i=0; i<nVar; i++){
			bitVariable[i] = (unsigned)ceil(log((fabs(lInf[i]) + fabs(lSup[i])) + 1) / log(2));
			nAlelo += bitVariable[i];
		}
	}
	else 
		nAlelo = nVar;
}

/*Reserva la memoria necesaria para las variables del programa*/
void reservaMemoria(void)
{
	unsigned i;

	poblacion = memPoblacion;			/*Poblacion regular*/
	bestSocialExp = memPoblacion;		/*Mejores experiencias sociales*/
	bestIndividualExp = memPoblacion;	/*Mejores experiencias individuales*/

	/*Reserva memoria para cada partícula de la poblacion*/
	for(i=0; i<tPob; i++){
		reservaMemoriaParticula(&poblacion[i],0);
		reservaMemoriaParticula(&bestIndividualExp[i],1);
		reservaMemoriaParticula(&bestSocialExp[i],1);
	}
	
	/*Reserva memoria para almacenar a la mejor y a la peor particula de la generacion*/
	reservaMemoriaParticula(&(Gen.mejor),1);
	reservaMemoriaParticula(&(Gen.peor),1);

	/*Reserva memoria para almacenar a la mejor y a la peor particula de la corrida*/
	reservaMemoriaParticula(&(Corr.mejor),1);
	reservaMemoriaParticula(&(Corr.peor),1);

	/*Memoria necesaria para el manejo de las matrices*/
	reservaMemoriaMatriz();

}

/*Reserva memoria para los elementos de la particula*/
void reservaMemoriaParticula(particula *par, unsigned rec)
{
	par->vi = memVi;
	par->cromX = memCromX;
}


/*Libera toda la memoria reservada a lo largo del programa una vez terminadas las corridas*/
void liberaMemoria(void)
{
	unsigned i;

	/*Libera la memoria reservada para la tabla de verdad*/
	liberaMemoriaTT();

	/*Libera la memoria de cada particula en las poblaciones*/
	for(i=0; i<tPob; i++){
		liberaMemoriaParticula(&poblacion[i],0);
		liberaMemoriaParticula(&bestSocialExp[i],1);
		liberaMemoriaParticula(&bestIndividualExp[i],1);
	}

	/*Libera la memoria de las poblaciones*/
	free(poblacion);
	free(bestSocialExp);
	free(bestIndividualExp);

	/*Libera la memoria de la mejor y de la peor particula de la generacion*/
	liberaMemoriaParticula(&(Gen.mejor),1);
	liberaMemoriaParticula(&(Gen.peor),1);

	/*Reserva memoria para almacenar a la mejor y a la peor particula de la corrida*/
	liberaMemoriaParticula(&(Corr.mejor),1);
	liberaMemoriaParticula(&(Corr.peor),1);

	/*Libera la memoria asociada a las variables*/
	free(lInf);
	free(lSup);
	free(bitVariable);

	/*Libera la memoria usada para el manejo de las matrices*/
	reservaMemoriaMatriz();
}

/*Libera la memoria reservada para cada partícula*/
void liberaMemoriaParticula(particula *par, unsigned rec)
{
	free(par->vi);
	free(par->cromX);
}

/*Programa principal del Particle Swarm*/
void pSwarm(unsigned corr)
{
	unsigned gen, semilla = initrandom(0);
	char nfileGen[sizeNameFile];

	nameFileCorr(corr,nfileGen);		/*Obtiene el nombre del archivo de estadisticas de la corrida*/
	poblacionInicial();					/*Poblacion inicial de la corrida*/
	estEncCorr(nfileGen,semilla);		/*Encabezado de la estadistica de la corrida*/
	for(gen=0; gen<nGen; gen++){
		iniEstadistica(&Gen);			/*Se inicializa la información para la estadística de la generacion*/
		aptitudPoblacion(gen);			/*Calcula los datos de cada particula*/
		//mutacion();
		//iniEstadistica(&Gen);			/*Se inicializa la información para la estadística de la generacion*/
		//aptitudPoblacion(gen);		/*Calcula los datos de cada particula*/
		infoCorrida(nfileGen,gen);		/*Datos que se obtuvieron en la generacion*/
		algoritmoPSO(gen);				/*Algoritmo principal*/
		mutacion();
	}
	estPieCorr(nfileGen);				/*Resultados finales de la corrida*/
}

/*Genera el nombre del archivo para la estadistica de la corrida actual*/
void nameFileCorr(unsigned corr, char *nombre)
{
	char num[10];

	strcpy(nombre,nfGen);
	dtoa(corr,0,num);
	strcat(nombre,corr == 0 ? "0" : num);
	strcat(nombre,".csv");
}

/*Genera la poblacion aleatoria inicial de individuos*/
void poblacionInicial(void)
{
	unsigned i,j; 

	for(i=0; i<tPob; i++)
		for(j=0; j<nAlelo; j++)
			poblacion[i].cromX[j] = representacion == BINARIA ? flip(0.5) : rndIR(lInf[j],(lSup[j]+1)*cardinalidad)%(lSup[j]+1);
}

/*Calcula las caracteristicas de cada particula*/
void aptitudPoblacion(unsigned gen)
{
	unsigned i,j;
	int par, par0, par1;

	for(i=0; i<tPob; i++){
		/*Obtiene la aptitud de la particula*/
		aptitudParticula(&poblacion[i]);

		/*Busca a la mejor particula local que se ha encontrado para cada particula*/
		if(gen == 0 || poblacion[i].aptitud > bestIndividualExp[i].aptitud)
			copiaParticula(poblacion[i],&(bestIndividualExp[i]));
	}

	for(i=0; i<tPob; i++){
		/*Busca a la mejor partícula que se ha encontrado en el vecindario de la particula*/
		par = par0 = par1 = i;
		for(j=0; j<tVec; j+=2){
			par0 = par0 + 1 < (int)tPob ? par0 + 1 : 0;
			par1 = par1 - 1 < 0 ? (int)tPob - 1 : par1 - 1;
			if(bestIndividualExp[par0].aptitud > bestIndividualExp[par].aptitud ) par = par0;
			if(bestIndividualExp[par1].aptitud > bestIndividualExp[par].aptitud ) par = par1;
		}
		copiaParticula(bestIndividualExp[par],&bestSocialExp[i]);

		/*Se calculan los datos para las estadisticas de la generacion*/
		Gen.mediaAptitud += bestSocialExp[i].aptitud/tPob;
		Gen.cuadradosAptitud += pow(bestSocialExp[i].aptitud,2)/tPob;

	   /*Se almacena a la mejor y a la peor particula de la generacion*/
		if(i == 0){
			copiaParticula(bestSocialExp[i],&(Gen.mejor));
			copiaParticula(bestSocialExp[i],&(Gen.peor));
		}
		else{
			if(bestSocialExp[i].aptitud > Gen.mejor.aptitud)
				copiaParticula(bestSocialExp[i],&(Gen.mejor));
			if(bestSocialExp[i].aptitud < Gen.peor.aptitud)
				copiaParticula(bestSocialExp[i],&(Gen.peor));
		}
	}
}

/*Obtiene la aptitud de la partícula*/
void aptitudParticula(particula *par)
{
	/*Se evalua la matriz y se obtienen el número de salidas iguales y el número de WIRE's*/	
	evalua(par->cromX, &(par->numIgual));

	/*Pone a cero la variable que cuenta el numero de WIRE's*/
	par->numNoGates = 0;

	/*Si existen violaciones*/
	if(par->numIgual >= numTSalidas){
		/*Hay que contar el numero total de compuertas que estan en la solucion
		y restarlo del numero total de compuertas en la matriz, es decir, tomar 
		a todos los no involucrados como WIRE's*/
		par->numGates = cuentaCompuertas();
		par->numNoGates = tMat - par->numGates;
	}
	par->aptitud = (double)(par->numIgual + par->numNoGates);
}

/*Copia los datos de la particula in en out*/
void copiaParticula(particula in,particula *out)
{
	unsigned i;
	
	out->numGates = in.numGates;
	out->numIgual = in.numIgual;
	out->numNoGates = in.numNoGates;
	out->aptitud = in.aptitud;
	
	/*Copia el cromosoma y las velocidades de cada alelo en el cromosoma*/
	for(i=0; i<nAlelo; i++){
		out->cromX[i] = in.cromX[i];
		out->vi[i] = in.vi[i];  
	}
}


void mutacion(void)
{
	unsigned i, j, flag = 0;

	for(i=0; i<tPob; i++)
		//if(poblacion[i].aptitud <= bestSocialExp[i].aptitud )
			for(j=0; j<nAlelo; j++){
				if(flip(pMut)){
					switch(representacion){
						case BINARIA:
							poblacion[i].cromX[j] = poblacion[i].cromX[j] ? 0 : 1;
							break;
						case ENTERAA:
						case ENTERAB:
							poblacion[i].cromX[j] = rndIR(lInf[j],(lSup[j]+1)*cardinalidad)%(lSup[j]+1);
							break;
					}
				}
			}
				
}

/*Se calculan los datos para las estadisticas de la corrida*/
void infoCorrida(char *nfile, unsigned gen)
{	
	/*Calcula la media y los datos de la varianza*/
	Corr.mediaAptitud += Gen.mediaAptitud/nGen;
	Corr.cuadradosAptitud += Gen.cuadradosAptitud/nGen;

	/*Se almacena el mejor y el peor individuo de la corrida*/
	if(gen == 0){
		copiaParticula(Gen.mejor,&(Corr.mejor));
		copiaParticula(Gen.peor,&(Corr.peor));
		/*Si se encontró un resultado mejor, se escribe al archivo y se despliega*/
		estGen(nfile,gen);				
		printf("\nGeneracion: %04d -- Aptitud: %f",gen,Gen.mejor.aptitud);
		Corr.generacion = gen;
	}
	else{
		if(Gen.mejor.aptitud > Corr.mejor.aptitud){
			copiaParticula(Gen.mejor, &(Corr.mejor));
			/*Si se encontró un resultado mejor, se escribe al archivo y se despliega*/				
			Corr.generacion = gen;
		}
		estGen(nfile,gen);				
		printf("\nGeneracion: %04d -- Aptitud: %f",gen,Gen.mejor.aptitud);

		if(Gen.peor.aptitud < Corr.peor.aptitud)
			copiaParticula(Gen.peor, &(Corr.peor));
	}
}

/*Algoritmo de Particle Swarm Optimization*/
void algoritmoPSO(unsigned gen)
{
	unsigned i,d;
	double phi1p, phi2p, vMaxNorm;

	for(i=0; i<tPob; i++){
		/*Calcula la predisposicion individual*/
		for(d = 0; d < nAlelo; d++){
			phi1p = rndF()*phi1;
			phi2p = rndF()*phi2;

			if(gen){
				poblacion[i].vi[d] += phi1p * ((double)bestIndividualExp[i].cromX[d] - (double)poblacion[i].cromX[d]); 
				poblacion[i].vi[d] += phi2p * ((double)bestSocialExp[i].cromX[d] - (double)poblacion[i].cromX[d]);
			}
			else{
				//if(representacion == BINARIA){
					if(!bestSocialExp[i].cromX[d] && !bestIndividualExp[i].cromX[d])
						poblacion[i].vi[d] = -vMax;
					else
						if(bestSocialExp[i].cromX[d] && bestIndividualExp[i].cromX[d])
							poblacion[i].vi[d] = vMax;
						else
							poblacion[i].vi[d] = 0;
				//}
				//else
				//	poblacion[i].vi[d] = rndFR(-vMax,vMax);

				poblacion[i].vi[d] += phi1p * ((double)bestIndividualExp[i].cromX[d] /*- (double)Gen[gen].mejorInd.cromosoma[d]*/); 
				poblacion[i].vi[d] += phi2p * ((double)bestSocialExp[i].cromX[d] /*- (double)Gen[gen].mejorInd.cromosoma[d]*/);

			}

			poblacion[i].vi[d] = 
				poblacion[i].vi[d] > vMax ? vMax : 
				poblacion[i].vi[d] < -vMax ? -vMax : poblacion[i].vi[d];

			vMaxNorm = sigmoide(poblacion[i].vi[d]);
			switch(representacion){
				case BINARIA:
					poblacion[i].cromX[d] = flip(vMaxNorm);
					break;
				case ENTERAA:
					poblacion[i].cromX[d] = flip(vMaxNorm) ? bestSocialExp[i].cromX[d] : poblacion[i].cromX[d];
					break;
				case ENTERAB:
					poblacion[i].cromX[d] = flip(vMaxNorm) ? bestSocialExp[i].cromX[d] : flip(1 - vMaxNorm) ? bestIndividualExp[i].cromX[d] : poblacion[i].cromX[d];
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