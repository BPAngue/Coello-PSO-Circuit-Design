/********************************************************/
/*					CINVESTAV - IPN						*/
/*			Departamento de Ingenería Eléctrica			*/
/*					Sección Computación					*/
/*														*/	
/*				   Computacion Evolutiva				*/
/*														*/
/*					Erika Hernandez Luna				*/
/*			 eluna@computacion.cs.cinvestav.mx			*/
/*					  2 / agosto / 2003					*/
/*														*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Archivo: estadisticas.c								*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Descripción: Biblioteca para generar las estadisticas*/
/* por generacion y por	corrida del particle swarm		*/
/* optimization. Las estadisticas incluyen a la mejor   */
/* y a la peor particula, la media de la aptitud, la 	*/
/* mejor y la peor aptitud, asi como la desviacion		*/
/* estandar, la	varianza y la cadena cromosomica de la	*/
/* mejor particula encontrada por generacion.			*/
/********************************************************/

#include "estadisticas.h"
#include "psomatrixcircuit.h"
#include "circuitos.h"
#include "matrixpso.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/*Inicializa la informacion de la estadística de la corrida o de la generacion*/
void iniEstadistica(estadistica *est)
{
	est->cuadradosAptitud = 0;
	est->mediaAptitud  = 0;
}

/*Genera el encabezado de la estadistica global de todas las corridas*/
void estEncGlobal(char *nfile)
{
	FILE *pf = fopen(nfile,"w");

	fprintf(pf,"Centro de Investigación y de Estudios Avanzados del IPN\n");
	fprintf(pf,"Diseño de Circuitos Lógicos Combinatorios usando Optimización mediante Cúmulos de Partículas\n");
	fprintf(pf,"Departamento de Ingeniería Eléctrica - Sección Computación \n");
	fprintf(pf,"Erika Hernández Luna  eluna@computacion.cs.cinvestav.mx\n\n");
	fprintf(pf,"Estadisticas de los resultados obtenidos en todas las corridas del programa\n\n");

	fprintf(pf,"     Archivo Procesado:,   %s\n",nfEntrada);
	fprintf(pf,"   Tamaño de Poblacion:,   %u\n",tPob);
	fprintf(pf,"Numero de Generaciones:,   %u\n",nGen);
	fprintf(pf,"                 Phi 1:, %.3f\n",phi1);
	fprintf(pf,"                 Phi 2:, %.3f\n",phi2);
	fprintf(pf,"                  Vmax:, %.3f\n",vMax);
	fprintf(pf,"Porcentaje de mutación:, %.3f\n",pMut*100);
	fprintf(pf,"  Número de Compuertas:,   %u\n",numGates);
	fprintf(pf,"          Cardinalidad:,   %u\n",cardinalidad);
	fprintf(pf,"Renglones en la matriz:,   %u\n",numReng);
	fprintf(pf," Columnas en la matriz:,   %u\n",numCols);	
	fprintf(pf,"            Vecindario:,   %u\n",tVec);
	switch(representacion){
		case BINARIA:
			fprintf(pf,"        Representacion:, Binaria\n");
			break;
		case ENTERAA:
			fprintf(pf,"        Representacion:, Entera A\n");
			break;
		case ENTERAB:
			fprintf(pf,"        Representacion:, Entera B\n");
	}

	fprintf(pf,"\n");
	fprintf(pf,"Corrida,Media Aptitud,Generacion,");
	fprintf(pf,"Mejor Cromosoma,Aptitud,");
	fprintf(pf,"Peor Cromosoma,Aptitud,");
	fprintf(pf,"Desviacion Estandar,Varianza, Expresion\n");
	fclose(pf);
}

/*Genera la estadistica de cada corrida en el archivo global de corridas*/
void estCorr(char *nfile, unsigned corr)
{
	unsigned i;
	double varianza;
	char mejor[sizeCrom], peor[sizeCrom], expBooleana[sizeExpBooleana];
	FILE *pf = fopen(nfile,"a+");

	varianza = Corr.cuadradosAptitud - Corr.mediaAptitud*Corr.mediaAptitud;
	cadenaCrom(Corr.mejor, mejor);
	cadenaCrom(Corr.peor, peor);

	fprintf(pf,"%u,%.8f,%u,%s,%.8f,%s,%.8f,%.8f,%.12f",
		corr,
		Corr.mediaAptitud,
		Corr.generacion,
		mejor,
		Corr.mejor.aptitud,
		peor,
		Corr.peor.aptitud,
		sqrt(varianza),
		varianza);

	for(i=0; i<numSalidas; i++){
		expresion(Corr.mejor.cromX, expBooleana, i);
		fprintf(pf,",%s",expBooleana);
	}
	fprintf(pf,"\n");
	fclose(pf);

}

/*Genera el encabezado de la estadistica de la corrida*/
void estEncCorr(char *nfile,unsigned semilla)
{
	FILE *pf = fopen(nfile,"w");

	fprintf(pf,"Centro de Investigación y de Estudios Avanzados del IPN\n");
	fprintf(pf,"Diseño de Circuitos Lógicos Combinatorios usando Optimización mediante Cúmulos de Partículas\n");
	fprintf(pf,"Departamento de Ingeniería Eléctrica - Sección Computación \n");
	fprintf(pf,"Erika Hernández Luna  eluna@computacion.cs.cinvestav.mx\n\n");
	fprintf(pf,"Resultados obtenidos en una corrida del programa\n\n");
	
	fprintf(pf,"     Archivo Procesado:,   %s\n",nfEntrada);
	fprintf(pf,"   Tamaño de Poblacion:,   %u\n",tPob);
	fprintf(pf,"Numero de Generaciones:,   %u\n",nGen);
	fprintf(pf,"                 Phi 1:, %.3f\n",phi1);
	fprintf(pf,"                 Phi 2:, %.3f\n",phi2);
	fprintf(pf,"                  Vmax:, %.3f\n",vMax);
	fprintf(pf,"Porcentaje de mutación:, %.3f\n",pMut*100);
	fprintf(pf,"  Número de Compuertas:,   %u\n",numGates);
	fprintf(pf,"          Cardinalidad:,   %u\n",cardinalidad);
	fprintf(pf,"Renglones en la matriz:,   %u\n",numReng);
	fprintf(pf," Columnas en la matriz:,   %u\n",numCols);	
	fprintf(pf,"            Vecindario:,   %u\n",tVec);
	fprintf(pf,"               Semilla:,   %u\n",semilla);
	switch(representacion){
		case BINARIA:
			fprintf(pf,"        Representacion:, Binaria\n");
			break;
		case ENTERAA:
			fprintf(pf,"        Representacion:, Entera A\n");
			break;
		case ENTERAB:
			fprintf(pf,"        Representacion:, Entera B\n");
	}

	fprintf(pf,"\n");
	fprintf(pf,"Generacion,Media Aptitud,Mejor Aptitud,Peor Aptitud,");
	fprintf(pf,"Mejor Cromosoma");
	fprintf(pf,"\n");
	fclose(pf);
}

/*Genera la estadistica de una generaciones de la corrida*/
void estGen(char *nfile, unsigned gen)
{
	char mejor[sizeCrom];
	FILE *pf = fopen(nfile,"a+");

	cadenaCrom(Gen.mejor,mejor);
	fprintf(pf,"%u,%.8f,%.8f,%.8f,%s,",
		gen,
		Gen.mediaAptitud,
		Gen.mejor.aptitud,
		Gen.peor.aptitud,
		mejor);
	fprintf(pf,"\n");
	fclose(pf);
}

/*Genera los resultados finales de la estadistica de la corrida*/
void estPieCorr(char *nfile)
{
	unsigned i;
	char mejor[sizeCrom];
	cadBooleana expBooleana;
	FILE *pf = fopen(nfile,"a+");

	cadenaCrom(Corr.mejor, mejor);

	fprintf(pf,"\n");
	fprintf(pf,"Media Aptitud:, %.6f",Corr.mediaAptitud);
	fprintf(pf,"\n\n");
	
	fprintf(pf,"   Mejor Individuo\n\n");
	fprintf(pf,"              Aptitud:, %.8f\n",Corr.mejor.aptitud);
	fprintf(pf,"          Violaciones:,   %d\n",numTSalidas - Corr.mejor.numIgual);
	fprintf(pf,"           Compuertas:,   %d\n",Corr.mejor.numGates);
	fprintf(pf,"             Cromosma:,   %s\n",mejor);
	fprintf(pf,"   Expresion Booleana:\n");
	for(i=0; i<numSalidas; i++){
		expresion(Corr.mejor.cromX, expBooleana, i);
		fprintf(pf,"%s\n",expBooleana);
	}
	fclose(pf);
}



/*Convierte el cromosoma de números reales o binario a una cadena de caracteres*/
void cadenaCrom(particula par, char *cadena)
{
	unsigned i;
	char numDouble[sizeNumDouble];

    strcpy(cadena,"'");
	for(i=0; i<nAlelo; i++){
		dtoa(par.cromX[i], 0 , numDouble);
		strcat(cadena,numDouble);
		if(representacion != BINARIA) strcat(cadena," ");
	}
}

/*Convierte un número real en una cadena de caracteres*/
void dtoa(double num, unsigned prec, char *cadena)
{
	int dec,sign, i=0;
	double num1 = num;
	char *aux = fcvt(num,prec,&dec,&sign);

	sign = sign != 0 ? 1 : 0;
	if(sign == 1) strcpy(cadena,"-"); else strcpy(cadena,"");

	while(dec <= 0){
		cadena[i++] = '0';
		dec++;
	}
	for(; i<dec; i++)
		cadena[i] = *(aux++);
	if(prec > 0){
		cadena[i++] = '.';
		cadena[i] = EndString;
		strcat(cadena,aux);
		return;
	}
	
	cadena[i] = EndString;
}
