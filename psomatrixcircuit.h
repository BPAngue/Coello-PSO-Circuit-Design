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
/* Archivo: psomatrixcircuit.h							*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Descripción: Archivo de encabezado de la biblioteca	*/
/* psomatrixcircuit.c.									*/	
/********************************************************/

#include <malloc.h>
#include <stdlib.h>

#ifndef PSWARM
#define PSWARM

/*Macros para la reserva de memoria*/
#define memPoblacion  (particula *)malloc(sizeof(particula)*tPob)			/*Reserva de memoria para las poblaciones*/
#define memVarX		  (unsigned *)malloc(sizeof(unsigned)*nVar)				/*Reserva de memoria para las variables del individuo*/
#define memCromX	  (unsigned *)malloc(sizeof(unsigned)*nAlelo)			/*Reserva de memoria para el cromosoma del individuo*/
#define memVi         (double *)malloc(sizeof(double)*nAlelo)				/*Reserva de memoria para las velocidades de cada uno de los alelos del cromosoma*/

#define sizeNameFile 200	/*Tamaño maximo de caracteres para el nombre de los archivos*/
#define BINARIA	2			/*Representación Binaria*/
#define ENTERAA	1			/*Representacion entera A*/
#define ENTERAB 0			/*Representación entera B*/

/*Estructura de una particula de la poblacion*/
typedef struct{					
	/***Informacion usada por el PSO***/
	double *vi;				/*Velocidades para cada uno de los alelos del cromosoma*/
	unsigned *cromX;		/*Es el cromosoma que va a estar evolucionando en el pso, el que tiene la probabilidad de mutacion de cada nodo*/

	/***Información del individuo***/
	unsigned numGates;		/*Número de compuertas utilizadas en la solución*/
	unsigned numNoGates;	/*Número de wire's en la solución*/
	unsigned numIgual;		/*Número de salidas iguales entre la tabla de verdad dada por el circuito y la tabla de verdad de entrada*/
	double aptitud;			/*Aptitud de la particula*/
}particula;

/*Variables con los nombres de los archivos*/
char nfEntrada[sizeNameFile];	/*Nombre del archivo de entrada de datos*/
char nfCorr[sizeNameFile];		/*Nombre del archivo para almacenar las corridas*/
char nfGen[sizeNameFile];		/*Nombre del archivo para almacenar las generacion*/

/*Variables introducidas por el usuario*/
unsigned tPob;					/*Tamaño de la poblacion*/
unsigned nGen;					/*Numero maximo de generaciones*/
unsigned nCorr;					/*Numero maximo de corridas*/
unsigned tVec;					/*Tamaño del vecindario*/
double phi1;					/*Constante de aceleración 1*/
double phi2;					/*Constante de aceleración 2*/
double vMax;					/*Constante de velocidad maxima y minima en el PSO*/
double pMut;					/*Porcentaje de mutación*/
unsigned representacion;		/*2 representación Binaria, 1 Entera A, 0 Entera B*/
unsigned cardinalidad;

/*Variables usadas por el programa*/
unsigned nVar;					/*Numero de variables del problema*/
unsigned nAlelo;				/*Numero de alelos en el cromosoma*/	
unsigned *lInf;					/*Limite inferior del rango de cada variable*/
unsigned *lSup;					/*Limite superior del rango de cada variable*/
unsigned *bitVariable;			/*Numero de bits que ocupa cada variable en el cromosoma*/

/*Variables con Poblaciones usadas en el programa*/
particula *poblacion;			/*Estructura que almacena una poblacion*/
particula *bestSocialExp;		/*Estrcutura que almacena la poblacion de las mejores experiencias sociales*/
particula *bestIndividualExp;	/*Estrcutura que almacena la poblacion de las mejores experiencias individuales*/

/*Funciones utilizadas por el programa*/
unsigned cargaParametros(char*);					/*Almacenan en la memoria los datos obtenidos del archivo de entrada*/
void imprimeParametros(void);						/*Imprime en pantalla los parametros de entrada al programa recibidos en el archivo*/
void iniVariables(void);							/*Inicializa las variables necesarias para el programa*/
void iniLimites(void);								/*Inicializa el rangos de las variables del cromosoma de la particula*/
void reservaMemoria(void);							/*Reserva la memoria necesaria para las variables del programa*/
void reservaMemoriaParticula(particula*,unsigned);	/*Reserva memoria para los elementos de la particula*/
void liberaMemoria(void);							/*Libera toda la memoria reservada a lo largo del programa una vez terminadas las corridas*/
void liberaMemoriaParticula(particula*,unsigned);	/*Libera la memoria reservada para cada partícula*/
void pSwarm(unsigned);								/*Programa principal del Particle Swarm*/
void nameFileCorr(unsigned,char*);					/*Genera el nombre del archivo para la estadistica de la corrida actual*/
void poblacionInicial(void);						/*Genera la poblacion aleatoria inicial de individuos*/
void aptitudPoblacion(unsigned);					/*Calcula las caracteristicas de cada particula*/
void aptitudParticula(particula*);					/*Obtiene la aptitud de la partícula*/
void copiaParticula(particula,particula*);			/*Copia los datos de la particula in en out*/
void mutacion(void);
void infoCorrida(char*,unsigned);					/*Se calculan los datos para las estadisticas de la corrida*/
void algoritmoPSO(unsigned);						/*Algoritmo de Particle Swarm Optimization*/
#endif
