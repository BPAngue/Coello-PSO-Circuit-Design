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
/* Archivo: estadisticas.h								*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Descripción: Archivo de encabezado de la biblioteca	*/
/* estadisticas.c.										*/
/********************************************************/

#include "psomatrixcircuit.h"

#ifndef ESTADISTICAS
#define ESTADISTICAS

#define sizeCrom 10000			/*Numero maximo de carateres de un cromosoma*/
#define sizeNumDouble 100		/*Numero maximo de digitos que contendra un numero double*/	

/*Estructura para almacenar la informacion estadistica de una generacion o de una corrida del programa*/
typedef struct{				
	double mediaAptitud;		/*Aptitud media de las partículas*/
	double cuadradosAptitud;	/*Cuadrado de las aptitudes, para la posterior obtención de la varianza*/
	unsigned generacion;		/*Generacion en la que se encontró el óptimo*/
	particula mejor;			/*Mejor partícula encontrada*/
	particula peor;				/*Peor partícula encontrada*/
}estadistica;

estadistica Gen;				/*Informacion de la estadistica de la generacion*/
estadistica Corr;				/*Informacion de la estadistica de la corrida*/

/*Funciones utilizadas para la obtención de las estadisticas*/
void iniEstadistica(estadistica*);			/*Inicializa la informacion de la estadística de la corrida o de la generacion*/
void estEncGlobal(char*);					/*Genera el encabezado de la estadistica global de todas las corridas*/
void estCorr(char*,unsigned);				/*Genera la estadistica de cada corrida en el archivo global de corridas*/
void estEncCorr(char*,unsigned);			/*Genera el encabezado de la estadistica de la corrida*/
void estGen(char*,unsigned);				/*Genera la estadistica de una generaciones de la corrida*/
void estPieCorr(char*);						/*Genera los resultados finales de la estadistica de la corrida*/
void cadenaCrom(particula,char*);			/*Convierte el cromosoma de números reales o binario a una cadena de caracteres*/
void dtoa(double,unsigned,char*);			/*Convierte un número real en una cadena de caracteres*/

#endif
