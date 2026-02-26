/********************************************************/
/*					CINVESTAV - IPN						*/
/*			Departamento de Ingenería Eléctrica			*/
/*					Sección Computación					*/
/*														*/	
/*				   Computacion Evolutiva				*/
/*														*/
/*					Erika Hernandez Luna				*/
/*			 eluna@computacion.cs.cinvestav.mx			*/
/*				      21 / enero / 2003					*/
/*														*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Archivo: random.h									*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Descripción: Archivo de encabezado de la biblioteca	*/
/* random.c.											*/
/********************************************************/

#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

#ifndef RANDOM_LIB
#define RANDOM_LIB

typedef struct eo{
	double valor;
	unsigned pos;
}elemordena;

elemordena *arrAuxOrdenamiento;			/*Arreglo auxiliar para el método de ordenamiento merge sort*/

/*Funciones de la biblioteca random*/
unsigned initrandom(unsigned);						/*Inicialización de la semilla para la generación de numeros aleatorios*/
double rndF(void);									/*Generacion de numeros aleatorios entre 0 y 1*/
unsigned rndI(unsigned);							/*Generacion de un numero entero aleatorio entre 0 y rng*/
double rndFR(double,double);						/*Generacion de un numero aleatorio real entre inf y sup*/
unsigned rndIR(unsigned,unsigned);					/*Genera un numero aleatorio positivo entre inf y sup*/
unsigned flip(double);								/*Devuelve true con probabilidad prob*/
void swap(int*,int*);								/*Intercambia el valor de dos números*/
void swapun(unsigned*,unsigned*);					/*Intercambia el valor de dos números sin signo*/
void barajeo(unsigned*,unsigned);					/*Barajeo de números*/
unsigned round(double);								/*Redondeo de números*/
int rndInt(int);									/*Genera un número entero no mayor que limite-1*/
int rndIntLimite(int,int);							/*Genera un número entero entre limInf y limSup-1*/
double sigmoide(double);							/*Funcion sigmoide que devuelve un valor entre 0 y 1 de acuerdo al parametro enviado*/
void algMergeSort(elemordena*,unsigned);			/*Algoritmo pricipal de merge sort, aqui se reserva la memoria necesaria para el arreglo auxiliar*/
void mergeSort(elemordena*,unsigned,unsigned);		/*Algoritmo merge sort para ordenamiento de datos*/
void merge(elemordena*,unsigned,unsigned,unsigned);	/*Función usada por el algoritmo merge sort que "mezcla" datos*/

#endif