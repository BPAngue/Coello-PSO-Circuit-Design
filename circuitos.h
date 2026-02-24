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
/* Archivo: circuitos.h									*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Descripción: Archivo de encabezado de la biblioteca	*/
/* circuitos.c.											*/
/********************************************************/

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#ifndef CIRCUITOS
#define CIRCUITOS

#define EndLine '\n'
#define EndString '\x0'

/*Funciones que pueden estar contenidas en el arbol*/
#define AND   0
#define OR    1
#define NOT   2
#define WIRE  3
#define XOR   4
#define NOT1  5
#define WIRE1 6
#define XOR1  7

#define DONTCARE 2			/*Condicion de no importa de la tabla de verdad de salida*/

unsigned numGates;			/*Numero de compuertas disponibles*/
unsigned numEntradas;		/*Numero de entradas en la tabla de verdad*/
unsigned numSalidas;		/*Numero de salidas definidas en la tabla de verdad*/
unsigned numRengTT;			/*Numero de renglones en la tabla de verdad*/
unsigned numTSalidas;		/*Numero total de salidas que debe cumplir el circuito, numSalidas*numRengTT*/
unsigned **entradaTT;		/*Entradas en la tabla de verdad*/
unsigned **salidaTT;		/*Salidas de la tabla de verdad*/

void generaTT(void);				/*Genera las entradas de la tabla de verdad del circuito*/
void cargaTT(FILE *);				/*Carga en memoria los datos de la tabla de verdad del archivo de entrada*/
void imprimeTT(void);				/*Despliega los datos de la tabla de de verdad obtenidos del archivo de entrada*/
void liberaMemoriaTT(void);			/*Libera la memoria reservada para la tabla de verdad del circuito*/
double leeNumero(FILE *);			/*Lee un número del archivo apuntado por pf*/
void leeCadena(FILE *pf,char*);		/*Lee una cadena del archivo apuntado por pf*/
unsigned esDigito(char);			/*Devuelve 1 en caso de que el caracter sea igual a 0-9, 0 en cualquier otro casi*/

#endif
