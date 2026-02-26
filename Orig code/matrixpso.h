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
/* matrixpso.c.											*/
/********************************************************/

#include <malloc.h>

#ifndef MATRIXPSO
#define MATRIXPSO

/*Define el tipo de datos de la matriz*/
#define matrix unsigned*	

/*Reserva memoria para almacenar las variables decodificadas de cada una de las celdas de la matriz*/
#define memValCelda (unsigned)malloc(sizeof(unsigned)*numReng*numCols)	

/*Tamaño maximo de la expresion booleana y su tipo de dato*/
#define sizeExpBooleana	4000			
typedef char cadBooleana[sizeExpBooleana];

unsigned numReng;		/*Numero de renglones de la matriz*/
unsigned numCols;		/*Numero de columnas de la matriz*/
unsigned *entrada1;		/*Valores decodificado de las entradas 1 de las compuertas para obtener la aptitud de la matriz (temporal)*/
unsigned *entrada2;		/*Valores decodificado de las entradas 2 de las compuertas para obtener la aptitud de la matriz (temporal)*/
unsigned *tipoG;		/*Valores decodificado del tipo de compuerta para obtener la aptitud de la matriz (temporal)*/
unsigned **inTT;		/*Almacena entrada y salida de cada etapa de evaluacion de la matriz de compuertas (temporal)*/
unsigned *cuentaG;		/*Almacena la cuenta de las compuertas de una matriz (temporal)*/
unsigned *salida;		/*Almacena la salida de una sola etapa de evaluacion  de la matriz de compuertas (temporal)*/
unsigned tMat;			/*Tamaño de la matriz usada = numReng*numCols*/


/*Funciones para el manejo de la matriz*/
void reservaMemoriaMatriz(void);				/*Reserva memoria para las variables globales de la decodificacion*/
void liberaMemoriaMatriz(void);					/*Libera memoria para las variables globales de la decodificacion*/
void decodifica(matrix);						/*Convierte los valores de cada una de las celdas de la matriz en las entradas1, entradas2 y tipo de compuerta*/
void evalua(matrix,unsigned*);					/*Evalua a la matriz y obtiene el número de violaciones al circuito pedido asi como el número de WIRE's*/
unsigned cuentaCompuertas(void);				/*Cuenta el número de compuertas totales involucradas en las soluciones*/
void cCompuerta(unsigned,int,unsigned*);		/*Cuenta el numero de compuertas involucradas en una solucion*/
void expresion(matrix,cadBooleana,unsigned);	/*Obtiene la expresión booleana que representa la matriz*/
void cadenaBooleana(unsigned,int,cadBooleana);	/*Procedimiento recursivo que va pegando a la cadena los caracteres que definirán la exp booleana*/
void imprimeMatriz(matrix M);					/*Imprime los datos de la matriz de compuertas*/

#endif