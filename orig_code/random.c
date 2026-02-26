/********************************************************/
/*					CINVESTAV - IPN						*/
/*			Departamento de Ingenería Eléctrica			*/
/*					Sección Computación					*/
/*														*/	
/*				   Computacion Evolutiva				*/
/*														*/
/*					Erika Hernandez Luna				*/
/*			 eluna@computacion.cs.cinvestav.mx			*/
/*				     21 / enero / 2003					*/
/*														*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Archivo: random.c									*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Descripción: Biblioteca para la generacion de		*/
/* numeros aleatorios con propiedades específicas,		*/
/* y para el barajeo aleatorio de números.				*/
/********************************************************/

#include "random.h"

/*Inicialización de la semilla para la generación de
  numeros aleatorios*/
unsigned initrandom(unsigned sem)
{    
	unsigned semilla = !sem ? (unsigned)(time(NULL)) : sem;
	srand(semilla);
	return semilla;
}

/*Generacion de numeros aleatorios entre 0 y 1*/
double rndF(void)
{
	return ((double)rand())/RAND_MAX;
}

/*Generacion de un numero entero aleatorio entre 0 y rng*/
unsigned rndI(unsigned rng)
{
	double val;
	if((val = rndF()*(rng)-1)<0) return 0;
	
    return (unsigned)(val);
}

/*Generacion de un numero aleatorio real entre inf y sup*/
double rndFR(double inf, double sup)
{
	return inf + rndF()*(sup - inf);
}

/*Genera un numero aleatorio positivo entre inf y sup*/
unsigned rndIR(unsigned inf, unsigned sup)
{
	return (unsigned)(inf + rndF()*(sup - inf));
}


/*Devuelve true con probabilidad prob*/
unsigned flip(double prob)
{
    if(rndF() <= prob) return(1);
    return(0);
}

/*Intercambia el valor de dos números*/
void swap(int *a, int *b)
{
	int temp = *a;

	*a = *b;
	*b = temp;
}

/*Intercambia el valor de dos números sin signo*/
void swapun(unsigned *a, unsigned *b)
{
	int temp = *a;

	*a = *b;
	*b = temp;
}

/*Barajeo de números*/
void barajeo(unsigned *arr, unsigned tam)
{
	int i;
	
	for(i=tam-1; i>0; i--)
		swapun(&(arr[i]),&(arr[rndI(i)]));
}

/*Redondeo de números*/
unsigned round(double a)
{
	double w = fabs(a) - floor(fabs(a));
	w = w<0.5 ? floor(fabs(a)) : ceil(fabs(a));
	w = a<0 ? w*(-1) : w;
	return((int)w);
}

/*Genera un número entero no mayor que limite-1*/
int rndInt(int limite)
{	
	double div = ((double)rand() / RAND_MAX);
	return (int)(div * limite);
}

/*Genera un número entero entre limInf y limSup-1*/
int rndIntLimite(int limInf, int limSup)
{
	if(limInf>limSup) swap(&limInf,&limSup);
	return rndInt(limSup - limInf) + limInf;
}


/*Funcion sigmoide que devuelve un valor entre 0 y 1 de acuerdo al parametro enviado*/
double sigmoide(double val)
{
	return 1/(1 + exp(-val));
}

/*Algoritmo pricipal de merge sort, aqui se reserva la memoria necesaria para el arreglo auxiliar*/
void algMergeSort(elemordena *datos,unsigned tam)
{
	arrAuxOrdenamiento = (elemordena *)malloc(sizeof(elemordena)*tam);
	mergeSort(datos, 0, tam-1);
	free(arrAuxOrdenamiento);
}

/*Algoritmo merge sort para ordenamiento de datos*/
void mergeSort(elemordena *datos, unsigned inicio, unsigned fin)
{
	unsigned mitad;

	if(inicio < fin){
		mitad = (fin + inicio) / 2;
		mergeSort(datos, inicio, mitad);
		mergeSort(datos, mitad + 1, fin);
		merge(datos, inicio, mitad, fin);
	}
}

/*Función usada por el algoritmo merge sort que "mezcla" datos*/
void merge(elemordena *datos, unsigned inicio, unsigned mitad, unsigned fin)
{
	unsigned i,j = inicio, k = mitad + 1;

	for(i=inicio; i<=fin; i++){
		if((datos[j].valor >= datos[k].valor && j<=mitad) || k>fin){
			arrAuxOrdenamiento[i].valor = datos[j].valor;
			arrAuxOrdenamiento[i].pos = datos[j++].pos;
		}
		else{
			arrAuxOrdenamiento[i].valor = datos[k].valor;
			arrAuxOrdenamiento[i].pos = datos[k++].pos;
		}
	}
			
	for(i=inicio; i<=fin; i++){
		datos[i].valor = arrAuxOrdenamiento[i].valor;
		datos[i].pos = arrAuxOrdenamiento[i].pos;
	}

}