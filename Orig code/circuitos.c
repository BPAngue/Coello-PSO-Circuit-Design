/********************************************************/
/*					CINVESTAV - IPN						*/
/*			Departamento de Ingenería Eléctrica			*/
/*					Sección Computación					*/
/*														*/	
/*				   Computacion Evolutiva				*/
/*														*/
/*					Erika Hernandez Luna				*/
/*			 eluna@computacion.cs.cinvestav.mx			*/
/*					 2 / agosto / 2003	  				*/
/*														*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Archivo: circuitos.c									*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Descripción: Biblioteca para cargar en memoria un	*/
/* archivo con el diseño del circuito lógico a optimizar*/
/********************************************************/

#include "circuitos.h"

/*Genera las entradas de la tabla de verdad del circuito*/
void generaTT(void)
{
	unsigned i, j;

	for(i=0; i<numRengTT; i++){
		entradaTT[i] = (unsigned *)malloc(sizeof(unsigned)*numEntradas);
		for(j=0; j<numEntradas; j++)
			entradaTT[i][j] = (i>>j)&1;
	}
}

/*Carga en memoria los datos de la tabla de verdad del archivo de entrada*/
void cargaTT(FILE *pf)
{
	int j;
	unsigned i;

	numRengTT = 1;
	numEntradas = (unsigned)leeNumero(pf);
	numRengTT <<= numEntradas;
	numSalidas = (unsigned)leeNumero(pf);
	numTSalidas = numSalidas * numRengTT;

	entradaTT = (unsigned **)malloc(sizeof(unsigned *)*numRengTT);
	salidaTT = (unsigned **)malloc(sizeof(unsigned *)*numRengTT);
	
	generaTT();

	for(i=0; i<numRengTT; i++){
		salidaTT[i] = (unsigned *)malloc(sizeof(unsigned)*numSalidas);
		for(j=0; (unsigned)j<numSalidas; j++)
			salidaTT[i][j] = (unsigned)leeNumero(pf);
	}
}

/*Despliega los datos de la tabla de de verdad obtenidos del archivo de entrada*/
void imprimeTT(void)
{
	int j;
	unsigned i;

	printf("      Numero de entradas de la tabla de verdad: %d\n",numEntradas);
	printf("       Numero de salidas de la tabla de verdad: %d\n",numSalidas);
	printf("     Numero de renglones de la tabla de verdad: %d\n",numRengTT);

	printf("\nTabla de verdad\n\n");
	for(i=0; i<numEntradas; i++) printf("E%d ",i);
	for(i=0; i<numSalidas; i++)	printf("S%d ",i);
	
	for(i=0; i<numRengTT; i++){
		printf("\n");
		for(j=numEntradas-1; j>=0; j--)
			printf(" %d ",entradaTT[i][j]);
		for(j=0; (unsigned)j<numSalidas; j++)
			if(salidaTT[i][j] == DONTCARE)
				printf(" * ");
			else	
				printf(" %d ",salidaTT[i][j]);
	}
	printf("\n\n");
}

/*Libera la memoria reservada para la tabla de verdad del circuito*/
void liberaMemoriaTT(void)
{
	unsigned i;

	for(i=0; i<numRengTT; i++){
		free(entradaTT[i]);
		free(salidaTT[i]);
	}
	free(entradaTT);
	free(salidaTT);
}

/*Lee un número del archivo apuntado por pf*/
double leeNumero(FILE *pf)
{
	char car, num[100];
	unsigned i = 0;

    while(!feof(pf) && !esDigito(car = fgetc(pf)) && car != '.' && car != ';');
	if(car == ';'){ 
		while(car != EndLine && !feof(pf)) car = fgetc(pf);
		car = fgetc(pf);
	}

	while(!feof(pf) && (car != EndLine && car!=' ' && car!=';')){
		if(esDigito(car) || car=='.')num[i++] = car;
		car = fgetc(pf);
	}
	if(car == ';') while((car = fgetc(pf)) != EndLine && !feof(pf));
	
	num[i] = EndString;
	return atof(num);
}

/*Lee una cadena del archivo apuntado por pf*/
void leeCadena(FILE *pf, char *cadena)
{
	unsigned i = 0;
	char car;

	while((car = fgetc(pf)) != ';' && car!=EndLine && !isalpha(car));
	if(car == ';'){
		while(car != EndLine && !feof(pf)) car = fgetc(pf);
		car = fgetc(pf);
	}
	
	while(!feof(pf) && (car != EndLine && car!=' ' && car!=';' && car!='\t')){
		cadena[i++] = car;
		car = fgetc(pf);
	}
	if(car == ';') while((car = fgetc(pf)) != EndLine && !feof(pf));
	cadena[i] = EndString;
}

/*Devuelve 1 en caso de que el caracter sea igual a 0-9, 0 en cualquier otro caso*/
unsigned esDigito(char a)
{
	return a=='0' || a=='1' || a=='2' || a=='3' || a=='4' || a=='5' || a=='6' || a== '7' || a=='8' || a=='9' ? 1 : 0;
}
