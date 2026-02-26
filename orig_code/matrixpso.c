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
/* Descripción: Biblioteca para llevar a cabo el manejo	*/
/* de la matriz que representa a un circuito			*/
/********************************************************/

#include <math.h>
#include <string.h>
#include "matrixpso.h"
#include "circuitos.h"
#include "psomatrixcircuit.h"

/*Reserva memoria para las variables globales de la decodificacion*/
void reservaMemoriaMatriz(void)
{
	unsigned i;

	entrada1 = memVarX;
	entrada2 = memVarX;
	tipoG = memVarX;
	salida = (unsigned *)malloc(sizeof(unsigned)*numReng);
	cuentaG = (unsigned *)malloc(sizeof(unsigned)*tMat);
	inTT = (unsigned **)malloc(sizeof(unsigned *)*numRengTT);

	for(i=0; i<numRengTT; i++)
		inTT[i] = (unsigned *)malloc(sizeof(unsigned)*numReng);
}

/*Libera memoria para las variables globales de la decodificacion*/
void liberaMemoriaMatriz(void)
{
	unsigned i;

	free(entrada1);
	free(entrada2);
	free(tipoG);
	free(salida);
	free(cuentaG);
	
	for(i=0; i<numRengTT; i++)
		inTT[i] = (unsigned *)malloc(sizeof(unsigned)*numReng);
	free(inTT);
}

/*Decodifica la matriz en los valores de las variables*/
void decodifica(matrix M)
{
	int k;
	unsigned i,j,num,in;

	for(i=0,j=0,in=0; i<nVar; i++){

		if(representacion != BINARIA)
			num = M[i];
		else{
			num = 0;
                          			for (k=bitVariable[i]-1; k>=0; k--,j++)
				if(M[j]) num += (unsigned)pow(2,k);
		}

		switch(i%3){
			case 0:	/*Entrada 1*/
				entrada1[in] = (num%numReng);
				break;
			case 1:	/*Entrada 2*/
				entrada2[in] = (num%numReng);
				break;
			case 2:	/*Tipo*/
				tipoG[in++] = (num%numGates);
		}
	}
}

/*Evalua a la matriz y obtiene el número de violaciones al circuito pedido asi como el número de WIRE's*/
void evalua(matrix M, unsigned *numIgual)
{
	unsigned i,j,k,in;

	*numIgual = 0;
	decodifica(M);

	for(i=0; i<numRengTT; i++)
		for(j=0; j<numReng; j++)
			inTT[i][j] = j<numEntradas ? entradaTT[i][j] : inTT[i][j-numEntradas];	

	for(i=0; i<numRengTT; i++){
		for(j=0,in=0; j<tMat; j++){
			switch(tipoG[j]) {
				case AND : 
					salida[in] = inTT[i][entrada1[j]] & inTT[i][entrada2[j]];
					break;
				case OR : 
					salida[in] = inTT[i][entrada1[j]] | inTT[i][entrada2[j]];
					break;
				case NOT : 
				case NOT1:
					salida[in] = inTT[i][entrada1[j]] ? 0 : 1;
					break;
				case WIRE :
				case WIRE1:
					salida[in] = inTT[i][entrada1[j]];
					break;
				case XOR : 
				case XOR1 : 
					salida[in] = inTT[i][entrada1[j]] ^ inTT[i][entrada2[j]];
					break;
			}

			if(!((j+1) % numReng)){
				in = 0;
		        for (k=0; k<numReng; k++) 
					inTT[i][k] = salida[k];
			}
			else in++;
		}

		for(j=0; j<numSalidas;j++)
			if(salida[j] == salidaTT[i][j]) (*numIgual)++; 
	}
}

/*Cuenta el número de compuertas totales involucradas en las soluciones*/
unsigned cuentaCompuertas(void)
{
	unsigned i, compuertas = 0; 
	unsigned salidaI = numReng*(numCols-1);
	
	for(i=0; i<tMat; i++) cuentaG[i] = 0;
	for(i=salidaI; i<salidaI+numSalidas; i++)
		cCompuerta(i,numCols-1,cuentaG);
	for(i=0; i<tMat; i++) if(cuentaG[i]) compuertas++;

	return compuertas;
}

/*Cuenta el numero de compuertas involucradas en una solucion*/
void cCompuerta(unsigned salida, int columna, unsigned *cuentaG)
{
	if(columna>=0 && salida>=0){
		if(tipoG[salida] != WIRE && tipoG[salida] != WIRE1)
			if(tipoG[salida] != AND && tipoG[salida] != OR)
			   cuentaG[salida] = 1;
			else
				if(entrada1[salida] != entrada2[salida])
					cuentaG[salida] = 1;

		cCompuerta(numReng*(columna-1) + entrada1[salida],columna-1,cuentaG);
		cCompuerta(numReng*(columna-1) + entrada2[salida],columna-1,cuentaG);  
	}
}

/*Obtiene la expresión booleana que representa la matriz*/
void expresion(matrix M, cadBooleana expBool, unsigned salida)
{
	strcpy(expBool,"");
	decodifica(M);
	cadenaBooleana(numReng*(numCols-1) + salida,numCols-1,expBool);
}

/*Procedimiento recursivo que va pegando a la cadena los caracteres que definirán la exp booleana*/
void cadenaBooleana(unsigned celda, int col, cadBooleana cadena)
{
	unsigned nuevae1, nuevae2;
	char chrgate[7], chrvar[2];

	chrvar[1] = EndString;
	switch(tipoG[celda]){
			case AND:
				strcpy(chrgate,"(AND1 ");
				break;
			case OR:
				strcpy(chrgate,"(OR1 ");
				break;
			case NOT:
			case NOT1:
				strcpy(chrgate,"(NOT1 ");
				break;
			case XOR:
			case XOR1:
				strcpy(chrgate,"(XOR1 ");
				break;
			case WIRE:
			case WIRE1:
				strcpy(chrgate,"(WIRE ");
				break;
	}

	strcat(cadena,chrgate);
	if(col<=0){
		chrvar[0] = entrada1[celda] % numEntradas + (int)'A';
		strcat(cadena,chrvar);
		if(tipoG[celda] != WIRE && tipoG[celda] != WIRE1 && tipoG[celda] != NOT && tipoG[celda] != NOT){
			chrvar[0] = entrada2[celda] % numEntradas + (int)'A';
			strcat(cadena," ");
			strcat(cadena,chrvar);
		}
	}
	else{
		nuevae1 = numReng*(col-1) + entrada1[celda];
		nuevae2 = numReng*(col-1) + entrada2[celda];

		cadenaBooleana(nuevae1,col-1,cadena);
		if(tipoG[celda] != WIRE && tipoG[celda] != WIRE1 && tipoG[celda] != NOT && tipoG[celda] != NOT){
			strcat(cadena," ");
			cadenaBooleana(nuevae2,col-1,cadena);
		}
	}
	strcat(cadena,")");
}

/*Imprime los datos de la matriz de compuertas*/
void imprimeMatriz(matrix M)
{
	int dec,sign;
	unsigned i, j;
	char cadena[5000];
	
	printf("\n");
	decodifica(M);
	for(i=0; i<numReng; i++){
		strcpy(cadena,"");
		for(j=0; j<numCols; j++){
			switch(tipoG[j*numReng + i]) {
					case AND : 
						strcat(cadena,"AND(");
						break;
					case OR: 
						strcat(cadena,"OR(");
						break;
					case NOT : 
						strcat(cadena,"NOT(");
						break;
					case WIRE : 
						strcat(cadena,"WIRE(");
						break;
					case XOR : 
						strcat(cadena,"XOR(");
						break;
					case NOT1:
						strcat(cadena,"NOT(");
						break;
					case WIRE1:
						strcat(cadena,"WIRE(");
						break;
					case XOR1 : 
						strcat(cadena,"XOR(");
						break;

			}

			strcat(cadena, entrada1[j*numReng + i] == 0 ? "0" : fcvt(entrada1[j*numReng + i],0,&dec,&sign));
			strcat(cadena, " ");
			strcat(cadena, entrada2[j*numReng + i] == 0 ? "0" : fcvt(entrada2[j*numReng + i],0,&dec,&sign));
			strcat(cadena,")");
			if(j!=numCols-1) strcat(cadena,",");
		}
		printf("%s\n",cadena);
	}
}
