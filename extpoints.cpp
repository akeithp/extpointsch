/**
 * Magister en Informatica
 * Universidad Austral de Chile
 * Computacion de Alto Rendimiento
 * Prof. Dr. Héctor Ferrada
 * 
 * Alan Keith Paz
 * Due Date: 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <random>

using namespace std;

#define PRINT 0
#define TEST 1
//#define NTHREADS 4

uint NTHREADS;	// number of threads
ulong REPET = 1;
ulong REPETTEST = 1;

//ELIMINAR
/*struct point{
	double x;
	double y;
};*/

float xup, xdown, xleft, xright, xri_up, xup_le, xle_do, xdo_ri;
float yup, ydown, yleft, yright, yri_up, yup_le, yle_do, ydo_ri;

typedef struct{
	float *X, *Y;	// array of float points X and Y
	bool NORMAL;	// flag probability: 1 = NORMAL, 0 = UNIFORM
	float sigma;	// for normal distribution probability function
	float mean;		// for normal distribution probability function
	ulong n;
	float minx;		// for uniform distribution
	float maxx;		// for uniform distribution 
	float miny;		// for uniform distribution 
	float maxy;		// for uniform distribution 
	//point mean;
	/*point up;
	point down;
	point right;
	point left;
	point *Q1;	// cuadrante 1
	point *Q2;	// cuadrante 2
	point *Q3;	// cuadrante 3
	point *Q4;	// cuadrante 4*/

} pointSet;

void genArrays(pointSet *points);
void divideArrays(pointSet *points);
void runEPS(pointSet *points);
void runEPSP(pointSet *points);
void filteredPercent(pointSet *points);
void testPoints(pointSet *points); //comprobar correctitud de los puntos encontrados (secuencial vs paralelo)

int main(int argc, char const *argv[]){
	double t1,t2;
	float avgTime;
	char aFile[400];
	char str[100];
	uint i, j;
	bool flagP;
	pointSet *points = new pointSet();
	flagP = atoi(argv[1]);
	points->n = atoi(argv[2]);
	points->NORMAL = atoi(argv[3]);
	if(points->NORMAL){
		points->mean = atof(argv[4]);
		points->sigma = atof(argv[5]);
	}
	else{
		points->minx = atof(argv[4]);
		points->maxx = atof(argv[5]);
		points->miny = atof(argv[6]);
		points->maxy = atof(argv[7]);
	}
	/*
	 * Parallel Section
	 * @argv[6] number of threads
	 */
	if(flagP){
		if(points->NORMAL) 
			NTHREADS = atoi(argv[6]);
		else 
			NTHREADS = atoi(argv[8]);
		omp_set_num_threads(NTHREADS);
	}
	
	//REPET
	//CARPETA DE RESULTADOS
	
	//REPETTEST = atol(argv[7]);
	// OJO: Corregir comprobación de parámetros de entrada y mensaje de error
	/*if(flagP && points->NORMAL && argc < 7){
		cout << "Execution Error! call: ./prog <PARALLEL> <n> <NORMAL flag> [<mean>] [<sigma>] <threads> <REPETS FOR TEST" << endl;
		exit(EXIT_FAILURE);
	}*/
	if(PRINT){
		cout << "Parameters..." << endl;
		cout << "n = " << points->n << endl;
		cout << "NORMAL flag = " << points->NORMAL << endl;
		cout << "Number of Threads = " << NTHREADS << endl;
		if (points->NORMAL){
			cout << "mean = " << points->mean << endl;
			cout << "sigma = " << points->sigma << endl;
		}
	}
	//float t; // = clock();
	
	genArrays(points);
	
	//t = (clock() - t)/CLOCKS_PER_SEC;	// seconds
	//cout << "Construction Time = " << t << " Seconds" << endl;
	
	/*
	 * Serial version
	 */
	 if(!flagP){
		//t = clock();
		for(i=0; i<REPETTEST; i++){
			t1=omp_get_wtime(); 
			for(j=0; j<REPET; j++){
				runEPS(points);
			}
			t2=omp_get_wtime();
			avgTime = (t2 - t1)/REPET;
			//t = (clock() - t)/CLOCKS_PER_SEC;	// seconds
			//cout << "Serial Time = " << t/REPET << " Seconds" << endl;
			if(PRINT){
				cout << "Average CPU time per execution, Secuential = " << avgTime << endl; //revisar unidad tiempos
				cout << "Primary Extreme Points: UP(" << xup << "," << yup << ")" << "; DOWN(" << xdown << "," << ydown <<"); LEFT(" << xleft << "," << yleft << "); RIGHT(" << xright << "," << yright << ")" << endl;
				cout << "Secondary Extreme Points: RIGHT_UP(" << xri_up <<"," << yri_up << "); UP_LEFT(" << xup_le << "," << yup_le << "); LEFT_DOWN(" << xle_do << "," << yle_do << "); DOWN_RIGHT(" << xdo_ri << "," << ydo_ri << ")" << endl;
			}
			strcpy(aFile, "./RESULTS/");
			strcpy(str, "");
			sprintf(str, "EPS%i", points->NORMAL); //%ld, points->n
			strcat(aFile, str);
			cout << "Resume File: " << aFile << endl;
			FILE *fp = fopen(aFile, "a+" );
			if (points->NORMAL){
				// [n] [REPETTEST] [avg eps-time/exec] [esperanza] [varianza]
				fprintf(fp, "%ld %ld %f %f %f\n", points->n, REPETTEST, (avgTime*1000000.0), points->mean, points->sigma);
			}else{
				// [n] [REPETTEST] [nOcc] [avg bs-time/exec]
				fprintf(fp, "%ld %ld %f\n", points->n, REPETTEST, (avgTime*1000000.0));
			}
			fclose(fp);
		}
	}
	
	/*
	 * Parallel version
	 */
	 if(flagP){
		 for(i=0; i<REPETTEST; i++){
			t1=omp_get_wtime(); 
			for(j=0; j<REPET; j++){
				runEPSP(points);
			}
			t2=omp_get_wtime();
			avgTime = (t2 - t1)/REPET;
			if(PRINT){
				cout << "Average CPU time per execution, Parallel = " << avgTime << endl;
				cout << "Primary Extreme Points: UP(" << xup << "," << yup << ")" << "; DOWN(" << xdown << "," << ydown <<"); LEFT(" << xleft << "," << yleft << "); RIGHT(" << xright << "," << yright << ")" << endl;
				cout << "Secondary Extreme Points: RIGHT_UP(" << xri_up <<"," << yri_up << "); UP_LEFT(" << xup_le << "," << yup_le << "); LEFT_DOWN(" << xle_do << "," << yle_do << "); DOWN_RIGHT(" << xdo_ri << "," << ydo_ri << ")" << endl;
			}
			strcpy(aFile, "./RESULTSP/");
			strcpy(str, "");
			sprintf(str, "EPSParallel%i%d", points->NORMAL, NTHREADS); //%ld, points->n
			strcat(aFile, str);
			cout << "Resume File: " << aFile << endl;
			FILE *fp = fopen(aFile, "a+" );
			if (points->NORMAL){
				// [n] [REPETTEST] [avg epsP-time/exec] [esperanza] [varianza]
				fprintf(fp, "%ld %ld %f %f %f\n", points->n, REPETTEST, (avgTime*1000000.0), points->mean, points->sigma);
			}else{
				// [n] [REPETTEST] [nOcc] [avg bs-time/exec]
				fprintf(fp, "%ld %ld %f\n", points->n, REPETTEST, (avgTime*1000000.0));
			}
			fclose(fp);
			if(TEST){
				testPoints(points);
			}
		}
		
	}
	filteredPercent(points);
	return 0;
}

// generar arrays
void genArrays(pointSet *points){
	ulong i;
	float numX, numY;

	if (points->NORMAL){ // DISTRIBUCION NORMAL
		points->X = new float[points->n]; // +1?
		points->Y = new float[points->n]; // +1?
		default_random_engine generator;
		normal_distribution<double> distribution(points->mean, points->sigma);	// (mean, stddev)

		for (i=0; i<points->n; i++){
			numX = distribution(generator);
	    	numY = distribution(generator);
	    	points->X[i] = numX;
	    	points->Y[i] = numY;
		}
	}else{ // DISTRIBUCION UNIFORME
		points->X = new float[points->n]; // +1?
		points->Y = new float[points->n]; // +1?
		//points->X[0] = (float)rand()/RAND_MAX*(maxi-mini)+mini; //parametro de entrada, definir
		//points->Y[0] = (float)rand()/RAND_MAX*(maxi-mini)+mini;
		float varx = points->maxx - points->minx;
		float vary = points->maxy - points->miny;
		for (i=0; i<points->n; i++){
			points->X[i] = (float)rand()/RAND_MAX*varx+points->minx;
			points->Y[i] = (float)rand()/RAND_MAX*vary+points->miny;
		}
	}

}

/**
 * Divide los arrays por cuadrante
 */
void divideArrays(pointSet *points){
	//ulong i, q1, q2, q3, q4;
	//q1 = q2 = q3 = q4 = 0;
	/*for(i=0; i<points->n;i++){
		if( (points->X[i].x >= 0) && (points->X[i].y >= 0) ){
			points->Q1[q1] = points[i];
			q1++;
		}
		if( (points->X[i].x >= 0) && (points->X[i].y < 0) ){
			points->Q2[q2] = points[i];
			q2++;
		}
		if( (points->X[i].x < 0) && (points->X[i].y < 0) ){
			points->Q3[q3] = points[i];
			q3++;
		}
		if( (points->X[i].x < 0) && (points->X[i].y >= 0) ){
			points->Q4[q4] = points[i];
			q4++;
		}			
	}*/
}

/**
 *  Busca los 4 puntos extremos de un array
 */
void runEPS(pointSet *points){
	ulong i;
	xup = xdown = xleft = xright = points->X[0]; // xri_up = xup_le = xle_do = xdo_ri = points->X[0];
	yup = ydown = yleft = yright = points->Y[0]; // yri_up = yup_le = yle_do = ydo_ri = points->Y[0];
	bool flagRU, flagUL, flagLD, flagDR; // flags para saber si se encontró el pto extremo secundario
	//points->left = points->right = points->up = points->down = points->X[0];
	
	for(i=1;i<points->n;i++){
		// Find the rightest point
		if(points->X[i] >= xright){
		// Qué pasa cuando hay más de un punto con la misma coordenada máxima?
			if(points->X[i] != xright){
				xright = points->X[i];
				yright = points->Y[i];
			}
			else{
				//Calcular distancia euclidiana a la media y elegir el con mayor distancia?
			}
		}
		// Find the leftest point
		if(points->X[i] < xleft){	// Qué pasa cuando hay más de un punto con la misma coordenada mínima?
			xleft = points->X[i];
			yleft = points->Y[i];
		}
		// Find the upest point
		if(points->Y[i] > yup){	// Qué pasa cuando hay más de un punto con la misma ordenada máxima?
			yup = points->Y[i];
			xup = points->X[i];
		}
		// Find downest point
		if(points->Y[i] < ydown){	// Qué pasa cuando hay más de un punto con la misma ordenada mínima?
			ydown = points->Y[i];
			xdown = points->X[i];
		}
	}
	// Calcular puntos extremos secundarios
	
	// Creo que es mucho más eficiente ordenar los puntos primero por coordenada x y así no evaluar todos los puntos?
	
	// Cuadrante izquierdo superior
	float a, b, dist, denom;
	a = (yup-yleft)/(xup-xleft);
	b = -(yup-yleft)/(xup-xleft) + yleft;
	denom = sqrt(a*a+1);
	dist = 0;
	for(i=0; i<points->n; i++){ // podemos mejorarlo haciendo que solo corra hasta cuando x es mayor que el centro.
		if(points->X[i] < ((points->Y[i]-b)/a) && points->Y[i] > (a*points->X[i]+b)){
			float disti = abs(a*points->X[i]-points->Y[i]+b)/denom;
			if(disti > dist){
				xup_le = points->X[i];
				yup_le = points->Y[i];
				dist = disti;
				flagUL = true;
			}
		}
	}
	// Si no hay puntos sobre la recta lo igualamos al pto extremo izquierdo
	if(!flagUL){
		xup_le = xleft;
		yup_le = yleft;
	}
	
	// Cuadrante izquiero inferior
	a = (ydown-yleft)/(xdown-xleft);
	b = -(ydown-yleft)/(xdown-xleft) + yleft;
	denom = sqrt(a*a+1);
	dist = 0;
	for(i=0; i<points->n; i++){ // podemos mejorarlo haciendo que solo corra hasta cuando x es mayor que el centro.
		if(points->X[i] < ((points->Y[i]-b)/a) && points->Y[i] < (a*points->X[i]+b)){
			float disti = abs(a*points->X[i]-points->Y[i]+b)/denom;
			if(disti > dist){
				xle_do = points->X[i];
				yle_do = points->Y[i];
				dist = disti;
				flagLD = true;
			}
		}
	}
	// Si no hay puntos bajo la recta lo igualamos al pto extremo inferior
	if(!flagLD){
		xle_do = xdown;
		yle_do = ydown;
	}
	
	// Cuadrante inferior derecho
	a = (ydown-yright)/(xdown-xright);
	b = -(ydown-yright)/(xdown-xright) + yright;
	denom = sqrt(a*a+1);
	dist = 0;
	for(i=0; i<points->n; i++){ // podemos mejorarlo haciendo que solo corra hasta cuando x es mayor que el centro.
		if(points->X[i] > ((points->Y[i]-b)/a) && points->Y[i] < (a*points->X[i]+b)){
			float disti = abs(a*points->X[i]-points->Y[i]+b)/denom;
			if(disti > dist){
				xdo_ri = points->X[i];
				ydo_ri = points->Y[i];
				dist = disti;
				flagDR = true;
			}
		}
	}
	// Si no hay puntos bajo la recta lo igualamos al pto extremo derecho
	if(!flagDR){
		xdo_ri = xright;
		ydo_ri = yright;
	}
	
	// Cuadrante superior derecho
	a = (yup-yright)/(xup-xright);
	b = -(yup-yright)/(xup-xright) + yright;
	denom = sqrt(a*a+1);
	dist = 0;
	for(i=0; i<points->n; i++){ // podemos mejorarlo haciendo que solo corra hasta cuando x es mayor que el centro.
		if(points->X[i] > ((points->Y[i]-b)/a) && points->Y[i] > (a*points->X[i]+b)){
			float disti = abs(a*points->X[i]-points->Y[i]+b)/denom;
			if(disti > dist){
				xri_up = points->X[i];
				yri_up = points->Y[i];
				dist = disti;
				flagRU = true;
			}
		}
	}
	// Si no hay puntos sobre la recta lo igualamos al pto extremo superior
	if(!flagRU){
		xri_up= xup;
		yri_up = yup;
	}
}

/** 
 * Búsqueda de puntos extremos en paralelo
 */
void runEPSP(pointSet *points){
	ulong i, j;

	uint nnew = (int)points->n/NTHREADS;
	float xupi[NTHREADS], xdowni[NTHREADS], xlefti[NTHREADS], xrighti[NTHREADS];
	float yupi[NTHREADS], ydowni[NTHREADS], ylefti[NTHREADS], yrighti[NTHREADS];
	// Ejecutamos los for en paralelo y/o tomar subarreglos n/NTHREADS
	#pragma omp parallel shared(xupi, xdowni, xlefti, xrighti, yupi, ydowni, ylefti, yrighti)
	{
		#pragma omp for private(i,j)
		for(i=0; i < NTHREADS; i++){ //revisar para casos en que la división de la cantidad de puntos no es entera!
			xrighti[i] = xdowni[i] = xlefti[i] = xupi[i] = points->X[nnew*i];
			yrighti[i] = ydowni[i] = ylefti[i] = yupi[i] = points->Y[nnew*i];
			for(j=nnew*i+1; j < (i+1)*nnew; j++){
				{
					// Find the rightest point for the i-segment
					if(points->X[j] > xrighti[i]){
						//cout << "xrighti[ " << i << "] =" << xrighti[i]<< endl;
						xrighti[i] = points->X[j];
						yrighti[i] = points->Y[j];
						//cout << "X[" << j << "] = " << points->X[j] << endl;
						//cout << "Es el " << i << "mayor: " << points->X[j] << endl;
					}
					// Find the leftest point for the i-segment
					if(points->X[j] < xlefti[i]){
						xlefti[i] = points->X[j];
						ylefti[i] = points->Y[j];
					}
					// Find the upest point for the i-segment
					if(points->Y[j] > yupi[i]){	
						yupi[i] = points->Y[j];
						xupi[i] = points->X[j];
					}
					// Find downest point for the i-segment
					if(points->Y[j] < ydowni[i]){
						ydowni[i] = points->Y[j];
						xdowni[i] = points->X[j];
					}
				}
			}
		}
		#pragma omp barrier
		
		#pragma omp single
		{
			xright = xrighti[0];
			xleft = xlefti[0];
			yup = yupi[0];
			ydown = ydowni[0];
			for(i=1; i<NTHREADS; i++){
				//cout << xrighti[i] << endl;
				//#pragma omp sections
				//{
					//#pragma omp section
					if(xrighti[i] > xright){
						xright = xrighti[i];
						yright = yrighti[i];
					}
					//#pragma omp section
					if(xlefti[i] < xleft){
						xleft = xlefti[i];
						yleft = ylefti[i];
					}
					//#pragma omp section
					if(yupi[i] > yup){
						xup = xupi[i];
						yup = yupi[i];
					}
					//#pragma omp section
					if(ydowni[i] < ydown){
						xdown = xdowni[i];
						ydown = ydowni[i];
					}
				//}
			}
		}
		// PRIMERA PARTE, Búsqueda de puntos extremos debe ser subdividido el arreglo	
	}

	
	//SEGUNDA PARTE, Búsqueda de puntos extremos secundarios debe ser con for paralelos
	
	// Calcular puntos extremos secundarios
	
	// Creo que es mucho más eficiente ordenar los puntos primero por coordenada x y así no evaluar todos los puntos?
	

	//xri_up = xup_le = xle_do = xdo_ri = points->X[0];
	//yri_up = yup_le = yle_do = ydo_ri = points->Y[0];
	bool flagRU, flagUL, flagLD, flagDR; // flags para saber si se encontró el pto extremo secundario
	
	float a, b, dist, denom;
	
	#pragma omp parallel
	{
		omp_set_dynamic(1);
		#pragma omp sections private(i, a, b, denom, dist)
		{
			#pragma omp section
			{
				// Cuadrante izquierdo superior
				a = (yup-yleft)/(xup-xleft);
				b = -(yup-yleft)/(xup-xleft) + yleft;
				denom = sqrt(a*a+1);
				dist = 0;
				for(i=1; i<points->n; i++){ // podemos mejorarlo haciendo que solo corra hasta cuando x es mayor que el centro.
					if(points->X[i] < ((points->Y[i]-b)/a) && points->Y[i] > (a*points->X[i]+b)){
						float disti = abs(a*points->X[i]-points->Y[i]+b)/denom;
						if(disti > dist){
							xup_le = points->X[i];
							yup_le = points->Y[i];
							dist = disti;
							flagUL = true;
						}
					}
				}
						}
			#pragma omp section
			{
				// Cuadrante izquiero inferior
				a = (ydown-yleft)/(xdown-xleft);
				b = -(ydown-yleft)/(xdown-xleft) + yleft;
				denom = sqrt(a*a+1);
				dist = 0;
				for(i=1; i<points->n; i++){ // podemos mejorarlo haciendo que solo corra hasta cuando x es mayor que el centro.
					if(points->X[i] < ((points->Y[i]-b)/a) && points->Y[i] < (a*points->X[i]+b)){
						float disti = abs(a*points->X[i]-points->Y[i]+b)/denom;
						if(disti > dist){
							xle_do = points->X[i];
							yle_do = points->Y[i];
							dist = disti;
							flagLD = true;
						}
					}
				}
			}
			#pragma omp section
			{
				// Cuadrante inferior derecho
				a = (ydown-yright)/(xdown-xright);
				b = -(ydown-yright)/(xdown-xright) + yright;
				denom = sqrt(a*a+1);
				dist = 0;
				for(i=1; i<points->n; i++){ // podemos mejorarlo haciendo que solo corra hasta cuando x es mayor que el centro.
					if(points->X[i] > ((points->Y[i]-b)/a) && points->Y[i] < (a*points->X[i]+b)){
						float disti = abs(a*points->X[i]-points->Y[i]+b)/denom;
						if(disti > dist){
							xdo_ri = points->X[i];
							ydo_ri = points->Y[i];
							dist = disti;
							flagDR = true;
						}
					}
				}
			}
			#pragma omp section
			{
				// Cuadrante superior derecho
				a = (yup-yright)/(xup-xright);
				b = -(yup-yright)/(xup-xright) + yright;
				denom = sqrt(a*a+1);
				dist = 0;
				for(i=1; i<points->n; i++){ // podemos mejorarlo haciendo que solo corra hasta cuando x es mayor que el centro.
					if(points->X[i] > ((points->Y[i]-b)/a) && points->Y[i] > (a*points->X[i]+b)){
						float disti = abs(a*points->X[i]-points->Y[i]+b)/denom;
						if(disti > dist){
							xri_up = points->X[i];
							yri_up = points->Y[i];
							dist = disti;
							flagRU = true;
						}
					}
				}
			}
			
		}
		#pragma omp barrier
		// Si no hay puntos sobre la recta lo igualamos al pto extremo izquierdo
		if(!flagUL){
			xup_le = xleft;
			yup_le = yleft;
		}
		// Si no hay puntos sobre la recta lo igualamos al pto extremo inferior
		if(!flagLD){
			xle_do = xdown;
			yle_do = ydown;
		}
		// Si no hay puntos sobre la recta lo igualamos al pto extremo derecho
		if(!flagDR){
			xdo_ri = xright;
			ydo_ri = yright;
		}
		// Si no hay puntos sobre la recta lo igualamos al pto extremo superior
		if(!flagRU){
			xri_up = xup;
			yri_up = yup;
		}
	}
}

/** 
 * Cálculo de Porcentaje de Puntos Filtrados
 */
 void filteredPercent(pointSet *points){
	ulong i;
	ulong count = 0;
	float percent = 0;
	char aFile[400];
	char str[100];
	float a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, b3, b4, b5, b6,b7, b8;
	a1 = (yup-yup_le)/(xup-xup_le);
	b1 = -(yup-yup_le)/(xup-xup_le) + yup_le;
	a2 = (yleft-yup_le)/(xleft-xup_le);
	b2 = -(yleft-yup_le)/(xleft-xup_le) + yup_le;
	a3 = (yleft-yle_do)/(xleft-xle_do);
	b3 = -(yleft-yle_do)/(xleft-xle_do) + yle_do;
	a4 = (ydown-yle_do)/(xdown-xle_do);
	b4 = -(ydown-yle_do)/(xdown-xle_do) + yle_do;
	a5 = (ydown-ydo_ri)/(xdown-xdo_ri);
	b5 = -(ydown-ydo_ri)/(xdown-xdo_ri) + ydo_ri;
	a6 = (yright-ydo_ri)/(xright-xdo_ri);
	b6 = -(yright-ydo_ri)/(xright-xdo_ri) + ydo_ri;
	a7 = (yright-yri_up)/(xright-xri_up);
	b7 = -(yright-yri_up)/(xright-xri_up) + yri_up;
	a8 = (yup-yri_up)/(xup-xri_up);
	b8 = -(yup-yri_up)/(xup-xri_up) + yri_up;
	for(i=0; i<points->n; i++){
		if(points->X[i] < ((points->Y[i]-b1)/a1) && points->Y[i] > (a1*points->X[i]+b1)){
			count++;
		}
		else if(points->X[i] < ((points->Y[i]-b2)/a2) && points->Y[i] > (a2*points->X[i]+b2)){
			count++;
		}
		else if(points->X[i] < ((points->Y[i]-b3)/a3) && points->Y[i] < (a3*points->X[i]+b3)){
			count++;
		}
		else if(points->X[i] < ((points->Y[i]-b4)/a4) && points->Y[i] < (a4*points->X[i]+b4)){
			count++;
		}
		else if(points->X[i] > ((points->Y[i]-b5)/a5) && points->Y[i] < (a5*points->X[i]+b5)){
			count++;
		}
		else if(points->X[i] > ((points->Y[i]-b6)/a6) && points->Y[i] < (a6*points->X[i]+b6)){
			count++;
		}
		else if(points->X[i] > ((points->Y[i]-b7)/a7) && points->Y[i] > (a7*points->X[i]+b7)){
			count++;
		}
		else if(points->X[i] > ((points->Y[i]-b8)/a8) && points->Y[i] > (a8*points->X[i]+b8)){
			count++;
		}
	}
	//cout << count << endl;
	long resto = points->n - count;
	//cout << resto << endl;
	percent = (float)resto/(float)points->n;
	//cout << percent << endl;
	strcpy(aFile, "./RESULTS/");
	strcpy(str, "");
	sprintf(str, "EPSFilteredPerc"); //%ld, points->n
	strcat(aFile, str);
	cout << "Resume File: " << aFile << endl;
	FILE *fp = fopen(aFile, "a+" );
	if(NTHREADS > 0){
		if (points->NORMAL){
			// [n] [filtered percent] [esperanza] [varianza]
			fprintf(fp, "%ld %f %d %d\n", points->n, percent, points->NORMAL, NTHREADS);
		}else{
			// [n] [filtered percent] [avg bs-time/exec]
			fprintf(fp, "%ld %f %d %d\n", points->n, percent, points->NORMAL, NTHREADS);
		}
		fclose(fp);
	}
	else{
		if (points->NORMAL){
			// [n] [filtered percent] [esperanza] [varianza]
			fprintf(fp, "%ld %f %d %d\n", points->n, percent, points->NORMAL, 0);
		}else{
			// [n] [filtered percent] [avg bs-time/exec]
			fprintf(fp, "%ld %f %d %d\n", points->n, percent, points->NORMAL, 0);
		}
		fclose(fp);
	}
}

/**
 * Revisar la correctitud del algoritmo paralelo
 */
 void testPoints(pointSet *points){
	 float xups, xdowns, xlefts, xrights, xri_ups, xup_les, xle_dos, xdo_ris;
	 float yups, ydowns, ylefts, yrights, yri_ups, yup_les, yle_dos, ydo_ris;

	 runEPS(points);
	 xups=xup; yups=yup; xdowns=xdown; ydowns=ydown; xlefts=xleft; ylefts=yleft; xrights=xright; yrights=yright;
	 xri_ups=xri_up; xup_les=xup_le; xle_dos=xle_do; xdo_ris=xdo_ri;
	 yri_ups=yri_up; yup_les=yup_le; yle_dos=yle_do; ydo_ris=ydo_ri;
	 runEPSP(points);
	 if(xup != xups) cout << "ERROR, NORTH X coordinates aren't equal" << endl;
	 if(yup != yups) cout << "ERROR, NORTH Y coordinates aren't equal" << endl;
	 if(xdown != xdowns) cout << "ERROR, SOUTH X coordinates aren't equal" << endl;
	 if(ydown != ydowns) cout << "ERROR, SOUTH Y coordinates aren't equal" << endl;
	 if(xleft != xlefts) cout << "ERROR, WEST X coordinates aren't equal" << endl;
	 if(yleft != ylefts) cout << "ERROR, WEST Y coordinates aren't equal" << endl;
	 if(xright != xrights) cout << "ERROR, EAST X coordinates aren't equal" << endl;
	 if(yright != yrights) cout << "ERROR, EAST Y coordinates aren't equal" << endl;
	 if(xup_le != xup_les) cout << "ERROR, NORTHWEST X coordinates aren't equal" << endl;
	 if(yup_le != yup_les) cout << "ERROR, NORTHWEST Y coordinates aren't equal" << endl;
	 if(xdo_ri != xdo_ris) cout << "ERROR, SOUTHEAST X coordinates aren't equal" << endl;
	 if(ydo_ri != ydo_ris) cout << "ERROR, SOUTHEAST Y coordinates aren't equal" << endl;
	 if(xle_do != xle_dos) cout << "ERROR, SOUTHWEST X coordinates aren't equal" << endl;
	 if(yle_do != yle_dos) cout << "ERROR, SOUTHWEST Y coordinates aren't equal" << endl;
	 if(xri_up != xri_ups) cout << "ERROR, NORTHEAST X coordinates aren't equal" << endl;
	 if(yri_up != yri_ups) cout << "ERROR, NORTHEAST Y coordinates aren't equal" << endl;
 }
