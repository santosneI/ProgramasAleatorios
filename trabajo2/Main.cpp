/*
 * Main.cpp
 *
 *  Created on: 13 sept. 2018
 *      Author: arias
 */
    
#include <CImg.h>
#include <math.h>
#include <stdio.h>
#include <immintrin.h>
using namespace cimg_library;
using namespace std;

/**********************************
 * TODO
 * 	- Change the data type returned by CImg.srcImage to adjust the
 * 	requirements of your workgroup
 * 	- Change the data type of the components pointers to adjust the
 * 	requirements of your workgroup
 */
//funcion retorna minimo color, pasandole como parametro el puntero que apunta al primer elemento de cada color
	double valorMin(double *pXcomp, int height, int width){
		double min = 255;
		__m256d min1;
		__m256d packet;
		for (int i = 0; i < height * width;i+=4) {
			packet = _mm256_loadu_pd (pXcomp);
			min1 = _mm256_min_pd (packet,_mm256_set1_pd (min));
			min = _mm256_cvtsd_f64 (min1);
			pXcomp+=4;
		}
		return min;
}
//funcion retorna maximo de cada color, pasandole como parametro el puntero que apunta al primer element de cada color
	double valorMax(double *pXcomp, int height, int width){
			double max = 0;
			__m256d max1;
			__m256d packet;
			for (int i = 0; i < height * width;i+=4) {
			packet = _mm256_loadu_pd (pXcomp);
			max1 = _mm256_max_pd (packet,_mm256_set1_pd (max));
			max = _mm256_cvtsd_f64 (max1);
			pXcomp+=4;
		}
			return max;
		}
int main() {
	CImg<double> srcImage("bailarina.bmp"); // Open file and object initialization

	double *pRcomp, *pGcomp, *pBcomp; // Pointers to the R, G and B components
	double *pRnew, *pGnew, *pBnew;
	double *pdstImage; // Pointer to the new image pixels
	int width, height; // Width and height of the image
	int nComp; // Number of image components
	struct timespec tStart, tEnd;
	double dElapsedTimeS;

	/***************************************************
	 *
	 * Variables initialization.
	 * Preparation of the necessary elements for the algorithm
	 * Out of the benchmark time
	 *
	 */

	srcImage.display(); // If needed, show the source image
	width = srcImage.width(); // Getting information from the source image
	height  = srcImage.height();
	nComp = srcImage.spectrum(); // source image number of components
				// Common values for spectrum (number of image components):
				//  B&W images = 1
				//	Normal color images = 3 (RGB)
				//  Special color images = 4 (RGB and alpha/transparency channel)


	// Allocate memory space for the pixels of the destination (processed) image 
	pdstImage = (double *) malloc (width * height * nComp * sizeof(double));
	if (pdstImage == NULL) {
		printf("\nMemory allocating error\n");
		exit(-2);
	}

	// Pointers to the RGB arrays of the source image
	pRcomp = srcImage.data(); // pRcomp points to the R component
	pGcomp = pRcomp + height * width; // pGcomp points to the G component
	pBcomp = pGcomp + height * width; // pBcomp points to B component
	// Pointers to the RGB arrays of the destination image
	pRnew = pdstImage;
	pGnew= pRnew + height * width;
	pBnew= pGnew + height * width;


	/*********************************************
	 * Algorithm start
	 *
	 * Measure initial time
	 *
	 *	COMPLETE
	 *
	 */
	
	if (clock_gettime(CLOCK_REALTIME, &tStart))
	{
		printf("ERROR: clock_gettime:\n");
		exit(EXIT_FAILURE);
	}
	
	 //numero de pixeles para el bucle y valores maximos y minimos de los colores de la imagen
        int numPixels = height * width;
        double minRed = valorMin(pRcomp, height, width);
        double minGreen = valorMin(pGcomp, height, width);
        double minBlue = valorMin(pBcomp, height, width);

        double maxRed = valorMax(pRcomp, height, width);
        double maxGreen = valorMax(pGcomp, height, width);
        double maxBlue = valorMax(pBcomp, height, width);


	/************************************************
	 * Algorithm.
	 * In this example, the algorithm is a components exchange
	 *
	 * TO BE REPLACED BY YOUR ALGORITHM
	 */

	for(int i = 0; i < 11 ; i++){//bucle para que el algoritmo dure mas de 5 segundos
		
		//paquetes double para componentes de la imagen fuente
		__m256d rPacket, gPacket, bPacket;
		//paquetes double para componentes de la imagen destino
		__m256d rNewPacket, gNewPacket, bNewPacket; 
		//paquetes para realizar division
		__m256d rSub, gSub, bSub;
		//variables de tipo double donde se restaran max y min de cada componente
		double rMsub, gMsub, bMsub;
		rMsub = maxRed - minRed;
		gMsub = maxGreen - minGreen;
		bMsub = maxBlue - minBlue; 
		for(int i = 0;i < numPixels ; i += 4){
			//Cargamos valores de las componentes en los paquetes
			rPacket = _mm256_loadu_pd (pRcomp);
			gPacket = _mm256_loadu_pd(pGcomp);
			bPacket = _mm256_loadu_pd(pBcomp);
			
			//rSub = (rPacket - minRed)
			rSub = _mm256_sub_pd (rPacket, _mm256_set1_pd (minRed));

			//hacemos lo mismo con demas componentes
			gSub = _mm256_sub_pd (gPacket, _mm256_set1_pd (minGreen));
			bSub = _mm256_sub_pd (bPacket, _mm256_set1_pd (minBlue));
			
			//rNewPacket = ((rPacket - minRed) / (maxRed - minRed)) * 255; == rNewPacket = (rSub / rMsub) * 255
			rNewPacket = _mm256_mul_pd (_mm256_div_pd (rSub,_mm256_set1_pd (rMsub)),_mm256_set1_pd (255));

			//hacemos lo mismo con las demas componentes
			gNewPacket = _mm256_mul_pd (_mm256_div_pd (gSub,_mm256_set1_pd (gMsub)),_mm256_set1_pd (255));
			bNewPacket = _mm256_mul_pd (_mm256_div_pd (bSub,_mm256_set1_pd (bMsub)),_mm256_set1_pd (255));
//			Aplicamos saturacion si es necesario
//			if(rNewPacket > 255){
//				rNewPacket = 255;
//			}
			//minimo entre valores y 255
			rNewPacket = _mm256_min_pd(rNewPacket, _mm256_set1_pd(255));
			//hacemos lo mismo con las demas componentes
			gNewPacket = _mm256_min_pd(gNewPacket, _mm256_set1_pd(255));
			bNewPacket = _mm256_min_pd(bNewPacket, _mm256_set1_pd(255));

			// Pointers to access the RGB packets as integers
			double *pR = (double *)&rNewPacket;
			double *pG =  (double *)&gNewPacket;
			double *pB = (double *)&bNewPacket;

			// Save the results in the destination arrays
			for (int k = 0; k < 4; k++)
			{
				*pRnew = (double)*pR;
				*pGnew = (double)*pG;
				*pBnew = (double)*pB;
				pRnew++; pGnew++; pBnew++;
				pR++; pG++; pB++;
			}
			// Prepare for the next 4 double
			pRcomp+=4; pGcomp+=4; pBcomp+=4;
		}
			//reiniciamos todos los punteros a su valor inicial
			// Pointers to the RGB arrays of the destination image
			pRcomp -= numPixels; pBcomp -= numPixels; pGcomp -= numPixels;
			pRnew -= numPixels; pBnew -= numPixels; pGnew -= numPixels;

	}
	/***********************************************
	 * End of the algorithm
	 *
	 * Measure the final time and calculate the time spent
	 *
	 * COMPLETE
	 *
	 */
	 if (clock_gettime(CLOCK_REALTIME, &tEnd))
	{
		printf("ERROR: clock_gettime:\n");
		exit(EXIT_FAILURE);
	}
	 dElapsedTimeS = (tEnd.tv_sec - tStart.tv_sec);
	 dElapsedTimeS += (tEnd.tv_nsec - tStart.tv_nsec) / 1e+9;
	 printf("Tiempo de Ejecucion del algoritmo: "); 
	 printf("%f\n",dElapsedTimeS);

	// Create a new image object with the calculated pixels
	// In case of normal color image use nComp=3,
	// In case of B&W image use nComp=1.
	CImg<double> dstImage(pdstImage, width, height, 1, nComp);
	// Store the destination image in disk
	dstImage.save("bailarina2-256-double.bmp");
	// Display the destination image
	dstImage.display(); // If needed, show the result image
	return(0);

}

