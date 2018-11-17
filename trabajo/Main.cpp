/*
 * Main.cpp
 *
 *  Created on: 13 sept. 2018
 *      Author: arias
 */
    
#include <CImg.h>
#include <math.h>
#include <stdio.h>
using namespace cimg_library;

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
		for (int i = 0; i < height * width;i++) {
			if(pXcomp[i] < min){
				min = pXcomp[i];
			}
		}

		return min;
	}
//funcion retorna maximo de cada color, pasandole como parametro el puntero que apunta al primer element de cada color
	double valorMax(double *pXcomp, int height, int width){
			double max = 0;
			for (int i = 0; i < height * width ;i++) {
				if(pXcomp[i] > max){
					max = pXcomp[i];
				}
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
	clock_gettime(CLOCK_REALTIME, &tStart);


	/************************************************
	 * Algorithm.
	 * In this example, the algorithm is a components exchange
	 *
	 * TO BE REPLACED BY YOUR ALGORITHM
	 */
/*
	for (int i = 0; i < height ; i++) {
		for (int j=0; j < width; j++)  {
			pRnew[i * width + j] = ((pRcomp[i * width +j] - ((*pRcomp - valorMin(pRcomp, height, width)) / (valorMax(pRcomp, height, width) - valorMin(pRcomp,height,width))) * 255;
			pGnew[i * width + j] = ((pGcomp[i * width +j] - ((*pGcomp - valorMin(pGcomp, height, width)) / (valorMax(pGcomp, height, width) - valorMin(pGcomp,height,width))) * 255;
			pBnew[i * width + j] = (((pBcomp[i * width +j]) - ((*pBcomp - valorMin(pBcomp, height, width)) / (valorMax(pBcomp, height, width) - valorMin(pBcomp,height,width))) * 255;
		}
	}*/

	int numPixels = height * width;
	double minRed = valorMin(pRcomp, height, width);
	double minGreen = valorMin(pGcomp, height, width);
	double minBlue = valorMin(pBcomp, height, width);

	double maxRed = valorMax(pRcomp, height, width);
	double maxGreen = valorMax(pGcomp, height, width);
	double maxBlue = valorMax(pBcomp, height, width);

	for(int i = 0;i < numPixels ; i++){
		*pRnew = ((*pRcomp - minRed) / (maxRed - minRed)) * 255;
		*pGnew = ((*pGcomp - minGreen) / (maxGreen - minGreen)) * 255;
		*pBnew = ((*pBcomp - minBlue) / (maxBlue - minBlue)) * 255;

		if(*pRnew > 255){
			*pRnew = 255;
		}

		if(*pGnew > 255){
			*pGnew = 255;
		}

		if(*pBnew > 255){
			*pBnew = 255;
		}

		pRnew++;
		pGnew++;
		pBnew++;
		pRcomp++;
		pGcomp++;
		pBcomp++;


	}

	/***********************************************
	 * End of the algorithm
	 *
	 * Measure the final time and calculate the time spent
	 *
	 * COMPLETE
	 *
	 */
	 clock_gettime(CLOCK_REALTIME, &tEnd);
	 dElapsedTimeS = (tEnd.tv_sec - tStart.tv_sec);
	 dElapsedTimeS += (tEnd.tv_nsec - tStart.tv_nsec) / 1e+9;
	 printf("Tiempo de Ejecucion del algoritmo: "); 
	 printf("%f\n",dElapsedTimeS);
		
	// Create a new image object with the calculated pixels
	// In case of normal color image use nComp=3,
	// In case of B&W image use nComp=1.
	CImg<double> dstImage(pdstImage, width, height, 1, nComp);
	// Store the destination image in disk
	dstImage.save("bailarina2.bmp"); 
	// Display the destination image
	dstImage.display(); // If needed, show the result image
	return(0);

}

