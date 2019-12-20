/*
* Name: Archit Jaiswal
* UTA ID: 1001543326
* This program creates a Mandelbrot set image at the input coordinates and it also computes using number of threads entered by the user
*/

#include "bitmap.h"
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <string.h>
#include <pthread.h>

int iteration_to_color( int i, int max );
int iterations_at_point( double x, double y, int max );
void compute_image(void* dataArr);

// create a struct to store all the data and pass in compute image through "pthread_create" call
typedef struct {
	struct bitmap *bm;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	int max;
	int imageStartHeight;
	int imageEndHeight;
} ImageData;

void show_help() // prints help to provvide details of input arguements
{
	printf("Use: mandel [options]\n");
	printf("Where options are:\n");
	printf("-x <coord>   X coordinate of image center point. (default=0)\n");
	printf("-y <coord>   Y coordinate of image center point. (default=0)\n");
	printf("-s <scale>   Scale of the image in Mandlebrot coordinates. (default=4)\n");
  printf("-m <max>     The maximum number of iterations per point. (default=1000)\n");
	printf("-W <pixels>  Width of the image in pixels. (default=500)\n");
	printf("-H <pixels>  Height of the image in pixels. (default=500)\n");
	printf("-n <threads> Number of threads to compute the image. (default=1)\n");
	printf("-o <file>    Set output file. (default=mandel.bmp)\n");
	printf("-h           Show this help text.\n");
	printf("\nSome examples are:\n");
	printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
	printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
	printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000\n\n");
}


int main( int argc, char *argv[] )
{

	char c;
    int i;

	// These are the default configuration values used
	// if no command line arguments are given.
	const char *outfile = "mandel.bmp";
	double xcenter = 0;
	double ycenter = 0;
	double scale = 4;
	int    image_width = 500;
	int    image_height = 500;
	int    max = 1000;
	int numberOfThreads = 1;


	// For each command line argument given,
	// override the appropriate configuration value.

	while((c = getopt(argc,argv,"x:y:s:W:H:m:n:o:h"))!=-1) {
		switch(c) {
			case 'x':
				xcenter = atof(optarg);
				break;
			case 'y':
				ycenter = atof(optarg);
				break;
			case 's':
				scale = atof(optarg);
				break;
			case 'W':
				image_width = atoi(optarg);
				break;
			case 'H':
				image_height = atoi(optarg);
				break;
			case 'm':
				max = atoi(optarg);
				break;
			case 'n':
			  numberOfThreads = atoi(optarg);
				break;
			case 'o':
				outfile = optarg;
				break;
			case 'h':
				show_help();
				exit(1);
				break;
		}
	}

	// Display the configuration of the image.
	printf("mandel: x=%lf y=%lf scale=%lf max=%d width=%d height=%d threads=%d outfile=%s\n"
							 ,xcenter,ycenter,scale,max,image_width,image_height,numberOfThreads,outfile);

    //printf("%d is the numberOfThreads\n", numberOfThreads);
    ImageData dataArr[numberOfThreads]; //creates a struct of all the configuration described by the user and keeps it in the struct
	// "dataArr" struct will be passed to "compute_image" via "pthread_create"

	pthread_t imageProcessThread[numberOfThreads]; // creates the number of thread required by the user


	// Create a bitmap of the appropriate size.
	struct bitmap *bm = bitmap_create(image_width,image_height);

  // stores the input parameters of all different inputs in "compute_image"
	for( i = 0; i < numberOfThreads; i++)
	{
		//printf("First:%d     Second:%d\n",i*(image_height/numberOfThreads), (i+1)*(image_height/numberOfThreads) );
		dataArr[i].bm = bm;
		dataArr[i].xmin = xcenter-scale;
		dataArr[i].xmax = xcenter+scale;
		dataArr[i].ymin = ycenter-scale;
		dataArr[i].ymax = ycenter+scale;
		dataArr[i].max = max;
		dataArr[i].imageStartHeight = (int) i*(image_height/numberOfThreads);

		// if the image is not evenly distributable among all the threads then give all the remaining part to last thread
		if( i == (numberOfThreads - 1 ) )
		{
			//printf("The actual width of the image is: %d\n", bitmap_width(bm));
		  dataArr[i].imageEndHeight = image_height;
		}
		else{
      //printf("It comes in else at i = %d\n", i);
			dataArr[i].imageEndHeight =(int) ((i+1)*(image_height/numberOfThreads));
		}

		// Fill it with a dark blue, for debugging
		//bitmap_reset(bm,MAKE_RGBA(0,0,255,0));

		if( pthread_create(&imageProcessThread[i], NULL, compute_image,(void*) &dataArr[i]) ) // Creates the thread and passes it to compute_image with the respective input parameters
		{
      perror("Error creating thread: ");
      exit( EXIT_FAILURE );
    }
	}

  // It waits till all the threads are done processing
  for (i =0; i < numberOfThreads; i++)
	{
		//printf("Waiting for thread#: %d\n",i );
		if(pthread_join(imageProcessThread[i], NULL))
		{
			perror("Problem with pthread_join: ");
		}
	}

	// Save the image in the stated file.
	if(!bitmap_save(bm,outfile)) {
		fprintf(stderr,"mandel: couldn't write to %s: %s\n",outfile,strerror(errno));
		return 1;
	}


	return 0;
}

/*
Compute an entire Mandelbrot image, writing each point to the given bitmap.
Scale the image to the range (xmin-xmax,ymin-ymax), limiting iterations to "max"
*/

void compute_image(void* dataArr)
{
	ImageData* arr =  (ImageData*) dataArr; // casting the void pointer back to original type
	int i,j;

	int width = bitmap_width(arr->bm); //this is same as the width of full
	int height = bitmap_height(arr->bm); // this is same as height of full image

	//printf("imageStartHeight: %d    endHeight: %d\n", arr->imageStartHeight, arr->imageEndHeight );

	// For every pixel in the image...

	for( j=arr->imageStartHeight; j <arr->imageEndHeight; j++) { // this keeps the start and stop controll of an image

		for(i=0;i<width;i++) {

			// Determine the point in x,y space for that pixel.
			double x = arr->xmin + i*(arr->xmax - arr->xmin)/width;
			double y = arr->ymin + j*(arr->ymax - arr->ymin)/height;

			// Compute the iterations at that point.
			int iters = iterations_at_point(x,y,arr->max);

			// Set the pixel in the bitmap.
			bitmap_set(arr->bm,i,j,iters);
		}
	}
}

/*
Return the number of iterations at point x, y
in the Mandelbrot space, up to a maximum of max.
*/

int iterations_at_point( double x, double y, int max )
{
	double x0 = x;
	double y0 = y;

	int iter = 0;

	while( (x*x + y*y <= 4) && iter < max ) {

		double xt = x*x - y*y + x0;
		double yt = 2*x*y + y0;

		x = xt;
		y = yt;

		iter++;
	}

	return iteration_to_color(iter,max);
}

/*
Convert a iteration number to an RGBA color.
Here, we just scale to gray with a maximum of imax.
Modify this function to make more interesting colors.
*/

int iteration_to_color( int i, int max )
{
	int gray = 255*i/max;
	return MAKE_RGBA(gray,gray,gray,0);
}
