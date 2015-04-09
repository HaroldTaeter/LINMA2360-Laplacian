//#include <stdio.h>
//#include <stdlib.h>
#include <time.h>
//#include <math.h>
#include <fstream>
#include <string>
#include <iostream>
#include <queue>
#include <new>


//#define _edges300_
//#define _edges1000_
//#define _edges3000_
#define _edges9000_

//#include "laplacian.cpp"

using namespace std;
void Solve(char *FileNameL,char *FileNameb);

int main(void)
{
       #ifdef _edges300_
       char *FileNameL="edges300.txt";
       char *FileNameb="datab300.txt";
       #endif
  	#ifdef _edges1000_
       char *FileNameL="edges1000.txt";
       char *FileNameb="datab1000.txt";
       #endif

       #ifdef _edges3000_
       char *FileNameL="edges3000.txt";
       char *FileNameb="datab3000.txt";
       #endif

	#ifdef _edges9000_
       char *FileNameL="edges9000.txt";
       char *FileNameb="datab9000.txt";
       #endif
       struct timespec start, inter1, finish;
       double elapsed;

       clock_gettime(CLOCK_MONOTONIC, &start);
       Solve(FileNameL,FileNameb);
       clock_gettime(CLOCK_MONOTONIC, &finish);

	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	printf("time_elapsed for total execution: %f\n",elapsed);

       //exit(0);
       printf("Execution OVER \n");
}

