//#include <stdio.h>
//#include <stdlib.h>
#include <time.h>
//#include <math.h>
#include <fstream>
#include <string>
#include <iostream>
#include <queue>
#include <new>

//#include "laplacian.cpp"

using namespace std;
void Solve(char *FileNameL,char *FileNameb);

int main(void)
{
       char *FileNameL="dataL.txt";
       char *FileNameb="datab.txt";
       
       struct timespec start, inter1, finish;
       double elapsed;
	
       clock_gettime(CLOCK_MONOTONIC, &start);
       Solve(FileNameL,FileNameb);
       clock_gettime(CLOCK_MONOTONIC, &finish);
 
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	printf("time_elapsed: %f\n",elapsed);
       
       //exit(0);
       printf("Execution OVER \n");
}

