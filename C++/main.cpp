//#include <stdio.h>
//#include <stdlib.h>
//#include <time.h>
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
       Solve(FileNameL,FileNameb);
       //exit(0);
       printf("Execution OVER \n");
}

