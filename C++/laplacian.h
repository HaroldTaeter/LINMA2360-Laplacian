
/*
 *  laplacian.h
 */

#ifndef _LAPLACIAN_H_
#define _LAPLACIAN_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <queue>
#include <iostream>

/* //////////////////////////// Les STRUCTURES ///////////////////////////////////// */

typedef struct Node Node;
struct Node {
int degree;
int indice;
struct Node **voisins;// Assez redondant ..
struct Edge **incidentes;
};

typedef struct Edge Edge;
struct Edge {
int indice;
double weight;
double f;
Node *a;
Node *b;
double stretch;
struct Chemin *edgeChemin;
};

typedef struct Tree Tree;
struct Tree {
Node *nodeSource;
int *predecessor;
Edge **edgesTree;
double stretch;
};

typedef struct Chemin Chemin;
struct Chemin{
int size;
Edge **theChemin;
};

typedef struct Problem Problem;
struct Problem {
int nNode;
int nEdge;
Node *nodes;
Edge *edges;
Edge *edgesSorted;
Edge **edgesOffTree;
double **Weights;
double *b;
Tree theTree;
double *theCumulatedProba;
int *parent;
};

/* //////////////////////////// Les fonctions ///////////////////////////////////// */
double 			testKPL(Edge *edgeCurrent);
int 			inTree(Edge *edgea,Problem *theProblem);
int 			find(int i, Problem *theProblem);
int 			uni(int i,int j, Problem *theProblem);
int 			findIndex(int a,int b,Problem *theProblem);
Problem 		*createProblem(char *FileNameL, char *FileNameb);
Chemin* 		findCycle(Edge *edgeCurrent,Problem *theProblem);
void 			Solve(char *FileName);
void 			edgeSort(Problem *theProblem);
void 			Kruskal(Problem *theProblem);
int 			edgeCompare( const void *edgea, const void *edgeb);
void			setFlow(int indexNodeA, int indexNodeB, Problem *theProblem);
Chemin* 		findPath(int IndexNodeA, int IndexNodeB, Problem *theProblem);
Edge* 			RandomPicking(Problem *theProblem, Edge **edgesOffTree);
int 			iterationsK(Problem *theProblem, double eps);
double* 		probaCompute(Problem *theProblem);
double        		stretchEdge(Edge *edgeCurrent, Chemin *Chemin);
double          	probabilityEdge(Edge *edgeCurrent, Problem *theProblem);
double*         	probaCompute(Problem *theProblem);
int             	iterationsK(Problem *theProblem, double eps);
void            	CycleUpdate(Edge *edgeCurrent, Problem *theProblem);
Edge*           	RandomPicking(Problem *theProblem, Edge **edgesOffTree);
double*         	InducedVoltages2(Problem *theProblem);
void                	stretchsAndChemins(Problem *theProblem);
void 			setFlow2(Problem *theProblem, double *ksi);
void 			SetPredecessor(Problem *theProblem);
void 			SetNodes(Problem *theProblem);
void 			FreeLaplacian(Problem *theProblem);
#endif
