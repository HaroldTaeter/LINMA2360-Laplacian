
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
};

typedef struct Tree Tree;
struct Tree {
Node *nodeSource;
int *predecessor;
Edge **edgesTree;
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
double **Weights;
double *b;
Tree theTree;
Edge **edgesOffTree;
double *theCumulatedProba;
};

/* //////////////////////////// Les fonctions ///////////////////////////////////// */

int 			inTree(Problem *theProblem, int edgeIndex, int *treeIndex);
int 			find(int i);
int 			uni(int i,int j);
int 			findIndex(int a,int b,Problem *theProblem);
Problem 		*createProblem(char *FileNameL, char *FileNameb);
Chemin* 		findCycle(Edge *edgeCurrent,Problem *theProblem);
void 			Solve(char *FileName);
void 			edgeSort(Problem *theProblem);
void 			Kruskal(Problem *theProblem);
int 			edgeCompare( const void *edgea, const void *edgeb);
void			setFlow(int indexNodeA, int indexNodeB, Problem *theProblem);
Chemin* 		findPath(int IndexNodeA, int IndexNodeB, Problem *theProblem);
//Chemin* 		DFS(Edge *edgeCurrent,Problem *theProblem);
double          stretchEdge(Edge *edgeCurrent, Chemin *Chemin);
double          stretchTree(Problem *theProblem);
double          probabilityEdge(Edge *edgeCurrent, Problem *theProblem);
double*         probaCompute(Problem *theProblem);
int             iterationsK(Problem *theProblem, double eps);
void            CycleUpdate(Edge *edgeCurrent, Problem *theProblem);
Edge*           RandomPicking(Problem *theProblem, Edge **edgesOffTree);
double*         InducedVoltages(Problem *theProblem);
double*         InducedVoltages2(Problem *theProblem);
#endif
