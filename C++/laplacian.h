
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
struct Node **voisins;// un tableau de pointeurs avec les voisins
struct Edge **incidentes;    // tableau edges incidentes TODO
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
Tree theTree;
};

/* //////////////////////////// Les fonctions ///////////////////////////////////// */

int 			inTree(Problem *theProblem, int edgeIndex, int *treeIndex);
int 			find(int i);
int 			uni(int i,int j);
int 			findIndex(int a,int b,Problem *theProblem);
Problem 		*createProblem(char *FileName);
Chemin* 		findCycle(Edge *edgeCurrent,Problem *theProblem);
void 			Solve(char *FileName);
void 			edgeSort(Problem *theProblem);
void 			Kruskal(Problem *theProblem);
int 			edgeCompare( const void *edgea, const void *edgeb);
//Chemin* 		DFS(Edge *edgeCurrent,Problem *theProblem);
//int 			edgeCompare( const Edge *edgea, const Edge *edgeb);
#endif
