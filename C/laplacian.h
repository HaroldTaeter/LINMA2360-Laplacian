
/*
 *  laplacian.h
 */

#ifndef _LAPLACIAN_H_
#define _LAPLACIAN_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* //////////////////////////// Les STRUCTURES ///////////////////////////////////// */
typedef struct Node Node;
struct Node {
int degree;// = nChild 
int indice;
struct Node **voisins;// un tableau de pointeurs avec les voisins
struct Edge **incidentes;    // tableau edges incidentes TODO
}; 

typedef struct Tree Tree;
struct Tree {
int nChild;
Node *node;
struct Tree **Child;// un tableau de pointeurs avec les voisins
};

typedef struct Edge Edge;
struct Edge {
int indice; 
double weight;
double f; // le flot ? on pourra surement se passer de stocker cela là. Voir vecteur de section 5 
Node *a;
Node *b;
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
void 		addToTree(Problem *theProblem, Tree *arbre, int *treeIndex);
int 		inTree(Problem *theProblem, int edgeIndex, int *treeIndex);
int 		find(int i);
int 		uni(int i,int j);
int 		findIndex(int a,int b,Problem *theProblem);
Problem 	*createProblem(char *FileName);
void 		Solve(char *FileName);
void 		edgeSort(Problem *theProblem);
int* 		Kruskal(Problem *theProblem);
int 		edgeCompare( const void *edgea, const void *edgeb);
//int edgeCompare( const Edge *edgea, const Edge *edgeb);
#endif
