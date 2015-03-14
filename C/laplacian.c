/*
 *  laplacian.c
 *  
 */

#include "laplacian.h"

static int parent[9];

void Solve(char *FileName)
{
	Problem *theProblem= createProblem(FileName);
	
	/*printf(" node 1: indice=%d \n ",theProblem->nodes[1].indice);
	printf(" node 1: voisins[0].indice=%d \n ",theProblem->nodes[1].voisins[1]->indice);
	printf("edge1.b: indice= %d \n ",theProblem->edges[1].b->indice);*/
	
	int *treeIndex=Kruskal(theProblem);
	addToTree(theProblem, &theProblem->theTree, treeIndex);
	
	// node 0
	printf("theProblem->theTree.node->indice: %d \n  ", theProblem->theTree.node->indice);
	printf("theProblem->theTree.node->Child[0]->nChild: %d \n  ", theProblem->theTree.Child[0]->nChild);	


}





void exploreDFS(Tree *tree, int nodeIndex)
{
	if( tree->node->indice == nodeIndex )
	{// on est arrivé !
		
	}
	else if (tree->nChild > 0)
	{
		explore(tree->Child[0], nodeIndex);
	}
	// else : on est en bas, on fait rien

}



Problem *createProblem(char *FileName)
{
//int j=0;
//size_t count;
//char *line = malloc(100);

Problem *theProblem = malloc(sizeof(Problem));
int i,trash;    
    FILE* file = fopen(FileName,"r");

    trash = fscanf(file, "Number of nodes %d \n", &theProblem->nNode);
    int nNode=theProblem->nNode;
    printf("nNode: %d \n", nNode);
    theProblem->nodes=malloc(nNode*sizeof(Node));
	    
    double** theArray;  
    theArray = (double**) malloc(nNode*sizeof(double*));  
    for (i = 0; i < nNode; i++)  
       theArray[i] = (double*) malloc(nNode*sizeof(double)); 
    	    
    theProblem->Weights=theArray;
    printf("on est là\n");
    
  /*  for(i = 0; i < nNode; i++)
  {
      for(j = 0; j < nNode; j++) 
      {
  //Use lf format specifier, %c is for character
       if (!fscanf(file, "%lf", &theProblem->Weights[i][j])) 
       {
       	printf("break\n");
        break;
       }
      // mat[i][j] -= '0'; 
       printf("%lf \n",theProblem->Weights[i][j]); //Use lf format specifier, \n is for new line
      }
  }*/
    
    /* while(getline(&line, &count, file)!=-1) {
        for (; count > 0; count--, j++)
            sscanf(line, "%lf", &theProblem->Weights[i][j]);
        i++;
    }*/
    fclose(file);
   /* for (i = 0; i < nNode; i++)
    {
    	for(j=0; j < nNode; j++)
    	{
    		fscanf(" .. ", theProblem->Weights[i][j]);    	
    	}
    	//ici il faut aller à la ligne
    }*/
	
//int nNode=4;
int nEdge=4;
theProblem->nNode=nNode;
theProblem->nEdge=nEdge;
theProblem->nodes=malloc(nNode*sizeof(Node));
theProblem->edges=malloc(4*sizeof(Edge));

//NODE0
theProblem->nodes[0].indice=0;
theProblem->nodes[0].degree=1;
theProblem->nodes[0].voisins=malloc(1*sizeof(Node));
theProblem->nodes[0].voisins[0]=&theProblem->nodes[1];
//NODE1
theProblem->nodes[1].indice=1;
theProblem->nodes[1].degree=3;
theProblem->nodes[1].voisins=malloc(3*sizeof(Node));
theProblem->nodes[1].voisins[0]=&theProblem->nodes[0];
theProblem->nodes[1].voisins[1]=&theProblem->nodes[2];
theProblem->nodes[1].voisins[2]=&theProblem->nodes[3];
//NODE2
theProblem->nodes[2].indice=2;
theProblem->nodes[2].degree=2;
theProblem->nodes[2].voisins=malloc(2*sizeof(Node));// degree*sizeof(Node)
theProblem->nodes[2].voisins[0]=&theProblem->nodes[1];
theProblem->nodes[2].voisins[1]=&theProblem->nodes[3];
//NODE3
theProblem->nodes[3].indice=3;
theProblem->nodes[3].degree=2;
theProblem->nodes[3].voisins=malloc(2*sizeof(Node));// degree*sizeof(Node)
theProblem->nodes[3].voisins[0]=&theProblem->nodes[1];
theProblem->nodes[3].voisins[1]=&theProblem->nodes[2];

//EDGE0
theProblem->edges[0].indice=0;
theProblem->edges[0].weight=12.0;
theProblem->edges[0].f=0.0;
theProblem->edges[0].a=&theProblem->nodes[0];
theProblem->edges[0].b=&theProblem->nodes[1];
//EDGE1
theProblem->edges[1].indice=1;
theProblem->edges[1].weight=1.0;
theProblem->edges[1].f=0.0;
theProblem->edges[1].a=&theProblem->nodes[1];
theProblem->edges[1].b=&theProblem->nodes[2];
//EDGE2
theProblem->edges[2].indice=2;
theProblem->edges[2].weight=7.0;
theProblem->edges[2].f=0.0;
theProblem->edges[2].a=&theProblem->nodes[1];
theProblem->edges[2].b=&theProblem->nodes[3];
//EDGE3
theProblem->edges[3].indice=3;
theProblem->edges[3].weight=20.0;
theProblem->edges[3].f=0.0;
theProblem->edges[3].a=&theProblem->nodes[2];
theProblem->edges[3].b=&theProblem->nodes[3];

//////////// Incidentes Edges /////////////////

theProblem->nodes[0].incidentes=malloc(1*sizeof(Edge));
theProblem->nodes[1].incidentes=malloc(3*sizeof(Edge));
theProblem->nodes[2].incidentes=malloc(2*sizeof(Edge));
theProblem->nodes[3].incidentes=malloc(2*sizeof(Edge));

theProblem->nodes[0].incidentes[0]=&theProblem->edges[0];

theProblem->nodes[1].incidentes[0]=&theProblem->edges[0];
theProblem->nodes[1].incidentes[1]=&theProblem->edges[1];
theProblem->nodes[1].incidentes[2]=&theProblem->edges[2];

theProblem->nodes[2].incidentes[0]=&theProblem->edges[1];
theProblem->nodes[2].incidentes[1]=&theProblem->edges[3];

theProblem->nodes[3].incidentes[0]=&theProblem->edges[2];
theProblem->nodes[3].incidentes[1]=&theProblem->edges[3];


//////////// Adjacency Matrix //////////////////

theProblem->Weights[0][0]=0.0;
theProblem->Weights[1][1]=0.0;
theProblem->Weights[2][2]=0.0;
theProblem->Weights[3][3]=0.0;

theProblem->Weights[0][1]=12.0;
theProblem->Weights[1][0]=12.0;
theProblem->Weights[0][2]=0.0;
theProblem->Weights[2][0]=0.0;
theProblem->Weights[0][3]=0.0;
theProblem->Weights[3][0]=0.0;

theProblem->Weights[1][2]=1.0;
theProblem->Weights[2][1]=1.0;
theProblem->Weights[1][3]=7.0;
theProblem->Weights[3][1]=7.0;

theProblem->Weights[2][3]=20.0;
theProblem->Weights[3][2]=20.0;

//printf("Weights[1][2]= %f \n", theProblem->Weights[0][1]);

Tree arbreNew;
arbreNew.nChild=0;
arbreNew.node=&theProblem->nodes[0];
arbreNew.Child=NULL;
theProblem->theTree= arbreNew;
/*theProblem->theTree.nChild=0;
theProblem->theTree.node=&theProblem->nodes[0];
theProblem->theTree->Child=NULL;*/ // a virer si cela fonctionne

return theProblem;

}


int edgeCompare( const void * edgea, const void * edgeb)
{
	return (  (*(Edge*)edgea).weight- (*(Edge*)edgeb).weight);
}

void edgeSort(Problem *theProblem)
{
	qsort(theProblem->edges, theProblem->nEdge, sizeof(Edge), edgeCompare);
}

int* Kruskal(Problem *theProblem)
{
	/*edgeSort(theProblem);
	int i; 
	for(i=0; i<theProblem->nEdge;i++)
	{// on remet les indices comme il faut après avoir trié
		theProblem->edges[i].indice=i;
	}
	//int *indexEdge=malloc(sizeof(int)*(theProblem->nNode-1));
	*/
	
    int i,j,k,a,b,u,v,ne=1;
    int n = theProblem->nNode;
    int min,mincost=0;
    int Wmax=999;
    int *treeIndex= malloc((theProblem->nNode-1)*sizeof(int));
    int compteur=0;
    
    while(ne <theProblem->nNode)
    	 {
    		for(i=0,min=Wmax;i<n;i++)
    		{
    			for(j=0;j < n;j++)
    			{
    				if(theProblem->Weights[i][j] < min && theProblem->Weights[i][j] >0 )
    				{
    					min=theProblem->Weights[i][j];
    					a=u=i;
    					b=v=j;
    				}
    			}
    		}
    	u=find(u);
    	v=find(v);
    	if(uni(u,v))
    	{
    		
    		treeIndex[compteur]=findIndex(a,b,theProblem);
    		ne++;
    		//printf("%d edge (%d,%d) =%d\n",ne++,a,b,min);
    		compteur++;
    		mincost +=min;
    	}
    		theProblem->Weights[a][b]=theProblem->Weights[b][a]=Wmax;
    	}
    	
    	/*printf("treeIndex[0]=%d \n ", treeIndex[0]);
    	printf("treeIndex[1]=%d \n ", treeIndex[1]);
    	printf("treeIndex[2]=%d \n ", treeIndex[2]);*/
    	return treeIndex;
}


int find(int i)
{
	while(parent[i])
	i=parent[i];
	return i;
}
int uni(int i,int j)
{
    	if(i!=j)
    	{
    		parent[j]=i;
    		return 1;
    	}
    	return 0;
}
    
int findIndex(int a, int b,Problem *theProblem)
{
  	int i; 
    	for(i=0; i<theProblem->nEdge; i++)
    	{
    		if( ( theProblem->edges[i].a->indice==a || theProblem->edges[i].b->indice==a) && 
    			( theProblem->edges[i].a->indice==b || theProblem->edges[i].b->indice==b) )	
    		{	
    			return theProblem->edges[i].indice;
    		}
    	}
}
/////////////////////////////////////////////////////////////////////////////////
void 	addToTree(Problem *theProblem, Tree *arbre, int *treeIndex)
{
	Node *nodeLocal=arbre->node;	
	int i;
	arbre->Child=malloc(5*sizeof(Tree));// TODO take care de ce degree max de 5 dans l'arbre
	for(i=0;i<nodeLocal->degree; i++)
	{
		if (inTree(theProblem, nodeLocal->incidentes[i]->indice,treeIndex)==1)
		{
			Tree arbre1;
			arbre1.nChild=0;
			if( nodeLocal->incidentes[i]->a->indice == nodeLocal->indice )
			{
				arbre1.node=nodeLocal->incidentes[i]->b;
			}
			else if( nodeLocal->incidentes[i]->b->indice == nodeLocal->indice )
			{
				arbre1.node=nodeLocal->incidentes[i]->a;
			}
			
			else{ printf("bug: addToTree \n"); }
			arbre1.Child=NULL;
			arbre->Child[arbre->nChild]=&arbre1;
			arbre->nChild++;			
		}	
	
	}
	for(i=0; i<arbre->nChild; i++)
	{
		addToTree(theProblem, arbre->Child[i],treeIndex);	
	}

}

int inTree(Problem *theProblem, int edgeIndex, int *treeIndex)
{// doit dire si l'edge est dans l'arbre (return 0 ou 1) et mettre à jour somehow le tableau edgeIndex 
//  pour qu'on ajoute jamais 2 fois une edge :)
	int i; 
	
	for(i=0; i< theProblem->nNode-1; i++ )
	{
		if ( edgeIndex == treeIndex[i] )
		{
			treeIndex[i]=-1;// TODO changer, c'est super lent car on va check plein de fois des trucs inutiles
			return 1;
		}
	}
	return 0;
}


////////////////////////// HAROLD'S PRATICE /////////////////////////////

int stretchCompute(Edge *edgeCurrent, Edge **chemin)
{
	/*
		edgeCurrent est un pointeur vers une edge dont on veut calculer le stretch.
		chemin est un tableau dont les éléments sont des pointeurs de type Edge. 
		logiquement, le chemin appartient au tree et relie les deux extrémités de edgeCurrent.
		La fonction devrait faire genre 10 lignes max j'ai l'impression :-) 	
		
		Je ne sais pas si c'est eaxctement cette fonction là qu'on utilisera mais surement une version similaire
		(peut être arguents un peu différent selon comme cela se passe aux autres etapes). 
		Ca devrait te permettre de voir un peu les structures qu'on a déjà et chipoter pas mal avec des pointeurs :p
		Après tu peux regarder ce qu'il faut pour la suite et faire une fonction qui calcule les autres quantités necessaires genre
		la résistance d'un chemin ou des trucs comme ça (voir article :-). Tiens moi au courant souvent et hésite pas 
		à me demander de te débloquer souvent, tu devrais attraper le truc pluc vite en faisant comme ça ;-)
	*/

	/*
	
	
	
	
	*/


}







