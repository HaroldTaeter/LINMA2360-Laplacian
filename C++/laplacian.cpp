/*
 *  laplacian.c
 *
 */

#include "laplacian.h"

using namespace std;

static int parent[9];


/*
Solve:
Main function that runs the algorithm after initialization of the data
@pre: FileName containing the system to solve
@post: void for the moment
*/
void Solve(char *FileName)
{
	Problem *theProblem= createProblem(FileName);
	printf("OK: Problem created \n");
	Kruskal(theProblem);
	int i;
//	for(i=0; i< theProblem->nNode-1; i++)
//	{
//		printf("weight of edge %d in tree is %f \n", i, theProblem->theTree.edgesTree[i]->weight);
//	}

//	printf(" predecessor[0]=%d \n ",theProblem->theTree.predecessor[0]);
//	printf(" predecessor[1]=%d \n ",theProblem->theTree.predecessor[1]);
//	printf(" predecessor[2]=%d \n ",theProblem->theTree.predecessor[2]);
//	printf(" predecessor[3]=%d \n ",theProblem->theTree.predecessor[3]);

	//Chemin* findCycle(Edge *edgeCurrent,Problem *theProblem)

	//Chemin *wayTest=findCycle( &theProblem->edges[3], theProblem );

	//printf("ok way \n");
	//printf("wayTest.size= %d \n", wayTest->size);
//	for(i=0; i<wayTest->size; i++)
//	{
//		printf("edge[%d].indice = %d \n", i, wayTest->theChemin[i]->indice);
//	}
}

void setFlow(int indexNodeA, int indexNodeB, Problem *theProblem)
{
	Chemin* path= findPath(indexNodeA, indexNodeB, theProblem);
	int i;
	for(i=0; i< path->size; i++)
	{
		path->theChemin[i]->f=1.0;
	}// TODO il faut s'occuper du sens du flow I guess

}



		  /////////////////////////////////////////////////////////////////////////
		 //////////////////////// Initialisation /////////////////////////////////
		/////////////////////////////////////////////////////////////////////////

/*
createProblem
Reads the data file and initialises all the data structure needed during the computations to come.
@pre: The FileName
@post: Returns a pointeur to a Problem structure containing all the data.
*/
Problem *createProblem(char *FileName)
{

Problem *theProblem = new Problem[1];
int i,j,trash;

FILE* fichier = NULL;
int length;
fichier = fopen(FileName,"r+");

	if (fichier != NULL)
	{
		fseek(fichier,16,SEEK_SET);
		fscanf(fichier,"%d",&theProblem->nNode);
		int nNode=theProblem->nNode;

		fseek(fichier,17,SEEK_CUR);
		fscanf(fichier,"%d",&theProblem->nEdge);

		printf("Nombre d'edges %d\n",theProblem->nEdge);



		double **tableau = (double**) malloc(nNode*sizeof(double*));
		for (i = 0; i < nNode; i++)
		       tableau[i] = (double*) malloc(nNode*sizeof(double));


		for(int row = 0; row < nNode; row++)
		{
			for(int column = 0; column < nNode; column++)
			{
				fseek(fichier,1,SEEK_CUR);
				fscanf(fichier,"%lf",&tableau[row][column]);
				tableau[row][column]= - tableau[row][column];
				if(row == column)	tableau[row][row]=0.0;

			}
		}
		fclose(fichier);
		theProblem->Weights=tableau;
	}
	else
	{
		printf("Impossible d'ouvrir le fichier data.txt");
	}

	for(i=0;i<theProblem->nNode; i++)
	{
		for(j=0; j<theProblem->nNode; j++)
		{
			printf("%f", theProblem->Weights[i][j]);
		}
	printf("\n");

	}


//int nNode=theProblem->nNode;
theProblem->nodes=new Node[theProblem->nNode];
theProblem->edges=new Edge[theProblem->nEdge];
int nEdgeCurrent=0; // à mettre à jour !
int nVoisins=0;
int *nVoisinsCurrent= new int[theProblem->nNode];// nombre de voisins déjà ajouté au node[i]

// pre calcul du nombre de voisins= barbare :-)
for(i = 0; i<theProblem->nNode; i++)
{
	nVoisinsCurrent[i]=0;
	for(j=0; j<theProblem->nNode; j++)
	{
		if(theProblem->Weights[i][j]!=0.0)
		{
			nVoisins=nVoisins+1;
		}

	}
	// init node[i]
  	theProblem->nodes[i].indice=i;
  	theProblem->nodes[i].degree=nVoisins;
  	theProblem->nodes[i].voisins = new Node*[nVoisins];
 	theProblem->nodes[i].incidentes = new Edge*[nVoisins];
	nVoisins=0;
}

for(i = 0; i < theProblem->nNode; i++)
  {
      for(j = 0; j < i; j++)
      {
  	     if(theProblem->Weights[i][j]!=0.0)
  	     {

  	     	// create edge[i,j] OK
  	     	theProblem->edges[nEdgeCurrent].indice=nEdgeCurrent;
  	     	theProblem->edges[nEdgeCurrent].a=&theProblem->nodes[i];
  	     	theProblem->edges[nEdgeCurrent].b=&theProblem->nodes[j];
  	     	theProblem->edges[nEdgeCurrent].weight=theProblem->Weights[i][j];
  	     	theProblem->edges[nEdgeCurrent].f=0.0;
  	     	theProblem->edges[nEdgeCurrent].SET_FLOW=false;


  	     	// update node[i]
  	     	theProblem->nodes[i].voisins[nVoisinsCurrent[i]]=&theProblem->nodes[j];
  	     	theProblem->nodes[i].incidentes[nVoisinsCurrent[i]]=&theProblem->edges[nEdgeCurrent];
  	     	nVoisinsCurrent[i]++;

  	     	// update node[j]
  	     	theProblem->nodes[j].voisins[nVoisinsCurrent[j]]=&theProblem->nodes[i];
  	     	theProblem->nodes[j].incidentes[nVoisinsCurrent[j]]=&theProblem->edges[nEdgeCurrent];
  	     	nVoisinsCurrent[j]++;

  	     	nEdgeCurrent++;

  	     }
      }
  }

Tree arbreNew;
arbreNew.edgesTree=NULL;
theProblem->theTree= arbreNew;

return theProblem;

}

		  /////////////////////////////////////////////////////////////////////////
		 //////////////////////// KRUSKAL & Cie //////////////////////////////////
		/////////////////////////////////////////////////////////////////////////

/*int edgeCompare( const void * edgea, const void * edgeb)
{
	return (  (*(Edge*)edgea).weight- (*(Edge*)edgeb).weight);
}

void edgeSort(Problem *theProblem)
{
	qsort(theProblem->edges, theProblem->nEdge, sizeof(Edge), edgeCompare);
}*/


/*
function Kruskal
Applies Kruskal's algorithm to find a (min weight) spanning tree. Initialises the Tree structure of the Problem
@pre: The problem
@post: void
*/
void Kruskal(Problem *theProblem)
{
    int i,j,k,a,b,u,v,ne=1;
    int n = theProblem->nNode;
    int min,mincost=0;
    int Wmax=999;
    theProblem->theTree.edgesTree=new Edge*[theProblem->nNode-1];

  //  int *treeIndex= malloc((theProblem->nNode-1)*sizeof(int));
    //  Edge **treeIndex= new Edge*[theProblem->nNode-1];// ((theProblem->nNode-1)*sizeof(int));
    int compteur=0;
//    delete[] treeIndex;
//
//    x = new int; => delete x;
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
    		//treeIndex[compteur]=findIndex(a,b,theProblem);
    		theProblem->theTree.edgesTree[compteur]=&theProblem->edges[findIndex(a,b,theProblem)];
    		ne++;
    		compteur++;
    		mincost +=min;
    	}
    		theProblem->Weights[a][b]=theProblem->Weights[b][a]=Wmax;
    	}

    	/// Maintenant on rempli predecessor en faisant DFS
    	theProblem->theTree.nodeSource = new Node[1];
    	theProblem->theTree.predecessor= new int[theProblem->nNode-1];
    	theProblem->theTree.predecessor[0]=-1;// le premier node n'a pas de predecessor, on prend racine = node d'indice zero

    	queue<int> myqueue;
	myqueue.push(0);
    	int *added= new int[theProblem->nNode];
    	added[0]=1;
    	for(i=1; i<theProblem->nNode; i++) added[i]=0;

    	while(!myqueue.empty())
	{// il boucle à l'infini car on remet sur la stack des nodes qui ont déjà été dessus !
		int currentNode=myqueue.front();// on va spliter ce node
        	myqueue.pop();
        	for (i = 0; i < theProblem->nodes[currentNode].degree ; i++)
            	{
            		if( added[theProblem->nodes[currentNode].voisins[i]->indice]==0)
            		{
			myqueue.push(theProblem->nodes[currentNode].voisins[i]->indice);
			theProblem->theTree.predecessor[theProblem->nodes[currentNode].voisins[i]->indice]=currentNode;
            		added[theProblem->nodes[currentNode].voisins[i]->indice]=1;
            		}
            	}
    	}

}

/*
Function find
Needed in Kruskal's algorithm
@pre:
@post:
*/
int find(int i)
{
	while(parent[i])
	i=parent[i];
	return i;
}
/*
Function uni
Needed in Kruskal's algorithm
@pre:
@post:
*/
int uni(int i,int j)
{
    	if(i!=j)
    	{
    		parent[j]=i;
    		return 1;
    	}
    	return 0;
}

/*
Function findIndex
Returns the index of an edge that links 2 nodes. Given the index of the 2 nodes a and b
@pre: Index of 2 neighbour nodes
@post: Index of edge linking the two nodes
*/
int findIndex(int a, int b,Problem *theProblem)
{// TODO cette fonction est assez moche donc si à terme on s'en passe ça serait cool :-)
int i;
int index=0;
for(i=0; i < theProblem->nodes[a].degree; i++)
{
	if(theProblem->nodes[a].incidentes[i]->a->indice == b || theProblem->nodes[a].incidentes[i]->b->indice == b )
	{
		index = theProblem->nodes[a].incidentes[i]->indice;
		break;
	}
}

return index;
}
/////////////////////////////////////////////////////////////////////
/////////////////////////// DFS & Cie //////////////////////////////

/*
Renvoie le  chemin qui relie les nodes d'indices IndexNodeA et IndexNodeB
Attention, on suppose que l'edge A-B (si elle existe), n'est PAS dans l'arbre !!
*/
Chemin* findPath(int IndexNodeA, int IndexNodeB, Problem *theProblem)
{

	Chemin *path= new Chemin[1];

	path->theChemin = new Edge*[theProblem->nNode-1]; // faudra en enlever à la fin car il est trop long TODO

	int numberEdges = 0;
	Node *nodeA= &theProblem->nodes[nodeA];
	Node *nodeB= &theProblem->nodes[nodeB];

	int nextB=0;
	int currentB=IndexNodeA;
	int currentA=IndexNodeB;
	int *tabA=new int[theProblem->nNode-1];// tableau d'indice des nodes sur le chemin
	tabA[0]=currentA;
	int sizeB=0;
	int sizeA=1;

	int stopNode=-1;
	int i=0;

	int check=1;
	while(check)
	{
		int nextB=theProblem->theTree.predecessor[currentB];
		int nextA=theProblem->theTree.predecessor[currentA];

		path->theChemin[sizeB]=&theProblem->edges[findIndex(nextB, currentB,theProblem)];// add indice new edge en B
		sizeB++;
		// check tab en A: stop ou non
		tabA[sizeA]=nextA;
		sizeA++;
		for(i=0; i<sizeA; i++)
		{
			if(nextB == tabA[i])
			{
				check=0;
				stopNode=i;
				break;
			}
		}

		currentA=nextA;
		currentB=nextB;
	}

	int count=0;
	for(i=stopNode; i>=1; i--)
	{
		path->theChemin[sizeB+count]=&theProblem->edges[findIndex(tabA[i],tabA[i-1],theProblem)];
		count++;
	}
	delete[] tabA;
	path->size=count+sizeB;
	return path;

}


/*
Function findCycle
Returns the Chemin in the tree that goes from both ends of the edgeCurrent
@pre: A pointer to an off-tree Edge
@post: A pointer to a Chemin structure that goes from the Node edgeCurrent->b to the Node edgeCurrent->a
*/
Chemin* findCycle(Edge *edgeCurrent,Problem *theProblem)
{
	Chemin *path= new Chemin[1];

	path->theChemin = new Edge*[theProblem->nNode-1]; // faudra en enlever à la fin car il est trop long TODO
	Node *nodeA= edgeCurrent->a;
	Node *nodeB= edgeCurrent->b;

	if( (theProblem->theTree.predecessor[nodeA->indice]==nodeB->indice)
		||(theProblem->theTree.predecessor[nodeB->indice]==nodeA->indice) )
    {// si l'edge est dans le tree
        path->size=1;
        path->theChemin[0]=edgeCurrent;
        return path;
    }

	path=findPath(edgeCurrent->a->indice,edgeCurrent->b->indice, theProblem);
	return path;

}

////////////////////////// HAROLD'S PRATICE /////////////////////////////

// ON Continue Jean Pierre ! ///

/* A VERIFIER */
double stretchEdge(Edge *edgeCurrent, Chemin *Chemin)
{
    double stretch = 0.0;
    int length = Chemin->size;

    for(int i = 0; i < length; i++)
    {
        stretch = stretch + 1/(Chemin->theChemin[i]->weight);
    }
    stretch = stretch*(edgeCurrent->weight);

    return stretch;
}

/* A VERIFIER */
double stretchTree(Problem *theProblem)
{
    double stretch = 0.0;

    for(int i = 0; i < theProblem->nEdge; i++)
    {
        Chemin *theChemin = findCycle(&theProblem->edges[i], theProblem);
        stretch = stretch + stretchEdge(&theProblem->edges[i], theChemin);
    }
    return stretch;
}

/* A VERIFIER */
double probabilityEdge(Edge *edgeCurrent, Problem *theProblem)
{
    Chemin *theChemin = findCycle(edgeCurrent, theProblem);
    double stretchE = stretchEdge(edgeCurrent, theChemin);
    double stretchT = stretchTree(theProblem);
    double re = 1/(edgeCurrent->weight);
    double Re = re*(1+stretchE);
    int m = theProblem->nEdge;
    int n = theProblem->nNode;
    double CondNum = stretchT + m - 2*n + 2.0;
    double probability = (1/CondNum)*(Re/re);

    return probability;
}

/* A VERIFIER */
double* probaCompute(Problem *theProblem)
{
    int nEdgesOffTree = theProblem->nEdge - (theProblem->nNode-1);

    double *probabilities = new double[nEdgesOffTree];

    int edgesFound = 0;

    for(int i = 0; i < theProblem->nEdge; i++) // PROBABLEMENT OK mais un peu barbare (ça arrive :p)
    {
        bool entered = false;
        for(int j = 0; j < (theProblem->nNode-1); j++)
        {
            if(theProblem->edges[i].indice == theProblem->theTree.edgesTree[j]->indice)
            {
                entered = true;
                break;
            }
        }
        if(!entered)
        {
            probabilities[edgesFound] = probabilityEdge(&theProblem->edges[i], theProblem);
            edgesFound = edgesFound + 1;
        }
    }

    return probabilities;
}

/* A VERIFIER */
int iterationsK(Problem *theProblem, double eps)
{
    double stretchT = stretchTree(theProblem);
    int m = theProblem->nEdge;
    int n = theProblem->nNode;
    double CondNum = stretchT + m - 2*n + 2;
    int K = (int)ceil(CondNum*log(stretchT*CondNum/eps));

    return K;
}

/* A VERIFIER */
void CycleUpdate(Edge *edgeCurrent, Problem *theProblem)
{
    double Delta = (edgeCurrent->f)/(edgeCurrent->weight);

    Chemin *myChemin = findCycle(edgeCurrent, theProblem);
    Edge **edgesChemin = myChemin->theChemin;

    double stretchE = stretchEdge(edgeCurrent, myChemin);
    double re = 1/(edgeCurrent->weight);
    double Re = re*(1+stretchE);

    for(int i = 0; i < myChemin->size; i++)
    {
        Delta = Delta + (edgesChemin[i]->f)/(edgesChemin[i]->weight);
    }

    edgeCurrent->f = edgeCurrent->f - Delta/Re;

    int indiceCurrent = edgeCurrent->b->indice;
    for(int i = 0; i < myChemin->size; i++)
    {
        if(indiceCurrent == edgesChemin[i]->a->indice)
        {
            edgesChemin[i]->f = edgesChemin[i]->f - Delta/Re;
            indiceCurrent = edgesChemin[i]->b->indice;
        }
        else
        {
            edgesChemin[i]->f = edgesChemin[i]->f + Delta/Re;
            indiceCurrent = edgesChemin[i]->a->indice;
        }
    }
}

/* A VERIFIER */
Edge* RandomPicking(Problem *theProblem, Edge **edgesOffTree)
{
    int nEdgesOffTree = theProblem->nEdge - (theProblem->nNode - 1);
    double *probabilities = probaCompute(theProblem);

    double *cumulatedProba = new double[nEdgesOffTree+1];
    cumulatedProba[0] = 0.0;
    for(int i = 1; i <= nEdgesOffTree; i++)
    {
        cumulatedProba[i] = cumulatedProba[i-1] + probabilities[i-1];
    }

    double value = ((double)rand() / ((double)RAND_MAX+1)); // Pas hyper top apparemment mais peut-etre suffisant ici

    for(int i = 1; i <= nEdgesOffTree; i++)
    {
        if( (cumulatedProba[i-1] <= value) && (value < cumulatedProba[i]) )
        {
            return edgesOffTree[i-1];
        }
    }
    // Je dois ajouter un delete qqpart ?
}

/*
	Tu as edgeCurrent qui t'as été donnée par la fonciton précédente. Tu calcule son stretch et tout si necessaire,
	son cycle via la fonction
	Chemin* DFS(Edge *edgeCurrent,Problem *theProblem);

	qui renvoie:
Specification:
Returns a Chemin structure that goes from the Node "edgeCurrent->b" to "edgeCurrent->a".
edgeCurrent is an edge that does not belong to theTree.


Un object/structure de type Chemin t'as cela sur slack si jamais :-)

Donc ici ton job c'est de parcourir toutes les edges du chemin (via le tableau d'edges de la struct Chemin) et de mettre le flot à jour.
Il est possible que la fonction ne fasse que 15 lignes, c'est un peu celle où on utilise quelques précédentes :-)
*/


// GO GO GO
//}




