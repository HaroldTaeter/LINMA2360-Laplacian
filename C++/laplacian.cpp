/*
 *  laplacian.cpp
 *
 */

#include "laplacian.h"

using namespace std;

/*
Solve:
Main function that runs the algorithm after initialization of the data
@pre: FileName containing the system to solve
@post: void for the moment
*/
void Solve(char *FileNameL,char *FileNameb)
{
	
	int i;
	Problem *theProblem= createProblem(FileNameL, FileNameb);
	printf("OK: Problem created \n");

	Kruskal(theProblem);	
	printf("OK: Kruskal \n");
	
	stretchsAndChemins(theProblem);
	printf("OK: Strech and Chemins \n");
	
			// RUN
	setFlow(0, theProblem->nNode-1, theProblem);
	printf("OK: Set flow \n");
		
	theProblem->theCumulatedProba=probaCompute(theProblem);

	int iter;
	int iterMax= iterationsK(theProblem, 0.001);
	printf("iterMax=%d \n", iterMax);
	for(iter=0; iter<iterMax; iter++)
	{
		//printf(" //////// Iteration number: %d /////////////\n", iter);
		Edge *edgeCurrent = RandomPicking(theProblem, theProblem->edgesOffTree);
		//printf(" edgeCurrent->indice = %d \n", edgeCurrent->indice);
		CycleUpdate(edgeCurrent, theProblem);
		//double Delta=testKPL(edgeCurrent);
		//printf("Delta after cycleUpdate = %f \n", Delta);
	}
	
    double *voltages = new double[theProblem->nNode];

    voltages = InducedVoltages(theProblem);

    for(int i = 0; i < theProblem->nNode; i++)
    {
        printf("Voltage on node %d = %f \n",i,voltages[i]);
    }
}

void setFlow(int indexNodeA, int indexNodeB, Problem *theProblem)
{
	Chemin* path= findPath(indexNodeA, indexNodeB, theProblem);
	int i;
	int indiceCurrent = indexNodeA;
	for(i = 0; i < path->size; i++)
	{
		if(indiceCurrent == path->theChemin[i]->a->indice)
		{
			path->theChemin[i]->f = 1.0;
			indiceCurrent = path->theChemin[i]->b->indice;
		}
		else
		{
			path->theChemin[i]->f = -1.0;
			indiceCurrent = path->theChemin[i]->a->indice;
		}
	}
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
Problem *createProblem(char *FileNameL, char *FileNameb)
{

Problem *theProblem = new Problem[1];
int i,j,trash;

FILE* fichier = NULL;
int length;
fichier = fopen(FileNameL,"r+");

	if (fichier != NULL)
	{
		fseek(fichier,16,SEEK_SET); // WINDOWS : 16
		fscanf(fichier,"%d",&theProblem->nNode);
		int nNode = theProblem->nNode;
		printf("Nombre de nodes %d\n",theProblem->nNode);

		fseek(fichier,17,SEEK_CUR); // WINDOWS : 19
		fscanf(fichier,"%d",&theProblem->nEdge);

		printf("Nombre d'edges %d\n",theProblem->nEdge);
// OK ??? TODO
//		fseek(fichier,1,SEEK_CUR);

		double **tableau = (double**) malloc(nNode*sizeof(double*));
		for (i = 0; i < nNode; i++)
		       tableau[i] = (double*) malloc(nNode*sizeof(double));

		for(int row = 0; row < nNode; row++)
		{
			for(int column = 0; column < nNode; column++)
			{
				fseek(fichier,1,SEEK_CUR);
				fscanf(fichier,"%lf",&tableau[row][column]);
				//printf("%f ", tableau[row][column]);
				tableau[row][column]= -tableau[row][column];
				if(row == column)	tableau[row][row]=0.0;
			}
			//printf("\n");
		}
		
		fclose(fichier);
		printf("\n");

		theProblem->Weights=tableau;
	}
	else
	{
		printf("Impossible d'ouvrir le fichier dataL.txt");
	}

//	for(i=0;i<theProblem->nNode; i++)
//	{
//		for(j=0; j<theProblem->nNode; j++)
//		{
//			printf(" %f ", theProblem->Weights[i][j]);
//		}
//	printf("\n");
//	}


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
arbreNew.nodeSource = &theProblem->nodes[0];
arbreNew.stretch=0.0;
theProblem->theTree= arbreNew;

theProblem->theCumulatedProba=NULL;
FILE* fichier2 = NULL;
fichier2 = fopen(FileNameb,"r+");
theProblem->b= new double[theProblem->nNode];

	if (fichier2 != NULL)
	{
		for(i=0; i<theProblem->nNode; i++)
		{
			fscanf(fichier2,"%lf \n",&theProblem->b[i]);
		}
	}
	else{printf("Impossible d'ouvrir le fichier datab.txt\n");}
	
theProblem->edgesSorted=new Edge[theProblem->nEdge];
theProblem->edgesOffTree = new Edge*[theProblem->nEdge-(theProblem->nNode-1)];
theProblem->parent = new int[theProblem->nNode];
return theProblem;
}

		  /////////////////////////////////////////////////////////////////////////
		 //////////////////////// KRUSKAL & Cie //////////////////////////////////
		/////////////////////////////////////////////////////////////////////////

int edgeCompare( const void * edgea, const void * edgeb)
{
	return (  -(*(Edge*)edgea).weight+ (*(Edge*)edgeb).weight);// max weight spanning tree
}
void edgeSort(Problem *theProblem)
{// TODO on pourrait avoir le tableau edgesSorted qui soit des pointeurs..
// mais ça rend les cast bien bordeliques donc pas pour tout de suite !!
	qsort(theProblem->edgesSorted, theProblem->nEdge, sizeof(Edge), edgeCompare);
}


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
   	theProblem->theTree.edgesTree=new Edge*[theProblem->nNode-1];

	for(i=0; i<theProblem->nEdge;i++)
	{
	theProblem->edgesSorted[i]=theProblem->edges[i]; 
	}
    	int compteurTree=0;
	int compteurOffTree=0;
	edgeSort(theProblem);	
        int iter;	
	for(iter=0; iter<theProblem->nEdge; iter++)
	{
		i=theProblem->edgesSorted[iter].a->indice;
		j=theProblem->edgesSorted[iter].b->indice;
		a=u=i;
		b=v=j;
		u=find(u, theProblem);
	    	v=find(v, theProblem);
    		if(uni(u,v, theProblem))
    		{	
    		theProblem->theTree.edgesTree[compteurTree]=&theProblem->edges[findIndex(a,b,theProblem)];
    		ne++;
    		compteurTree++;
    		}	
    		else{
    		theProblem->edgesOffTree[compteurOffTree]=&theProblem->edges[findIndex(a,b,theProblem)];
    		compteurOffTree++;
	    	}
    		if (ne==theProblem->nNode){ break;} 
	}
	    for(iter=ne-1; iter< theProblem->nEdge; iter++)
	    {
	    	theProblem->edgesOffTree[compteurOffTree]=&theProblem->edgesSorted[iter];
    		compteurOffTree++;
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
            		if( theProblem->nodes[currentNode].incidentes[i]->a->indice != currentNode  &&
            			added[theProblem->nodes[currentNode].incidentes[i]->a->indice] == 0 &&
            			inTree(theProblem->nodes[currentNode].incidentes[i], theProblem)==1 )
            		{
            			myqueue.push(theProblem->nodes[currentNode].incidentes[i]->a->indice);
            			theProblem->theTree.predecessor[theProblem->nodes[currentNode].incidentes[i]->a->indice]=currentNode;
	            		added[theProblem->nodes[currentNode].incidentes[i]->a->indice]=1;	
            		}
            		else if( theProblem->nodes[currentNode].incidentes[i]->b->indice != currentNode  &&
            			added[theProblem->nodes[currentNode].incidentes[i]->b->indice] == 0 &&
            			inTree(theProblem->nodes[currentNode].incidentes[i], theProblem)==1)
            		{
	      			myqueue.push(theProblem->nodes[currentNode].incidentes[i]->b->indice);
            			theProblem->theTree.predecessor[theProblem->nodes[currentNode].incidentes[i]->b->indice]=currentNode;
	            		added[theProblem->nodes[currentNode].incidentes[i]->b->indice]=1;
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
int find(int i, Problem *theProblem)
{
	while(theProblem->parent[i])
	i=theProblem->parent[i];
	return i;
}
/*
Function uni
Needed in Kruskal's algorithm
@pre:
@post:
*/
int uni(int i,int j, Problem *theProblem)
{
    	if(i!=j)
    	{
    		theProblem->parent[j]=i;
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

/*
Prend un pointeur sur une edge et renvoie 1 si elle est dans theTree, 0 sinon

*/
int inTree(Edge *edgea,Problem *theProblem)
{
	int i; 
	for(i=0; i<theProblem->nNode-1; i++)
	{
		if( theProblem->theTree.edgesTree[i]->indice == edgea->indice )
		{
			return 1; 
		}
	}
	
	return 0; 


}
/////////////////////////////////////////////////////////////////////
/////////////////////////// DFS & Cie //////////////////////////////

/*
Renvoie le  chemin qui relie les nodes d'indices IndexNodeA et IndexNodeB
Attention, on suppose que l'edge A-B (si elle existe), n'est PAS dans l'arbre !!
*/
Chemin* findPath(int IndexNodeA, int IndexNodeB, Problem *theProblem)
{
	Node *nodeA = &theProblem->nodes[IndexNodeA];
	Node *nodeB = &theProblem->nodes[IndexNodeB];

	int nextB = 0;
	int nextA = 0;
	int currentB = IndexNodeB;
	int currentA = IndexNodeA;
	int *tabA = new int[theProblem->nNode-1];// tableau d'indice des nodes sur le chemin
	int *tabB = new int[theProblem->nNode-1];// tableau d'indice des nodes sur le chemin
	tabA[0] = currentA;
	tabB[0] = currentB;
	int sizeB = 1;
	int sizeA = 1;

	int stopNodeA = -1;
	int stopNodeB = -1;
	int i = 0;

	int check = 1;

	while(check!=0)
	{
		int nextB = theProblem->theTree.predecessor[currentB];
		int nextA = theProblem->theTree.predecessor[currentA];

		if(nextA==-1 && nextB==-1)
        	{
            check = 0;
            stopNodeA = sizeA-1;
            stopNodeB = sizeB-1;
       		}
		else if(nextA==-1)
		{
			tabB[sizeB] = nextB;
			sizeB++;

			for(i=0; i<sizeA; i++)
			{
				if(nextB == tabA[i])
				{
					check = 0;
					stopNodeA = i;
					stopNodeB = sizeB-1;
					break;
				}
			}
			currentB = nextB;
		}
		else if(nextB==-1)
		{
			tabA[sizeA] = nextA;
			sizeA++;

			for(i=0; i<sizeB; i++)
			{
				if(nextA == tabB[i])
				{
					check = 0;
					stopNodeA = sizeA-1;
					stopNodeB = i;
					break;
				}
			}
			currentA = nextA;
		}
		else
        {
            tabB[sizeB] = nextB;
            sizeB++;
            tabA[sizeA] = nextA;
            sizeA++;

            for(i=0; i<sizeA; i++)
            {
                if(nextB == tabA[i])
                {
                    check = 0;
					stopNodeA = i;
					stopNodeB = sizeB-1;
					break;
                }
            }
            for(i=0; i<sizeB; i++)
            {
                if(nextA == tabB[i])
                {
                    check = 0;
					stopNodeA = sizeA-1;
					stopNodeB = i;
					break;
                }
            }

            currentA = nextA;
            currentB = nextB;
		}
		//printf("check... \n");
	}

	Chemin *path= new Chemin[1];
	path->theChemin = new Edge*[stopNodeA+stopNodeB];

	for(i = 1; i <= stopNodeB; i++)
    {
        path->theChemin[i-1] = &theProblem->edges[findIndex(tabB[i],tabB[i-1],theProblem)];
    }
    for(i = stopNodeA; i >= 1; i--)
    {
        path->theChemin[stopNodeB+(stopNodeA-i)] = &theProblem->edges[findIndex(tabA[i],tabA[i-1],theProblem)];
    }
    path->size = stopNodeA+stopNodeB;
    delete[] tabA;
    delete[] tabB;
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


/*
For an edge offTree, returns the value of Delta on the cycle ()
*/
double testKPL(Edge *edgeCurrent)
{
    double Delta = (edgeCurrent->f)/(edgeCurrent->weight);
    Edge **edgesChemin = edgeCurrent->edgeChemin->theChemin;  
    int indiceCurrent = edgeCurrent->b->indice;
    for(int i = 0; i < edgeCurrent->edgeChemin->size; i++)
    {   
        if(indiceCurrent == edgesChemin[i]->a->indice)
        {        
        Delta = Delta + (edgesChemin[i]->f)/(edgesChemin[i]->weight);
    	indiceCurrent = edgesChemin[i]->b->indice;
    	}
    	else
    	{
    	Delta = Delta - (edgesChemin[i]->f)/(edgesChemin[i]->weight);
    	indiceCurrent = edgesChemin[i]->a->indice;
    	}
    }
    
    return Delta;
}

void stretchsAndChemins(Problem *theProblem)
{
    double stretchTree = 0.0;

    for(int i = 0; i < (theProblem->nEdge - (theProblem->nNode - 1)); i++)
    {
        theProblem->edgesOffTree[i]->edgeChemin = findCycle(theProblem->edgesOffTree[i], theProblem);
        theProblem->edgesOffTree[i]->stretch = stretchEdge(theProblem->edgesOffTree[i], theProblem->edgesOffTree[i]->edgeChemin);
        stretchTree = stretchTree + theProblem->edgesOffTree[i]->stretch;
    }

    for(int i = 0; i < (theProblem->nNode - 1); i++)
    {
        theProblem->theTree.edgesTree[i]->stretch = 1;
     	
     	Chemin path;//= new Chemin[1];
     	path.size=1; 
     	path.theChemin=new Edge*[1];
     	path.theChemin[0]=theProblem->theTree.edgesTree[i];
     	theProblem->theTree.edgesTree[i]->edgeChemin=&path;// tout ça est fait un peu vite mais ma version devraitetre ok
        //theProblem->theTree.edgesTree[i]->edgeChemin->size = 1;// les trucs en commentaire c'est ta version, A Enlever normalement!
        //theProblem->theTree.edgesTree[i]->edgeChemin->theChemin[0] = theProblem->theTree.edgesTree[i];
    }
    stretchTree = stretchTree + (theProblem->nNode - 1);

    theProblem->theTree.stretch = stretchTree;
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
double probabilityEdge(Edge *edgeCurrent, Problem *theProblem)
{
    double stretchE = edgeCurrent->stretch;
    double stretchT = theProblem->theTree.stretch;    
    double re = 1/(edgeCurrent->weight);
    double Re = re*(1+stretchE);
    int m = theProblem->nEdge;
    int n = theProblem->nNode;
    double CondNum = stretchT + m - 2.0*n + 2.0;
    double probability = (1/CondNum)*(Re/re);

    return probability;
}

/* A VERIFIER */
double* probaCompute(Problem *theProblem)
{
    int nEdgesOffTree = theProblem->nEdge - (theProblem->nNode-1);

    double *probabilities = new double[nEdgesOffTree];
    double sumProba=0.0;
    int i; 
    for(i = 0; i < nEdgesOffTree; i++) // les probas correspondent à l'ordre de edgesOffTree
    {
            probabilities[i] = probabilityEdge(theProblem->edgesOffTree[i],theProblem);
            sumProba +=     probabilities[i];   
    }
    
    for(i = 0; i < nEdgesOffTree; i++) // TODO problem a régler ! La somme ne vaut pas 1 donc j'ai normalisé.
    {
            probabilities[i] = probabilities[i]/sumProba;

    }
	
//	for(i=0; i<nEdgesOffTree; i++)
//	{
//		printf(" edges number %d off the tree has indice= %d and proba=%f \n",i,theProblem->edgesOffTree[i]->indice, 
//			probabilities[i]);
//	}
//	
	
    double *cumulatedProba = new double[nEdgesOffTree+1];
    cumulatedProba[0] = 0.0;
    for(int i = 1; i <= nEdgesOffTree; i++)
    {
        cumulatedProba[i] = cumulatedProba[i-1] + probabilities[i-1];
    }
    return cumulatedProba;
}

/* A VERIFIER */
int iterationsK(Problem *theProblem, double eps)
{
    double stretchT = theProblem->theTree.stretch;   
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
    Edge **edgesChemin = edgeCurrent->edgeChemin->theChemin;
	
    double stretchE = edgeCurrent->stretch;
    double re = 1/(edgeCurrent->weight);
    double Re = re*(1+stretchE);
    
    int indiceCurrent = edgeCurrent->b->indice;
    for(int i = 0; i < edgeCurrent->edgeChemin->size; i++)
    {// Attention ici il faut prendre le sens du flot selon le 
    // sens de parcourt du cycle :-)        
        if(indiceCurrent == edgesChemin[i]->a->indice)
        {        
        Delta = Delta + (edgesChemin[i]->f)/(edgesChemin[i]->weight);
    	indiceCurrent = edgesChemin[i]->b->indice;
    	}
    	else
    	{
    	Delta = Delta - (edgesChemin[i]->f)/(edgesChemin[i]->weight);
    	indiceCurrent = edgesChemin[i]->a->indice;
    	}
    
    }
    if(Delta != 0.0)
    {
    edgeCurrent->f = edgeCurrent->f - Delta/Re;

    indiceCurrent = edgeCurrent->b->indice;
    for(int i = 0; i < edgeCurrent->edgeChemin->size; i++)
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
}

/* A VERIFIER */
Edge* RandomPicking(Problem *theProblem, Edge **edgesOffTree)
{
    int nEdgesOffTree = theProblem->nEdge - (theProblem->nNode - 1);

    double value = ((double)rand() / ((double)RAND_MAX+1)); // Pas hyper top apparemment mais peut-etre suffisant ici

    for(int i = 1; i <= nEdgesOffTree; i++)
    {
        if( (theProblem->theCumulatedProba[i-1] <= value) && (value < theProblem->theCumulatedProba[i]) )
        {
            return edgesOffTree[i-1];
        }
    }
}

/* A VERIFIER */
double* InducedVoltages(Problem *theProblem) // pas encore optimise
{
    double *voltages = new double[theProblem->nNode];
    for(int i = 0; i < theProblem->nNode; i++)
    {
        voltages[i] = 0.0;
    }

    for(int i = 0; i < theProblem->nNode; i++)
    {
        Chemin *myChemin = findPath(theProblem->theTree.nodeSource->indice, theProblem->nodes[i].indice, theProblem);
        // TODO: controler que findPath fait pas de la merde quand on lui envoie des indices proches etc
        for(int j = 0; j < myChemin->size; j++)
        {
            voltages[i] = voltages[i] + (myChemin->theChemin[j]->f)/(myChemin->theChemin[j]->weight);
        }
    }
    return voltages;
}

/* A VERIFIER */
double* InducedVoltages2(Problem *theProblem)
{
    double *voltages = new double[theProblem->nNode];
    for(int i = 0; i < theProblem->nNode; i++)
    {
        voltages[i] = 0.0;
    }

    Tree finalTree = theProblem->theTree; 
    Node *source = finalTree.nodeSource;

    queue<int> myqueue;
	myqueue.push(source->indice);
    int *visited = new int[theProblem->nNode];
    visited[0] = 1;
    for(int i = 1; i < theProblem->nNode; i++)
    {
        visited[i] = 0;
    }

    while(!myqueue.empty())
	{
		int currentNode = myqueue.front();
       	myqueue.pop();
       	for (int i = 0; i < theProblem->nodes[currentNode].degree ; i++)
        {
            if(visited[theProblem->nodes[currentNode].voisins[i]->indice]==0)
            {
                myqueue.push(theProblem->nodes[currentNode].voisins[i]->indice);
                Edge edgeCurrent = theProblem->edges[findIndex(currentNode, theProblem->nodes[currentNode].voisins[i]->indice, theProblem)];
                voltages[theProblem->nodes[currentNode].voisins[i]->indice] = voltages[currentNode] + (edgeCurrent.f)/(edgeCurrent.weight);
                visited[theProblem->nodes[currentNode].voisins[i]->indice] = 1;
            }
        }
    }

    return voltages;
}
