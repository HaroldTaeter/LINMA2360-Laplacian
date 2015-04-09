/*
 *  laplacian.cpp
 *
 */

//#define _PRINT_SOL_
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

	int i,j;
	Problem *theProblem= createProblem(FileNameL, FileNameb);
	printf("OK: Problem created \n");


	struct timespec startIN, finishIN;
       double elapsed;

       clock_gettime(CLOCK_MONOTONIC, &startIN);

	Kruskal(theProblem);
	printf("OK: Kruskal \n");

	SetNodes(theProblem);
	printf("OK: SetNodes \n");

	SetPredecessor(theProblem);
	printf("OK: Predecessor \n");

	stretchsAndChemins(theProblem);
	printf("OK: Strech and Chemins \n");

			// RUN
	setFlow(0, theProblem->nNode-1, theProblem);
	printf("OK: Set flow \n");

	theProblem->theCumulatedProba=probaCompute(theProblem);

	int iter;
	int iterMax= iterationsK(theProblem, 5.1);
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
	voltages = InducedVoltages2(theProblem);

	clock_gettime(CLOCK_MONOTONIC, &finishIN);

	elapsed = (finishIN.tv_sec - startIN.tv_sec);
	elapsed += (finishIN.tv_nsec - startIN.tv_nsec) / 1000000000.0;
	printf("time_elapsed_computation: %f\n",elapsed);

#ifdef _PRINT_SOL_
    for(int i = 0; i < theProblem->nNode; i++)
    {
        printf("Voltage on node %d = %f \n",i,voltages[i]);
    }
#endif

	FreeLaplacian(theProblem);

}

void FreeLaplacian(Problem *theProblem)
{
// TODO il y a plein d'autres trucs à free !!
delete [] theProblem;
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
int node1;
int node2;

	if (fichier == NULL) printf("Error: No data file !\n");

		fscanf(fichier,"%d \n",&theProblem->nNode);
		int nNode = theProblem->nNode;
		printf("Nombre de nodes %d\n",theProblem->nNode);

		fscanf(fichier,"%d \n",&theProblem->nEdge);
		printf("Nombre d'edges %d\n",theProblem->nEdge);

		theProblem->nodes=new Node[theProblem->nNode];
		theProblem->edges=new Edge[theProblem->nEdge];
		for(i=0; i<nNode; i++)
		{// TODO besoin de cette boucle ? ou elle peut aller dans SetPredecessor ?
		// je crois qu'on a besoin des indices avant predecessor (dans kruskal)
			theProblem->nodes[i].degree=0;
			theProblem->nodes[i].indice=i;
		}

		for(i=0; i<theProblem->nEdge; i++)
		{
			fscanf(fichier,"%d %d %lf\n",&node1, &node2, &theProblem->edges[i].weight);

			theProblem->edges[i].a=&theProblem->nodes[node1];
			theProblem->edges[i].b=&theProblem->nodes[node2];
			theProblem->edges[i].indice=i;
			theProblem->edges[i].f=0.0;// toutes les edges sont ok :-)
		}

		fclose(fichier);
Tree arbreNew;
arbreNew.edgesTree=NULL;
arbreNew.predecessor=NULL;
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
		fclose(fichier2);
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
	fill_n(theProblem->parent, theProblem->nNode, 0);// fonctionne mais pourquoi ?
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
    		theProblem->theTree.edgesTree[compteurTree]=&theProblem->edgesSorted[iter];
    		//theProblem->theTree.edgesTree[compteurTree]=&theProblem->edges[findIndex(a,b,theProblem)];
    		// est ce que findIndex ne renvoie pas edgesSorted[iter] par hasard ? :D TODO vérifier !!
    		ne++;
    		compteurTree++;
    		}
    		else{
    		theProblem->edgesOffTree[compteurOffTree]=&theProblem->edgesSorted[iter];
    		//theProblem->edgesOffTree[compteurOffTree]=&theProblem->edges[findIndex(a,b,theProblem)];
    		// il faudrait remplacer le findIndex ici car on a supprimé les incidentes hors Tree :p
    		// peut se faire à la bourin en parcourant toutes les edges :p
    		compteurOffTree++;
	    	}
    		if (ne==theProblem->nNode){ break;}
	}
	    for(iter=ne-1; iter< theProblem->nEdge; iter++)
	    {
	    	theProblem->edgesOffTree[compteurOffTree]=&theProblem->edgesSorted[iter];
    		compteurOffTree++;
	    }
}

/*
Initialises the Nodes structure but only based on the Tree edges.
*/
void SetNodes(Problem *theProblem)
{
 int i,j,node1,node2;;
 int *nVoisinsCurrent= new int[theProblem->nNode];// nombre de voisins déjà ajouté au node[i]

for(i=0; i<theProblem->nNode-1; i++)// precalcul du nombre de voisins (serait cool de l'éviter)
{
	node1=theProblem->theTree.edgesTree[i]->a->indice;
	node2=theProblem->theTree.edgesTree[i]->b->indice;

	theProblem->nodes[node1].degree++;
	theProblem->nodes[node2].degree++;
}

for(i=0; i<theProblem->nNode; i++)
{// initialisation pour chaque noeud
	nVoisinsCurrent[i]=0;
	theProblem->nodes[i].voisins= new Node*[theProblem->nodes[i].degree];
	theProblem->nodes[i].incidentes= new Edge*[theProblem->nodes[i].degree];
}
for(i=0; i<theProblem->nNode-1; i++)// boucle sur edgesTree
{
	node1=theProblem->theTree.edgesTree[i]->a->indice;
	node2=theProblem->theTree.edgesTree[i]->b->indice;
	// update node 1

	theProblem->nodes[node1].incidentes[nVoisinsCurrent[node1]]=theProblem->theTree.edgesTree[i];
	theProblem->nodes[node1].voisins[nVoisinsCurrent[node1]]=&theProblem->nodes[node2];
	nVoisinsCurrent[node1]++;
	// update node 2

	theProblem->nodes[node2].incidentes[nVoisinsCurrent[node2]]=theProblem->theTree.edgesTree[i];
	theProblem->nodes[node2].voisins[nVoisinsCurrent[node2]]=&theProblem->nodes[node1];
	nVoisinsCurrent[node2]++;
}
delete[] nVoisinsCurrent;
}


/*
Using the edgesTree, computes the predecessor array from the Tree (used in findPath)
*/
void SetPredecessor(Problem *theProblem)
{
	int i;
    	/// Maintenant on rempli predecessor en faisant DFS
	theProblem->theTree.nodeSource = &theProblem->nodes[0];// smart choice for source ? node central ?
    	theProblem->theTree.predecessor= new int[theProblem->nNode];
    	theProblem->theTree.predecessor[0]=-1;// le premier node n'a pas de predecessor, on prend racine = node d'indice zero

    	queue<int> myqueue;
	myqueue.push(theProblem->theTree.nodeSource->indice);
    	int *visited= new int[theProblem->nNode];
    	fill_n(visited, theProblem->nNode,0);
    	visited[theProblem->theTree.nodeSource->indice]=1;

    	while(!myqueue.empty())
	{
		int currentNode=myqueue.front();// on va spliter ce node
        	myqueue.pop();
        	for (i = 0; i < theProblem->nodes[currentNode].degree ; i++)
            	{// TODO check
            		if(visited[theProblem->nodes[currentNode].voisins[i]->indice] == 0)
            		{
            		myqueue.push(theProblem->nodes[currentNode].voisins[i]->indice);
            		theProblem->theTree.predecessor[theProblem->nodes[currentNode].voisins[i]->indice]=currentNode;
	            	visited[theProblem->nodes[currentNode].voisins[i]->indice]=1;
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
@post: Index of edge linking the two nodes (edge in the Tree)
*/
int findIndex(int a, int b,Problem *theProblem)
{// je crois que si on sait enlever cette fonction c'est en utilisant des pures listes d'adjacences mais ça revient un peu au même
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

     	Chemin path;
     	path.size=1;
     	path.theChemin=new Edge*[1];
     	path.theChemin[0]=theProblem->theTree.edgesTree[i];
     	theProblem->theTree.edgesTree[i]->edgeChemin=&path;
    }
    stretchTree = stretchTree + (theProblem->nNode - 1);

    theProblem->theTree.stretch = stretchTree;
}



////////////////////////// HAROLD'S PRATICE /////////////////////////////

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

    /*for(int i = 1; i <= nEdgesOffTree; i++)
    {
        if( (theProblem->theCumulatedProba[i-1] <= value) && (value < theProblem->theCumulatedProba[i]) )
        {
            return edgesOffTree[i-1];
        }
    }*/

    int beginIndex = 0;
    int endIndex = nEdgesOffTree;
    int midIndex = (int)round(0.5*(beginIndex+endIndex));
    while(beginIndex != endIndex-1)
    {
        if(theProblem->theCumulatedProba[midIndex] > value)
        {
            endIndex = midIndex;
            midIndex = (int)round(0.5*(beginIndex+endIndex));
        }
        else if(theProblem->theCumulatedProba[midIndex] <= value)
        {
            beginIndex = midIndex;
            midIndex = (int)round(0.5*(beginIndex+endIndex));
        }
    }
    return edgesOffTree[beginIndex];
}

/*
Semble OK: Renvoie le vecteur des voltages aux nodes sur base des flots.
*/
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
                if(edgeCurrent.a->indice==currentNode)
                {// sens de l'edge pour le courant
                voltages[theProblem->nodes[currentNode].voisins[i]->indice] = voltages[currentNode] - (edgeCurrent.f)/(edgeCurrent.weight);
                }
                else{
                voltages[theProblem->nodes[currentNode].voisins[i]->indice] = voltages[currentNode] + (edgeCurrent.f)/(edgeCurrent.weight);
                }

                visited[theProblem->nodes[currentNode].voisins[i]->indice] = 1;
            }
        }
    }

    return voltages;
}


//void setFlow2(Problem *theProblem, double *ksi)
//{
//// TODO stoquer les leaves quand on fait DFS avec Kruskal
//    Tree myTree = theProblem->theTree;
//    int nLeaves = 0;
//    Node **leaves = new Node*[theProblem->nNode]; // TROP GRAND
//    for(int i = 0; i<theProblem->nNode; i++)
//    {
//        if(theProblem->nodes[i]->degree == 1) // LA STRUCTURE DOIT CHANGER POUR AVOIR degreeTree
//        {					// le degree du node avec les edges de l'arbre
//            leaves[nLeaves] = theProblem->nodes[i];
//            nLeaves++;
//        }
//    }

//    double *edgesDone = new double[theProblem->nNode];
//    for(int i = 0; i < theProblem->nNode; i++)// nNode ?
//    {// edge qui mène au predecessor ?
//        edgesDone[i] = 0;
//    }

//    int nNewLeaves;
//    while(nLeaves > 0)
//    {
//        nNewLeaves = 0;
//        Node **Newleaves = new Node*[theProblem->nNode]; // TROP GRAND

//        for(int i = 0; i < nLeaves; i++)
//        {
//            int edgeIndex = findIndex(leaves[i]->indice,myTree.predecessor[leaves[i]->indice],theProblem);
//            if(theProblem->edges[edgeIndex].b->indice == leaves[i]->indice)
//            {
//                theProblem->edges[edgeIndex].f += ksi[leaves[i]->indice];
//                ksi[myTree.predecessor[leaves[i]->indice]] += -theProblem->edges[edgeIndex].f;
//            }
//            else if(theProblem->edges[edgeIndex].a->indice == leaves[i]->indice)
//            {
//                theProblem->edges[edgeIndex].f += -ksi[leaves[i]->indice];
//                ksi[myTree.predecessor[leaves[i]->indice]] += theProblem->edges[edgeIndex].f;
//            }
//            edgesDone[myTree.predecessor[leaves[i]->indice]] += 1;

//            if(theProblem->nodes[myTree.predecessor[leaves[i]->indice]].degreeTree - edgesDone[myTree.predecessor[leaves[i]->indice]] == 1
//               && myTree.predecessor[leaves[i]->indice] != -1)
//            {
//                Newleaves[nNewLeaves] = theProblem->nodes[myTree.predecessor[leaves[i]->indice]];
//                nNewLeaves++;
//            }
//        }
//        nLeaves = nNewLeaves;
//        leaves = Newleaves;

//        delete[] Newleaves;
//    }
//}


/////////////////////////////////////////////////////////////////////
/////////////////////////// Section 5 ///////////////////////////////




