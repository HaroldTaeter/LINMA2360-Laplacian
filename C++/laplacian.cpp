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
	Kruskal(theProblem);
	int i; 
	for(i=0; i< theProblem->nNode-1; i++)
	{
		printf("weight of edge %d in tree is %f \n", i, theProblem->theTree.edgesTree[i]->weight);
	}
//	printf(" predecessor[0]=%d \n ",theProblem->theTree.predecessor[0]);
//	printf(" predecessor[1]=%d \n ",theProblem->theTree.predecessor[1]);
//	printf(" predecessor[2]=%d \n ",theProblem->theTree.predecessor[2]);
//	printf(" predecessor[3]=%d \n ",theProblem->theTree.predecessor[3]);
	
	//Chemin* findCycle(Edge *edgeCurrent,Problem *theProblem)
	
	Chemin *wayTest=findCycle( &theProblem->edges[3], theProblem );
	printf("ok way \n");
	printf("wayTest.size= %d \n", wayTest->size);
	for(i=0; i<wayTest->size; i++)
	{
		printf("edge[%d].indice = %d \n", i, wayTest->theChemin[i]->indice);
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
		
		printf("La longueur est %d\n",nNode);

		double **tableau = (double**) malloc(nNode*sizeof(double*));
		for (i = 0; i < nNode; i++)
		       tableau[i] = (double*) malloc(nNode*sizeof(double));
		

		for(int row = 0; row < nNode; row++)
		{
			for(int column = 0; column < nNode; column++)
			{
				fseek(fichier,1,SEEK_CUR);
				fscanf(fichier,"%lf",&tableau[row][column]); 
			}
		}
		fclose(fichier);
		theProblem->Weights=tableau;
	}
	else
	{
		printf("Impossible d'ouvrir le fichier data.txt");
	}
int nNode=theProblem->nNode;
for(i = 0; i < nNode; i++)
  {
      for(j = 0; j < nNode; j++)
      {
  	     printf("%lf ",theProblem->Weights[i][j]); //Use lf format specifier, \n is for new line
      }
      printf("\n");
  }


int nEdge=4;
theProblem->nNode=nNode;
theProblem->nEdge=nEdge;
theProblem->nodes=new Node[nNode];
theProblem->edges=new Edge[nEdge];

//NODE0
theProblem->nodes[0].indice=0;
theProblem->nodes[0].degree=1;
theProblem->nodes[0].voisins=new Node*[1];
theProblem->nodes[0].voisins[0]=&theProblem->nodes[1];
//NODE1
theProblem->nodes[1].indice=1;
theProblem->nodes[1].degree=3;
theProblem->nodes[1].voisins=new Node*[3];
theProblem->nodes[1].voisins[0]=&theProblem->nodes[0];
theProblem->nodes[1].voisins[1]=&theProblem->nodes[2];
theProblem->nodes[1].voisins[2]=&theProblem->nodes[3];
//NODE2
theProblem->nodes[2].indice=2;
theProblem->nodes[2].degree=2;
theProblem->nodes[2].voisins=new Node*[2];// degree*sizeof(Node)
theProblem->nodes[2].voisins[0]=&theProblem->nodes[1];
theProblem->nodes[2].voisins[1]=&theProblem->nodes[3];
//NODE3
theProblem->nodes[3].indice=3;
theProblem->nodes[3].degree=2;
theProblem->nodes[3].voisins=new Node*[2];// degree*sizeof(Node)
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

//////////// Incidentes Edges /////////////////

theProblem->nodes[0].incidentes= new Edge*[1];
theProblem->nodes[1].incidentes= new Edge*[3];
theProblem->nodes[2].incidentes= new Edge*[2];
theProblem->nodes[3].incidentes= new Edge*[2];

theProblem->nodes[0].incidentes[0]=&theProblem->edges[0];

theProblem->nodes[1].incidentes[0]=&theProblem->edges[0];
theProblem->nodes[1].incidentes[1]=&theProblem->edges[1];
theProblem->nodes[1].incidentes[2]=&theProblem->edges[2];

theProblem->nodes[2].incidentes[0]=&theProblem->edges[1];
theProblem->nodes[2].incidentes[1]=&theProblem->edges[3];

theProblem->nodes[3].incidentes[0]=&theProblem->edges[2];
theProblem->nodes[3].incidentes[1]=&theProblem->edges[3];

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
//    delete x;
//    x = new int;
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
Function findCycle
Returns the Chemin in the tree that goes from both ends of the edgeCurrent
@pre: A pointer to an off-tree Edge
@post: A pointer to a Chemin structure that goes from the Node edgeCurrent->b to the Node edgeCurrent->a 
*/
Chemin* findCycle(Edge *edgeCurrent,Problem *theProblem)
{
	Chemin *path= new Chemin[1];
	path->theChemin = new Edge*[theProblem->nNode-1]; // faudra en enlever à la fin car il est trop long 
	
	int numberEdges = 0;
	Node *nodeA= edgeCurrent->a;
	Node *nodeB= edgeCurrent->b;
	int nextB=0; 
	int currentB=nodeB->indice;
	int currentA=nodeA->indice;
	int *tabA=new int[theProblem->nNode-1];// tableau d'indice des nodes sur le chemin
	int sizeB=0; 
	int sizeA=0; 
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
		if(stopNode >=0) break; 

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


// ON N'UTILISE PLUS CETTE FONCTION POUR L'INSTANT !
//Chemin* DFS(Edge *edgeCurrent,Problem *theProblem)
//{
///*
//Specification:
//Returns a Chemin structure that goes from the Node "edgeCurrent->b" to "edgeCurrent->a".
//edgeCurrent is an edge that does not belong to theTree.
//*/
//	
//	/// DFS ///
//	int i; 
//	int found=0;
//	int *generators = new int[theProblem->nNode];
//	//int *generators = malloc(sizeof(int)*theProblem->nNode);
//	int numberEdges = 0;
//	Node *nodeA= edgeCurrent->a;
//	Node *nodeB= edgeCurrent->b;
//	/// 1) DFS from nodeA tot nodeB ////
//	queue<int> myqueue;
//	myqueue.push(nodeA->indice);
//	
//	while(!myqueue.empty())
//	{
//        if(myqueue.front() == nodeB->indice)
//        {
//        	found =1;// check pour dire qu'on a trouvé le node destination 
//        	break;
//        } 
//	else
//	{
//		int currentNode=myqueue.front();// on va spliter ce node
//        	myqueue.pop();
//        	for (i = 0; i < theProblem->nodes[currentNode].degree ; i++)
//            	{	
//			myqueue.push(theProblem->nodes[currentNode].voisins[i]->indice);
//			generators[theProblem->nodes[currentNode].voisins[i]->indice]=currentNode;// TODO tc                 	
//            	}
//            }
//    	}
//    
//	// sanity check
//	if (found ==1)	
//	{
//		 cout << "on a found le chemin \n ";
//	}
//	else{
//		cout << "ERROR: chemin not found in DFS \n";
//	}
//	
//	/// 2) Recover path ////
//	
//	// la structure à remplir, (puis changer ordre pour les edges ?)
//	//Chemin *chemin=malloc(sizeof(Chemin));
//	Chemin *chemin=new Chemin[1];
//	chemin->size=0;
//	chemin->theChemin=NULL;
//	int *reverse=new int[theProblem->nEdge];
//	//malloc(theProblem->nEdge*sizeof(Edge));// on le rempli trop grand mais on va pas jusqu'au bout
//	
//	int k=nodeB->indice;						   	
//	while( generators[k] != nodeA->indice)
//	{// remplir reverse en parcourant le mapping
//		
//		// trouver l'indice de l'edge qui relie les nodes d'indices k et generators[k]:
//		for(i=0; i< theProblem->nodes[k].degree; i++ )
//		{
//			if( (theProblem->nodes[k].incidentes[i]->a->indice == k &&
//			     theProblem->nodes[k].incidentes[i]->b->indice == generators[k] ) || 
//			     (theProblem->nodes[k].incidentes[i]->a->indice == generators[k] &&
//			     theProblem->nodes[k].incidentes[i]->b->indice == k ) )
//			{
//				reverse[chemin->size]=theProblem->nodes[k].incidentes[i]->indice;
//				break; 	
//			}
//		}
//		chemin->size= chemin->size+1;	
//		k=generators[k];
//	}
//	chemin->theChemin=new Edge*[chemin->size];
//	//malloc(theChemin->size*sizeof(Edge));
//	for(i=0; i<chemin->size; i++)
//	{
//		chemin->theChemin[i]=&(theProblem->edges[reverse[i]]);// theChemin c'est des pointeurs et reverse c'est des indices
//	}
	
/* ON A PEUT ETRE PAS BESOIN DE RETOURNER LE CYCLE, DEPEND DES SIGNES de flots ETC, TODO à vérifier */	
//	chemin->size=numberEdges;
//	chemin->theChemin=malloc(numberEdges*sizeof(Edge));
//	for(i=0; i<numberEdges; i++)
//	{
//		chemin->theChemin[i]=reverse[numberEdges-1-i];
//	}
	
//	return chemin;
//}



////////////////////////// HAROLD'S PRATICE /////////////////////////////

double stretchEdge(Edge *edgeCurrent, Edge **chemin, int length)
{
    double stretch = 0.0;
    int i; 
    for(i = 0; i < length; i++)
    {
        stretch = stretch + 1/(chemin[i]->weight);
    }
    stretch = stretch*(edgeCurrent->weight);
    return stretch;

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
}

/* A METTRE A JOUR (en fontion de la structure theProblem)
double stretchTree(Problem *theProblem, Tree *arbre)
{
    nEdge = theProblem->nEdge;
    listEdges = theProblem->listEdges; A CHANGER

    double stretch = 0;

    for(i = 0; i < nEdge; i++)
    {
        edgeCurrent = listEdges[i];
        chemin = ...  A CHANGER
        length = ...  A CHANGER

        stretch = stretch + stretchEdge(edgeCurrent, chemin, length);
    }
    return stretch;
}
*/








