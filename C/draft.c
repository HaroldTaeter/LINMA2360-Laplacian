    #include<stdio.h>
    //#include<conio.h>
    #include<stdlib.h>
    int i,j,k,a,b,u,v,n,ne=1;
    int min,mincost=0,cost[9][9],parent[9];
    int find(int);
    int uni(int,int);
    int Wmax=999;
    int *treeIndex= malloc((nNode-1)*sizeof(int));
    int compteur=0;
    void main()
    {
    	 while(ne < nNode)
    	 {
    		for(i=1,min=Wmax;i<=n;i++)
    		{
    			for(j=1;j <= n;j++)
    			{
    				if(cost[i][j] < min)
    				{
    					min=cost[i][j];
    					a=u=i;
    					b=v=j;
    				}
    			}
    		}
    	u=find(u);
    	v=find(v);
    	if(uni(u,v))
    	{
    		printf("%d edge (%d,%d) =%d\n",ne++,a,b,min);
    		indiceEdge=findIndex(a,b,theProblem);
    		treeIndex[compteur]=indiceEdge;
    		compteur++;
    
    // TODO c'est ici qu'il faut choper l'indice de l'edge
    		mincost +=min;
    	}
    		cost[a][b]=cost[b][a]=Wmax;
    	}
    //getch();
    	}
    		int find(int i)
    		{
    			while(parent[i])
    			i=parent[i];
    			return i;
    		}
    ///////////////////////////////////////////////		
    int uni(int i,int j)
    {
    	if(i!=j)
    	{
    		parent[j]=i;
    		return 1;
    	}
    	return 0;
    }
    
    int findIndex(a,b,theProblem)
    {
    	int i; 
    	for(i=0; i<nEdge; i++)
    	{
    		if( ( theProblem->edges[i].a->indice==a || theProblem->edges[i].b->indice==a) && 
    			( theProblem->edges[i].a->indice==b || theProblem->edges[i].b->indice==b) )	
    		{	
    			return theProblem->edges[i].indice;
    		}
    	}
    }
