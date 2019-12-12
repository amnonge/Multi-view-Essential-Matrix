#include "mex.h"
#define NIL -1
//Part of the code was adopted from here: ttps://www.geeksforgeeks.org/articulation-points-or-cut-vertices-in-a-graph/
int minn(int a, int b)
{
    if (a<b){
        return a;
        
    }else{
        return b;
    }
}
void removeNode(int n, int V, int * numedges, int ** adj){
    int i; int foundn; int * temp;
    for (i = 0; i < V; i++){
        if (i != n){
            //remove it from the list
            //			int** cur = &adj[i];
            foundn = -1;
            
            for (int j = 0; j < numedges[i]; j++){
                if (adj[i][j] == n){
                    foundn = j;
                }
                else{
                    if (adj[i][j]  > n){
                        adj[i][j] = adj[i][j] - 1;
                    }
                    
                }
            }
            if (foundn >= 0){
                for (int j = foundn; j < numedges[i] - 1; j++){
                    adj[i][j] = adj[i][j + 1];
                }
                numedges[i] = numedges[i] - 1;
            }
            
            
        }
    }
    temp = adj[n];
    
    for (int i = n; i < V - 1; i++){
        adj[i] = adj[i + 1];
        numedges[i] = numedges[i + 1];
    }
    adj[V - 1] = temp;
    
}

void initiGraph(int V, int ** numedges, int *** adj)
{
    int i;
    *adj=(int**)malloc(sizeof(int*)*V);
    *numedges = (int*)malloc(sizeof(int)*V);
    for ( i = 0; i < V; ++i){
        (*adj)[i] = (int*)malloc(sizeof(int)*V);
        (*numedges)[i] = 0;
    }
    
}

void freeGraph(int origVertices, int * numedges, int ** adj)
{
    int i;
    
    for ( i = 0; i < origVertices; ++i){
        free(adj[i]);
        
    }
    free(numedges);
    free(adj);
    
}

void addEdge(int v, int w, int * numedges, int ** adj)
{
    
    adj[v][numedges[v]] = w;
    
    numedges[v] = numedges[v] + 1;
    
    adj[w][numedges[w]] = v;
    numedges[w] = numedges[w] + 1;
    
}

// A recursive function that find articulation points using DFS traversal
// u --> The vertex to be visited next
// visited[] --> keeps tract of visited vertices
// disc[] --> Stores discovery times of visited vertices
// parent[] --> Stores parent vertices in DFS tree
// ap[] --> Store articulation points
void APUtil(int u, bool visited[], int disc[],
        int low[], int parent[], bool ap[], int * numedges, int ** adj)
{
    // A static variable is used for simplicity, we can avoid use of static
    // variable by passing a pointer.
    int children; int i;
    static int time = 0;
    
    // Count of children in DFS Tree
    children = 0;
    
    // Mark the current node as visited
    visited[u] = true;
    
    // Initialize discovery time and low value
    disc[u] = low[u] = ++time;
    
    // Go through all vertices aadjacent to this
    //list<int>::iterator i;
    //for (i = adj[u].begin(); i != adj[u].end(); ++i)
    for ( i = 0; i < numedges[u]; i++)
    {
        int v = adj[u][i];  // v is current adjacent of u
        
        // If v is not visited yet, then make it a child of u
        // in DFS tree and recur for it
        if (!visited[v])
        {
            children++;
            parent[v] = u;
            APUtil(v, visited, disc, low, parent, ap, numedges, adj);
            
            // Check if the subtree rooted with v has a connection to
            // one of the ancestors of u
            low[u] = min(low[u], low[v]);
            
            // u is an articulation point in following cases
            
            // (1) u is root of DFS tree and has two or more chilren.
            if (parent[u] == NIL && children > 1)
                ap[u] = true;
            
            // (2) If u is not root and low value of one of its child is more
            // than discovery value of u.
            if (parent[u] != NIL && low[v] >= disc[u])
                ap[u] = true;
        }
        
        // Update low value of u for parent function calls.
        else if (v != parent[u])
            low[u] = min(low[u], disc[v]);
    }
}

// The function to do DFS traversal. It uses recursive function APUtil()
void AP(bool *ap, bool *visited, int *disc, int *low, int *parent, int * numedges, int ** adj,int V)
{
    int i;
    
    for ( i = 0; i < V; i++)
    {
        parent[i] = NIL;
        visited[i] = false;
        ap[i] = false;
    }
    
    // Call the recursive helper function to find articulation points
    // in DFS tree rooted with vertex 'i'
    for ( i = 0; i < V; i++)
        if (visited[i] == false)
            APUtil(i, visited, disc, low, parent, ap, numedges, adj);
    
    
}
bool checkAllCameras(int * tr, bool * mask, int camNum, int tripletNum, int curindTripletCanRemove, int * camInds)
{
    int i;
    for ( i = 0; i < camNum; i++){
        camInds[i] = 0;
    }
    for (i = 0; i < tripletNum; i++){
        if (mask[i] && i != curindTripletCanRemove){
            int ind1 = tr[3 * i];
            int ind2 = tr[3 * i + 1];
            int ind3 = tr[3 * i + 2];
            camInds[ind1 - 1] = 1;
            camInds[ind2 - 1] = 1;
            camInds[ind3 - 1] = 1;
            
            
            
        }
    }
    
    for (i = 0; i < camNum; i++){
        if (camInds[i] == 0){
            return false;
        }
    }
    return true;
    
}


void getAB(bool * mask,int * idds,int * tr,int * allCams,size_t numTriplets,size_t  numcam,int *mapTriangleIndices,int * edges,size_t numedgess,int *  mapCurGraphOrig,bool * mask_cutting ,bool * maskcc){
    
    int i; bool foundsomething; bool *ap; int *disc; int *low; int curind; bool condition1; int jj; int *parent;
    
    
    int origVertices;
    int **adj;    // A dynamic array of adjacency lists
    int *numedges;
    int  V; bool *visited; int temptry;
   
    V = numTriplets;
    initiGraph(V,&numedges, &adj);
    
    
    
    for ( i = 0; i < numedgess; i++){
        addEdge(edges[2 * i] - 1, edges[2 * i + 1] - 1, numedges, adj);
        
    }
       
    foundsomething = true;
    ap = (bool*)malloc(sizeof(bool)*numTriplets);
    
    visited = (bool*)malloc(sizeof(bool)*numTriplets);
    disc = (int*)malloc(sizeof(int)*numTriplets);
    low = (int*)malloc(sizeof(int)*numTriplets);
    parent = (int*)malloc(sizeof(int)*numTriplets);
    while (foundsomething == true)
    {
        AP(ap, visited, disc, low, parent, numedges, adj,V);
        for ( i = 0; i < numTriplets; i++){
            mask_cutting[i] = 1;
        }
        for ( i = 0; i < V; i++){
            mask_cutting[mapCurGraphOrig[i]] = !ap[i];
        }
        
        foundsomething = false;
        for ( i = 0; i<numTriplets; i++){
            curind = idds[i] - 1;
            if (mask[curind] && maskcc[curind] && mask_cutting[curind]){
                
                condition1 = checkAllCameras(tr, mask, numcam, numTriplets, curind, allCams);
                if (condition1){
                    removeNode(mapTriangleIndices[curind],V,numedges,adj);
                    V--;
                    mask[curind] = 0;
                    
                    
                    temptry = mapTriangleIndices[curind];
                    //we need to have a map to the current indices of each non removed triplet in the graph
                    for ( jj = curind; jj < numTriplets; jj++){
                        mapTriangleIndices[jj] = mapTriangleIndices[jj] - 1;
                    }
                    for ( jj = temptry; jj < numTriplets - 1; jj++){
                        mapCurGraphOrig[jj] = mapCurGraphOrig[jj + 1];
                    }
                    //std::cout << curind << std::endl;
                    foundsomething = true;
                    //update the backward inbds!!!!!!!!!!
                    i = numTriplets;
                    
                    
                    
                }
                else{
                    //WE can not remove this triplet because this is the only triplet that contain some camera.
                    //As a result we will never check it again
                    maskcc[curind] = 0;
                }
            }
        }
        
    }
    
    freeGraph(numTriplets, numedges, adj); 
    free(ap);
    free(visited);
    free(disc);
    free(low);
    free(parent);
    
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    int *idds;bool * mask;int * tr;int * allCams;int * mapTriangleIndices;
    int * edges;int *mapCurGraphOrig;bool * maskcc; bool *mask_cutting;
    
    double* output;
    
    size_t numTriplets,numCams,numEdges;
    
    
    
    idds = (int *) mxGetData(prhs[0]);;
    mask = (bool *) mxGetData(prhs[8]);;
    mapCurGraphOrig = (int *) mxGetData(prhs[9]);;
    maskcc = (bool *) mxGetData(prhs[10]);;
    mask_cutting = (bool *) mxGetData(prhs[11]);;
    //mask= mxGetPr(prhs[]);
    tr= (int *) mxGetData(prhs[1]);;
    allCams= (int *) mxGetData(prhs[2]);;
    numTriplets = (int)mxGetScalar(prhs[5]);
    numCams =(int) mxGetScalar(prhs[6]);
    numEdges =(int) mxGetScalar(prhs[7]);
    mapTriangleIndices= (int *) mxGetData(prhs[3]);;
    edges= (int *) mxGetData(prhs[4]);;
    
    
    /*  get the dimensions of the matrix input y */
    
    /*  set the output pointer to the output matrix */
    plhs[0] = mxCreateDoubleMatrix( (mwSize)numTriplets, (mwSize)1, mxREAL);
    
    /*  create a C pointer to a copy of the output matrix */
    
   
    /*  call the C subroutine */
    getAB(mask,idds,tr,allCams,numTriplets,numCams,mapTriangleIndices,edges,numEdges,mapCurGraphOrig,maskcc,mask_cutting);
    output = mxGetPr(plhs[0]);
    
    for(int i=0;i<numTriplets;i++)
    {
        output[i]=mask[i]*2;
    }
}
