#ifndef GRAPH_H
#define GRAPH_H

#include <map>
#include <string>
#include <vector>
#include <queue>

/* Prevents bypassing degree 2 vertex in U that has both neighbors in N
 * and the neighbors are nonadjacent.
 * Crucial for all algorithms that use tents in their analysis. */
//#define DO_NOT_BYPASS_DEG2_U_VERTICES_WITH_NEIGHBORS_IN_N

using std::vector;
using std::string;
using std::map;
using std::pair;
using std::queue;

const int STATUS_IN_SOLUTION = -1;
const int STATUS_DELETED = -2;
const int STATUS_NONE = 0;
const int STATUS_IN_SAFE = 1;

class Graph {
   private:
    int N, M;  // number of vertices and edges
    int deleted;
    bool changeInU;
    vector<vector<int>> neighbours;
    // Information whether the edge is a double edge
    vector<vector<bool>> neighDouble;

    // Number of double edges of given vertex
    vector<int> doubleNum;

    // Keeps the position of my vertex at the neighbour list
    // For faster deleting neighbours
    vector<vector<int>> posAtNeigh;

    // State of the given vertex
    vector<int> status;

    // Represents vertices numbers before cutting out induced subgraph
    vector<int> oldNumbers;

    // X - Vertices in solution
    vector<int> X;

    // U - Already safe vertices
    // vector<int> U;
 
    // Vertices with selfloops - have to be deleted
    vector<int> selfLoops;

    queue<int> deg2ToReduce;

    vector<int> vis;
    vector<int> ccNr;
    vector<vector<int>> ccVertices;

    // For find union
    vector<int> fuRepr;
    vector<int> fuWeight;

   public:
    Graph();

    // Build induced subgraph of full using vertices from ccNr
    Graph(Graph& full, int ccNr);
    void buildInducedSubgraph(Graph& full, int ccNr);
    void buildGraph(int N_, int M_, vector<pair<int, int> >& edges);
    void printGraph();

    int getN();

    void debugCheckReduced();
    void debugCheckDeg3Applicable();
    void populateDeg2Queue();

    // void add_vertex(int v);
    void pureDeleteVertex(int v);
    void deleteVertex(int v);
    void deleteNeighbour(int v, int which);
    int addEdge(int v, int w, bool isDouble);
    void removeEdge(int v, int w);
    //Merges vertex v into s (puts its neighbour in it)
    void merge(int s, int v);
    vector<int> getDoubleNeighbours(int v);

    int degree(int v);

    //Degree of a vertex including it double edges
    int degWithDouble(int v);
    bool isEmpty();

    // Returns pair <v_num, degree>, where v is not deleted nor safe vertex with max deg
    pair<int, int> getMaxDegV();

    // Same as above, but the vertex is neighbour of s
    pair<int, int> getMaxDegNeigh(int s);

    // Same as above, just takes under consideration also double edges
    pair<int, int> getMaxDegVWDouble();
    
    // Same as above, just takes under consideration just edges to undeletable part
    pair<int, int> getMaxDegToUndeletable();

    // Same as above, but sorts by (max double degree, max degree)
    pair<int, int> getMaxDegPreferDouble();

    bool connected(int v, int to);
    bool hasSelfLoop(int v);

    bool isDoubleEdge(int v, int to);
    bool hasDoubleToU(int v);

    bool deleteSelfLoops();
    //Deletes Vertices with deg <= 1
    bool deleteDeg1Vertices();
    void deleteDeg1Chain(int v);

    int doubleNr(int v);
    int doubleNr2(int v);

    int findNrInNeigh(int where, int me);
    void mergeDeg2Neigh(int cur);

    // Reduce deg2 Vertices from the queue
    bool reduceDeg2Queue();

    // Reduction from Kociumaka-Pilipczuk
    // For a deg3-vertex in Nwith 2 neighbors in U and one in N,
    // subdivide the edge to the vertex in N and put the new vertex in U.
    bool reduceTents();

    void putInSolution(int v);
    void purePutInSolution(int v);
    void putInSafe(int v);

    vector<int> getX();
    vector<int>& getStatuses();
    vector<vector<int>>& getNeighbours();

    int getXSize();

    bool isNone(int v);
    bool isInU(int v);
    bool isInX(int v);
    bool isDeleted(int v);
    bool checkStatus(int v, int status);
    int getStatus(int v);
    void setStatus(int v, int st);

    bool checkCycles();
    bool dfsCycles(int v, int from);

// Checks if some vertex has two neighbours in same CC in U (adding it to U would create a cycle)
// checks only for vertices in None and their neighbours in U
    bool checkSameCCNeighboursOneVertex(int v);
    bool checkSameCCNeighbours();


// Recalculates components, returns nr of CC with most vertices
    int recalcComponents();
    void dfsComp(int v, int nr);


    vector<pair<int, int>> createAllEdgeList();
// Creates a vector of all edges in the CC
    vector<pair<int, int>> createCCEdgeList(int nr);
    vector<int> getCC(int nr);

    int getNumberOfCCs();
    int getCCSize(int nr);
    int calcNotNoneInCC(int nr);
    bool isAnyNoneInCC(int nr);

// Provides number in the graph before cutting out the CC
    int getOriginalNr(int nr);
    void printOriginalNums();

// Debug print of components
    void printCC(int nr);
    void printCCs();

    int fFind(int x);
    void fUnion(int x, int y);
    pair<int, int> getFAUVals(int x);

    // Makes union with all neighbours of a vertex that belong to U
    void mergeWithNeighinU(int v);

    int allDoubleCnt();
    bool reduceTooManyDouble(int bound);

    // Used for tents reduction
    bool isPotentialTent(int v);

    // Simple degree-counting lower bound
    int degLowerBound();
};

#endif  // GRAPH_H
