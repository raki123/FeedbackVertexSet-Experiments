#ifndef HALFINT_H
#define HALFINT_H

#include <map>
#include <vector>
#include "Graph.h"

using std::map;
using std::pair;
using std::vector;

// O - free, I - on path, N - on leg , H - on cycle
enum class eWhat {O, I, N, H};

// T - cycle + end of leg
enum class vType {O, I, N, H, T, s};

const int NOEDGE = -1;
const double EPS = 1e-6;

class HalfIntRelax {
public:
	HalfIntRelax(Graph& g, int s);

	void findPacking();
	void transformToCover();
	void findCover();
	void coverDfs(int v, int edge, bool& was);
	vector<double> getCover();
	double getCoverSum();
    
    // Returns vertices from X^-1(1)
	vector<int> getX1();

	// Returns vertices reachable through X^-1(0) vertices
	vector<int> getReachableThrough0();

    // Returns bridges between s and tree CCs in G - X - s
	vector<int> getBridgesToTreeCC();
	bool isTreeCC(int v, int edge);

	void update();

	int getType(int v);
	bool isType(int v, vType t);
	bool findAugmentingPath(int v, int edge);

	void update_hit_flower();
	void update_hit_s();
	void update_hit_self();
	void zeroPrev();

	bool isLegEnd(int posAtStack);
	void clearUnreachableI();

	int findAndChange(int v, vType target, vType newT, eWhat oldE, eWhat newE);
	int dfsTypeChange(int v, int edge, vType target, vType newT, eWhat oldE, eWhat newE);
	bool reachesSThroughI(int v);
	bool reachSDfs(int v, int edge);

	int clearLeg(int v, int edge); 
	void dfsClearLeg(int v, int edge, vector<pair<int, int>>& st, int& legEnd, int start); 
	int countLegs(int v);
	void decomposeFlower(int v, int edge, bool includeBeg);

	vector<pair<int, int>> getCycle(int v, int edge);
	void getCycleDFS(int v, int edge, vector<pair<int, int>>& res, int start);
	int findNextLeg(int pos, vector<pair<int, int>>& cyc) const;
	void connectWithI(int curLeg, int nextLeg, vector<pair<int, int>>& cycle);
	void connectWithO(int curLeg, int nextLeg, vector<pair<int, int>>& cycle);
	int getEdgeNr(int v1, int v2);

	void xorStack(int beg, int end);

	void prepareNewVisit();
	bool visited(int v);
	void visit(int v);

// Debug functions:
	vType brutVType(int v);
	bool brutCheckAllTypes();
	bool checkCoverDfs(int v, int edge);
	bool checkCover();
    void printStackTrace();
	void printVTypes();
	void printVType(int i);
	void printETypes();
	void printETypes(int i);
	void printEType(int eNr);
	void printCover();
	void printAllDebug();
	void printPrev();

private: 
    void buildLocalGraph(Graph& g);
	int N, M;
	int s;
    // prev - edge from which we came to a given vertex
	vector<int> prev;
    vector<pair<int, int>> stos;
    vector<vector<int>> neighbours;
    vector<vector<int>> numberE;
	vector<vType> type;
	vector<eWhat> what;

	int curVis;
	vector<int> visi;
	vector<double> cover;
	double coverSum;
};

#endif