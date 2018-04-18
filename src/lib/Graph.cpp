#include "Graph.h"

#include <cassert>
#include <set>
#include <iostream>
#include <algorithm>

using namespace std;

Graph::Graph(): deleted(0), changeInU(true) {}


Graph::Graph(Graph& full, int ccNr) : deleted(0), changeInU(true) {
    buildInducedSubgraph(full, ccNr);
}

void Graph::buildInducedSubgraph(Graph& full, int ccNr) {
    vector<pair<int, int>> E = full.createCCEdgeList(ccNr);
    int n = full.getCCSize(ccNr);
    int m = E.size();

    vector<int> oldNums = full.getCC(ccNr);
    
    map<int, int> oldToNew;

    for(int i = 0; i < oldNums.size(); i++) {
        oldToNew[oldNums[i]] = i;
        assert(full.isNone(oldNums[i]) || full.isInU(oldNums[i]));
    }

    // Renumber edges
    for(int i = 0 ; i < E.size(); i++)  {
        int a = oldToNew[E[i].first];
        int b = oldToNew[E[i].second];

        E[i].first = a;
        E[i].second = b;
    }
    // Build subgraph from edges in the given CC
    buildGraph(n, m, E);

    // Fix oldnumbers values
    oldNumbers = oldNums;

    // putInSafe should recompute Find&Union correctly
    for(int i = 0; i < oldNums.size(); i++) {
        if(full.checkStatus(oldNums[i], STATUS_IN_SAFE))
            putInSafe(i);
    }

    // Fill reduce queue 
    for(int i = 0; i < N; i++) {
         if (degree(i) == 2 && (isInU(i) || isNone(i)))
            deg2ToReduce.push(i);
    }
}

void Graph::buildGraph(int N_, int M_, vector<pair<int, int> >& edges) {
    vector<map<int, int> > neighMap;
    N = N_;
    M = M_;

    neighbours.resize(N);
    status.resize(N, STATUS_NONE);
    neighDouble.resize(N);
    neighMap.resize(N);
    posAtNeigh.resize(N);
    doubleNum.resize(N, 0);

    fuRepr.resize(N, -1);
    fuWeight.resize(N, 0);

    oldNumbers.resize(N);

    for(int i = 0; i < N; i++)
        oldNumbers[i] = i;       

    for (int i = 0; i < edges.size(); i++) {
        int a = edges[i].first;
        int b = edges[i].second;
        if (a > b)
            swap(a, b);

        if (a == b)
            selfLoops.push_back(a);

        if (!neighMap[a].count(b)) {
            int pos = neighbours[a].size();
            posAtNeigh[b].push_back(pos);
            neighMap[a][b] = pos;
            neighbours[a].push_back(b);
            neighDouble[a].push_back(false);
            
            pos = neighbours[b].size();
            posAtNeigh[a].push_back(pos);
            neighMap[b][a] = pos;
            neighbours[b].push_back(a);
            neighDouble[b].push_back(false);
        }
        else {
            if(!neighDouble[a][neighMap[a][b]])
                doubleNum[a]++;

            if(!neighDouble[b][neighMap[b][a]])
                doubleNum[b]++;
                
            neighDouble[a][neighMap[a][b]] = true;
            neighDouble[b][neighMap[b][a]] = true;
        }
    }
}

void Graph::printGraph() {
    int cnt_n = 0, cnt_u = 0, cnt_d = 0, cnt_x = 0;
    for (int i = 0; i < N; i++) {
        cerr << i << " ";
        switch(getStatus(i)) {
            case STATUS_NONE: cerr << "N "; cnt_n++; break;
            case STATUS_IN_SOLUTION: cerr << "X "; cnt_x++; break;
            case STATUS_DELETED: cerr << "D "; cnt_d++; break;
            case STATUS_IN_SAFE: cerr << "U "; cnt_u++; break;
        }

        cerr << ": ";
        for (int j = 0; j < neighbours[i].size(); j++) {
            cerr << "(" << neighbours[i][j] << ", " << neighDouble[i][j] << ", " << posAtNeigh[i][j] << ") ";
        }
        cerr << "\n";
    }

    cerr << "Czy pusty: " << isEmpty() << " " << !checkCycles() << " ";
    cerr << "Num = (" << neighbours.size() << ", " << cnt_n + cnt_u + cnt_d + cnt_x << ") " <<
            "In N = " << cnt_n << " in X = (" << cnt_x   << ", " << X.size() << ") in D = (" << cnt_d <<
            ", " << deleted << ") in U = (" <<  cnt_u << ")\n";

        for(int i=0; i< X.size() ;i++){
            cerr << X[i] << " ";
        }
        cerr << "=X\n";

    cerr << "Deg2queue " << deg2ToReduce.size() << "\n";
}

int Graph::getN() {
    return N;
}

int Graph::degree(int v) {
    return neighbours[v].size();
}

int Graph::degWithDouble(int v) {
    return degree(v) + doubleNum[v];
}

int Graph::doubleNr2(int v) {
    return doubleNum[v];
}

pair<int, int> Graph::getMaxDegV() {
    int deg = 0, num = -1;
    for (int i = 0; i < neighbours.size(); i++)
        if (deg < degree(i) && checkStatus(i, STATUS_NONE)) {
            deg = degree(i);
            num = i;
        }

    return make_pair(num, deg);
}

pair<int, int> Graph::getMaxDegNeigh(int s) {
    int deg = 0, num = -1;
    for (int i = 0; i < neighbours[s].size(); i++) {
        int nei = neighbours[s][i];
        if (deg < degWithDouble(nei)) {
            deg = degWithDouble(nei);
            num = nei;
        }
    }
    return make_pair(num, deg);
}

pair<int, int> Graph::getMaxDegVWDouble() {
    int deg = 0, num = -1;
    for (int i = 0; i < neighbours.size(); i++)
        if (deg < degWithDouble(i) && checkStatus(i, STATUS_NONE)) {
            deg = degWithDouble(i);
            num = i;
        }

    return make_pair(num, deg);
}

pair<int, int> Graph::getMaxDegToUndeletable() {
    int degundel = -1, deg = -1, num = -1;
    for (int i = 0; i < neighbours.size(); i++) {
        if (!checkStatus(i, STATUS_NONE))
            continue;
        int cnt = 0;
        for (int j = 0; j < neighbours[i].size(); ++j) {
            int u = neighbours[i][j];
            if (checkStatus(u, STATUS_IN_SAFE)) {
                cnt++;
            }
        }
        if (degundel < cnt || (degundel == cnt && deg < degWithDouble(i))) {
            degundel = cnt;
            deg = degWithDouble(i);
            num = i;
        }
    }

    return make_pair(num, degundel);
}

pair<int, int> Graph::getMaxDegPreferDouble() {
    int ddeg = -1, deg = -1, num = -1;
    for (int i = 0; i < neighbours.size(); i++)
        if (checkStatus(i, STATUS_NONE) && (doubleNr2(i) > ddeg || doubleNr2(i) == ddeg && degWithDouble(i) > deg)) {
            deg = degWithDouble(i);
            ddeg = doubleNr2(i);
            num = i;
        }

    return make_pair(num, deg);
}


bool Graph::connected(int v, int to) {
    assert(v >= 0 && v < neighbours.size());
    for (int i = 0; i < neighbours[v].size(); i++)
        if (neighbours[v][i] == to)
            return true;

    return false;
}

bool Graph::hasSelfLoop(int v) {
    return connected(v, v);
}

bool Graph::deleteSelfLoops() {
    bool res = false;
    for(int i = 0; i < selfLoops.size(); i++) {
        int cur = selfLoops[i];
        if(isNone(cur)) {
            putInSolution(cur);
            res = true;
        }
    }

    selfLoops.clear();
    return res;
}

bool Graph::deleteDeg1Vertices() {
    bool res = false;
    for (int i = 0; i < neighbours.size(); i++)
        if (degree(i) <= 1 && (isNone(i) || checkStatus(i, STATUS_IN_SAFE))) {
            deleteDeg1Chain(i);
            res = true;
        }

    return res;
}

void Graph::deleteDeg1Chain(int v) {
    assert(isNone(v) || checkStatus(v, STATUS_IN_SAFE));
    if(degree(v) > 1)
        return;

    int neigh = -1;
    if (degree(v) > 0){
        neigh = neighbours[v][0];
        if(isDoubleEdge(v, 0)) {
            //is deg <= 1 and has double edge -> we delete its neighbour
            if (isNone(neigh)) {
                putInSolution(neigh);
                neigh = -1;
            }
            else {
                assert(isNone(v));
                putInSolution(v);
            }
        }
    }

    if (isNone(v) || checkStatus(v, STATUS_IN_SAFE)) {// Could have got deleted by chaining after deleting double neighbour
        deleteVertex(v);
        setStatus(v, STATUS_DELETED);
        deleted++;
    }
    
    if(neigh >= 0 && (isNone(neigh) || checkStatus(neigh, STATUS_IN_SAFE))) {
        if (degree(neigh) <= 1)
            deleteDeg1Chain(neigh);
        else if (degree(neigh) == 2)
            deg2ToReduce.push(neigh);
    }
}

// I will use it only in small cases (deg 2 vertices), no need to make 
// additional vector to keep it
int Graph::doubleNr(int v) {
    int res = 0;
    for(int i = 0; i < neighbours[v].size(); i++) {
        res += neighDouble[v][i];
    }
    return res;
}

int Graph::findNrInNeigh(int u, int me) {
    for(int i = 0; i < neighbours[u].size(); i++) {
        if (neighbours[u][i] == me) {
            return i;
        }
    }
    return -1;
}

void Graph::mergeDeg2Neigh(int cur) {
    assert(degree(cur) == 2);
    int u = neighbours[cur][0], w = neighbours[cur][1];

    if(isInU(u) && isInU(w))
        fUnion(u, w);

    // We create edge now by putting w/v in u and u/v in w
    int cur_in_u = posAtNeigh[cur][0];
    int cur_in_w = posAtNeigh[cur][1];

    neighbours[u][cur_in_u] = w;
    neighbours[w][cur_in_w] = u;

    // doubleNum is ok, there was no u-w edge

    neighDouble[u][cur_in_u] = false;
    neighDouble[w][cur_in_w] = false;

    posAtNeigh[u][cur_in_u] = cur_in_w;
    posAtNeigh[w][cur_in_w] = cur_in_u;

    pureDeleteVertex(cur);    
}

bool Graph::reduceDeg2Queue() {
    bool res = false;
    while(!deg2ToReduce.empty()) {
        int cur = deg2ToReduce.front();
        deg2ToReduce.pop();

        if(!isNone(cur) && !isInU(cur))
            continue;

        if (degree(cur) > 2) // may happen in LP
            continue;

        if(degree(cur) < 2) {
            res = true;
            deleteDeg1Chain(cur);
            continue;
        }


        int dNum = doubleNr(cur);
        // They both exist
        int u = neighbours[cur][0], w = neighbours[cur][1];

        if (dNum == 1) {
            // u is connected with a double edge
            if (neighDouble[cur][0] == false)
                swap(u, w);

            if (isNone(u)) {
                putInSolution(u);
            } else {
                assert(isNone(cur));
                putInSolution(cur);
            }
            res = true;
        }

        if (dNum == 0) {
            int nr_w_in_u = findNrInNeigh(u, w);

#ifdef DO_NOT_BYPASS_DEG2_U_VERTICES_WITH_NEIGHBORS_IN_N
            if (isInU(cur) && isNone(u) && isNone(w) && nr_w_in_u == -1)
                continue;
#endif

            res = true;
            // Check if they are connected
            if (isInU(u) && isInU(w) && isNone(cur) && fFind(u) == fFind(w)) {
                putInSolution(cur);
                continue;
            }

            // There was edge u -> w
            // Deleting v makes it double
            if(nr_w_in_u != -1) {
                assert(!(isInU(u) && isInU(w)));

                if (isInU(u) || isInU(w))
                    changeInU = true;

                res = true;

                if (!neighDouble[u][nr_w_in_u])
                    doubleNum[u]++;

                if (!neighDouble[w][posAtNeigh[u][nr_w_in_u]])
                    doubleNum[w]++;

                neighDouble[u][nr_w_in_u] = true;
                neighDouble[w][posAtNeigh[u][nr_w_in_u]] = true;

                //REMOVE him from neighbours too
                deleteNeighbour(u, posAtNeigh[cur][0]);
                deleteNeighbour(w, posAtNeigh[cur][1]);
                pureDeleteVertex(cur);
                if (degree(u) <= 2)
                    deg2ToReduce.push(u);
                if (degree(w) <= 2)
                    deg2ToReduce.push(w);
            }
            else { // There was no edge u->w
                // We create it now by putting w/v in u and u/v in w
                if (isInU(u) || isInU(w))
                    changeInU = true;

                mergeDeg2Neigh(cur);
            }
        }
    }
    if (res)
        checkSameCCNeighbours();
    return res;
}

void Graph::putInSolution(int v) {
    assert(isNone(v));
    setStatus(v, STATUS_IN_SOLUTION);
    deleteVertex(v);
    X.push_back(v);
}

void Graph::purePutInSolution(int v) {
    assert(isNone(v));
    assert(v >= 0 && v < neighbours.size());
    setStatus(v, STATUS_IN_SOLUTION);
    X.push_back(v);
    for(int i = 0; i < neighbours[v].size(); i++)
        deleteNeighbour(neighbours[v][i], posAtNeigh[v][i]);
    neighbours[v].clear();
    neighDouble[v].clear();
    posAtNeigh[v].clear();
    doubleNum[v] = 0;
}

//Maybe add merging vertices here?
void Graph::putInSafe(int v) {
    assert(isNone(v));
    changeInU = true;
    setStatus(v, STATUS_IN_SAFE);
    mergeWithNeighinU(v);
}

vector<int> Graph::getX() {
    return X;
}
vector<int>& Graph::getStatuses() {
    return status;
}
vector<vector<int>>& Graph::getNeighbours() {
    return neighbours;
}

int Graph::getXSize() {
    return X.size();
}

bool Graph::isInU(int v) {
    return checkStatus(v, STATUS_IN_SAFE);
}

bool Graph::isInX(int v) {
    return checkStatus(v, STATUS_IN_SOLUTION);
}

bool Graph::isDeleted(int v) {
    return checkStatus(v, STATUS_DELETED);
}

bool Graph::isNone(int v) {
    return checkStatus(v, STATUS_NONE);
}

bool Graph::checkStatus(int v, int status) {
    return getStatus(v) == status;
}

int Graph::getStatus(int v) {
    assert(v >= 0 && v < status.size());
    return status[v];
}

void Graph::setStatus(int v, int st) {
    assert(v >= 0 && v < status.size());
    status[v] = st;
}

//  X + DEL = neighbours.size()?
// U should get reduced completely if there is no N
bool Graph::isEmpty() {
    return (X.size() + deleted) == neighbours.size();
}

bool Graph::hasDoubleToU(int v) {
    for(int i = 0; i < neighbours[v].size(); i++) {
        int cur = neighbours[v][i];
        if(isInU(cur) && neighDouble[v][i])
            return true;
    }
    return false;
}

// Deletes vertex without affecting its neighbours
void Graph::pureDeleteVertex(int v) {
    assert(isNone(v) || isInU(v));
    assert(v >= 0 && v < neighbours.size());
    setStatus(v, STATUS_DELETED);
    deleted++;
    neighbours[v].clear();
    neighDouble[v].clear();
    posAtNeigh[v].clear();
    doubleNum[v] = 0;
}

// Status of vertex is changed by caller
void Graph::deleteVertex(int v) {
    assert(v >= 0 && v < neighbours.size());

    vector<int> d1Neighs;

    for (int i = 0; i < neighbours[v].size(); i++) {
        int cur_neigh = neighbours[v][i];

        // Deletes me at neighbour vectors
        deleteNeighbour(cur_neigh, posAtNeigh[v][i]);

        if (neighbours[cur_neigh].size() <= 1)
            d1Neighs.push_back(cur_neigh);

        if (neighbours[cur_neigh].size() == 2)
            deg2ToReduce.push(cur_neigh);
    }

    neighbours[v].clear();
    neighDouble[v].clear();
    posAtNeigh[v].clear();
    doubleNum[v] = 0;

    for (int i = 0; i < d1Neighs.size(); i++)
        if (isNone(d1Neighs[i])) // TODO
            deleteDeg1Chain(d1Neighs[i]);
}

void Graph::deleteNeighbour(int v, int which) {
    // Self loop -> don't mess up the vectors 
    if (v == neighbours[v][which])
        return;

    assert(v >= 0 && v < neighbours.size());
    assert(which >= 0 && which < neighbours[v].size());

    int del_neigh = neighbours[v][which];
    int last_neigh = neighbours[v].back();

    neighbours[v][which] = neighbours[v].back();
    neighbours[v].pop_back();


// Update doubleNum if deleted neighbour was double
    if(v != del_neigh && neighDouble[v][which]) {
        doubleNum[v]--;
        if(doubleNum[v] < 0) {
            cerr << v << " " << status[v]<< " "<< doubleNum[v] <<  " " << del_neigh << "\n";
            assert(false);

        }
    }

    neighDouble[v][which] = neighDouble[v].back();
    neighDouble[v].pop_back();

    posAtNeigh[v][which] = posAtNeigh[v].back();
    posAtNeigh[v].pop_back();

    if (del_neigh != last_neigh) {
        //we need to set pos at the neighbour correctly
        int myPos = posAtNeigh[v][which];
        posAtNeigh[last_neigh][myPos] = which;
    }   
}

int Graph::addEdge(int v, int w, bool isDouble) {
    int res = 0;
    if(neighbours[v].size() > neighbours[w].size())
        swap(v, w);

    for(int i = 0; i < neighbours[v].size(); i++) {
        if(neighbours[v][i] == w) {
            if(!neighDouble[v][i]) {
                doubleNum[v]++;
                doubleNum[w]++;
                res += 1;
            }
            neighDouble[v][i] = true;
            neighDouble[w][posAtNeigh[v][i]] = true;
            return res;
        }
    }

    int pos = neighbours[v].size();
    neighbours[v].push_back(w);
    posAtNeigh[w].push_back(pos);

    pos = neighbours[w].size();
    neighbours[w].push_back(v);
    posAtNeigh[v].push_back(pos);

    neighDouble[v].push_back(isDouble);
    neighDouble[w].push_back(isDouble);
     
    res = 1 + isDouble;

    if (isDouble) {
        doubleNum[v]++;
        doubleNum[w]++;  
    }

    if (isInU(v) && isInU(w))
        fUnion(v, w);

    return res;
}

void Graph::removeEdge(int v, int w) {
    int pos_in_v = findNrInNeigh(v, w), pos_in_w = findNrInNeigh(w, v);

    assert(pos_in_v != -1);
    assert(pos_in_w != -1);
    deleteNeighbour(v, pos_in_v);
    deleteNeighbour(w, pos_in_w);
}

void Graph::merge(int s, int v) {
    assert(isInU(s));
    vector<int> neighs = neighbours[v];
    vector<bool> dou = neighDouble[v];
    // Delete v from their neighbours
    for (int i = 0; i < neighbours[v].size(); i++) {
        int cur_neigh = neighbours[v][i];
        deleteNeighbour(cur_neigh, posAtNeigh[v][i]);
    }

    neighbours[v].clear();
    neighDouble[v].clear();
    posAtNeigh[v].clear();
    doubleNum[v] = 0;
    
    if (!isInU(v) && !isDeleted(v))
        putInSafe(v);

    // Add neighbours of v to s
    for(int i = 0; i < neighs.size(); i++)
        if(neighs[i] != s)
            addEdge(s, neighs[i], dou[i]);

    // s is in U -> double neighbours have to be put in solution
    vector<int> doubleFromS = getDoubleNeighbours(s);
    for(int i = 0; i < doubleFromS.size(); i++)
        if(!isInX(doubleFromS[i]))
            putInSolution(doubleFromS[i]);
}

vector<int> Graph::getDoubleNeighbours(int v) {
    vector<int> res;

    for(int i = 0; i < neighbours[v].size(); i++)
        if(neighDouble[v][i])
            res.push_back(neighbours[v][i]);

    return res;
} 

bool Graph::checkCycles() {
    vis.clear();
    vis.resize(N, 0);
    bool res = false;
    for (int i = 0; i < N; i++) {
        if(!vis[i] && !isInX(i)) {
            res |= dfsCycles(i, -1);
            if (res)
                return res;
        }
    }
    return res;
}

bool Graph::dfsCycles(int v, int from) {
    vis[v] = 1;
    bool ret = false;

    for(int i = 0; i < neighbours[v].size(); i++) {
        int akt = neighbours[v][i];

        if (vis[akt] && !(akt == from && !isDoubleEdge(v, i))) {
            return true;
        }

        if (!vis[akt] && !isInX(akt)) {
            ret |= dfsCycles(akt, v);
        }
    }

    return ret;
}

bool Graph::isDoubleEdge(int v, int pos) {
    return neighDouble[v].size() > pos && neighDouble[v][pos];
}

bool Graph::checkSameCCNeighboursOneVertex(int i) {
    bool res = false;
    if (isNone(i)) {
        set<int> s;
        int uCnt = 0;
        bool put = false;
        for (int j = 0; j < neighbours[i].size(); j++) {
            int curNeigh = neighbours[i][j];

            if (isInU(curNeigh) && isDoubleEdge(i, j)) {
                res = true;
                put = true;
                putInSolution(i);
                break;
            }

            if (isInU(curNeigh)) {
                uCnt++;
                s.insert(fFind(neighbours[i][j]));
                if (s.size() < uCnt && !put) {
                    res = true;
                    put = true;
                    putInSolution(i);
                    break;
                }
            }
        }
        if (s.size() < uCnt && !put) {
            res = true;
            putInSolution(i);
        }
    }

    return res;
}

bool Graph::checkSameCCNeighbours() {
    bool res = false;

    for(int i = 0; i < neighbours.size(); i++) {
        res |= checkSameCCNeighboursOneVertex(i);
    }

    changeInU = false;
    return res;
}

int Graph::getNumberOfCCs(){
    return ccVertices.size();
}

int Graph::getCCSize(int nr) {
    assert(nr >= 0 && nr < getNumberOfCCs());
    return ccVertices[nr].size();
}

int Graph::getOriginalNr(int nr) {
    return oldNumbers[nr];
}

void Graph::printOriginalNums() {
    for(int i = 0; i < neighbours.size(); i++) {
        cerr << i << " <- " << getOriginalNr(i) << "\n";
    }
}

int Graph::calcNotNoneInCC(int nr) {
    assert(nr >= 0 && nr < getNumberOfCCs());
    int cnt = 0;
    for(int i = 0; i < ccVertices[nr].size(); i++) {
        if(!isNone(ccVertices[nr][i]))
            cnt++;
    }
    return cnt;
}

bool Graph::isAnyNoneInCC(int nr) {
    return calcNotNoneInCC(nr) != getCCSize(nr);
}

int Graph::recalcComponents() {
    vis.clear();
    vis.resize(N, 0);
    ccNr.clear();
    ccNr.resize(N, -1);
    int nr = 0; 

    for(int i = 0; i < neighbours.size(); i++) {
        // Deleted vertices would only make mess
        if (!vis[i] && !isInX(i) && !isDeleted(i) && neighbours[i].size()) {
            dfsComp(i, nr);
            nr++;
        }
    }

    ccVertices.clear();
    ccVertices.resize(nr);
    int nr_max = 0, max_size = 0;

    for(int i = 0; i < neighbours.size(); i++) {
        if (isInX(i) || isDeleted(i) || !neighbours[i].size())
            continue;

        ccVertices[ccNr[i]].push_back(i);

        if(ccVertices[ccNr[i]].size() > max_size) {
            max_size = ccVertices[ccNr[i]].size();
            nr_max = ccNr[i];
        }
    }
    if (!(nr_max >= 0 && nr_max < nr)) {
        printGraph();
        cerr << "nr max = " << nr_max << " nr = " << nr << std::endl;
        assert(false);
    }
    return nr_max;
}

void Graph::dfsComp(int v, int nr) {
    vis[v] = 1;
    ccNr[v] = nr;
    assert(!isDeleted(v));

    for(int i = 0; i < neighbours[v].size(); i++) {
        int akt = neighbours[v][i];

        if (!vis[akt] && !isInX(akt)) {
            dfsComp(akt, nr);
        }
    }
}

vector<pair<int, int>> Graph::createAllEdgeList() {
     vector<pair<int, int> > res;

    for(int i = 0; i < N; i++) {
        int cur = i;
        for(int j = 0; j < neighbours[cur].size(); j++) {
            int nei = neighbours[cur][j];
            // We put only from the vertex with smaller number -> we put each edge only once.
            if(cur < nei && !isInX(cur) && !isInX(nei)) {
                res.push_back(make_pair(cur, nei));
                if(isDoubleEdge(cur, j))
                    res.push_back(make_pair(cur, nei));
            }
        }
    }

    return res;
}


vector<pair<int, int>> Graph::createCCEdgeList(int nr) {
    assert(nr >= 0 && nr < ccVertices.size());

    vector<pair<int, int> > res;

    for(int i = 0; i < ccVertices[nr].size(); i++) {
        int cur = ccVertices[nr][i];

        for(int j = 0; j < neighbours[cur].size(); j++) {
            int nei = neighbours[cur][j];
            // We put only from the vertex with smaller number -> we put each edge only once.
            if(cur < nei && !isInX(cur) && !isInX(nei)) {
                res.push_back(make_pair(cur, nei));
                if(isDoubleEdge(cur, j))
                    res.push_back(make_pair(cur, nei));
            }
        }
    }

    return res;
}

vector<int> Graph::getCC(int nr) {
    assert(nr >= 0 && nr < ccVertices.size());
    return ccVertices[nr];
}

void Graph::printCC(int nr) {
   cerr << nr << ": ";
   for(int j = 0; j < ccVertices[nr].size(); j++) {
            cerr << ccVertices[nr][j] << " ";
    }
    cerr << "\n";
}

void Graph::printCCs() {
    for(int i = 0; i < ccVertices.size(); i++) {
        printCC(i);
    }
}

int Graph::fFind(int x) {
    return (fuRepr[x] < 0 ? x : fuRepr[x] = fFind(fuRepr[x]));
}

void Graph::fUnion(int x, int y) {
    if ((x = fFind(x)) == (y = fFind(y)))
        return;
    if (fuWeight[x] > fuWeight[y])
        fuRepr[y] = x;
    else 
        fuRepr[x] = y;

    if(fuWeight[x] == fuWeight[y])
        fuWeight[y]++;
}

pair<int, int> Graph::getFAUVals(int x) {
    return make_pair(fuRepr[x], fuWeight[x]);
}

void Graph::mergeWithNeighinU(int v) {
    for(int i = 0; i < neighbours[v].size(); i++) {
        int akt = neighbours[v][i];
        if (isInU(akt)) {
            fUnion(v, akt);
        }
    }
}

int Graph::allDoubleCnt() {
    int res = 0;
    for(int i = 0; i < neighbours.size(); i++)
        res += doubleNr2(i);
    return res;
}

bool Graph::reduceTooManyDouble(int bound) {
    bool res = false;
    for(int i = 0; i < neighbours.size(); i++) {
        if(isNone(i) && doubleNr2(i) > bound) {
            res = true;
            bound--;
            putInSolution(i);
        }
    }
    return res;
}
void Graph::debugCheckReduced() {
    for(int i = 0; i < neighbours.size(); i++) {
        if (isNone(i) && degree(i) <= 2 && doubleNum[i] < 2) {
            cerr << "Not REDUCED: " << i << "\n";
            printGraph();
            assert(false);
        }
    }
}

void Graph::populateDeg2Queue() {
    for (int i = 0; i < neighbours.size(); ++i)
        if ((isNone(i) || isInU(i)) && degree(i) == 2)
            deg2ToReduce.push(i);
}

void Graph::debugCheckDeg3Applicable() {
    for (int i = 0; i < neighbours.size(); ++i) {
        if (isNone(i)) {
            assert(degree(i) == 3);
            assert(degWithDouble(i) == 3);
            assert(doubleNum[i] == 0);
            for (auto j : neighbours[i]) {
                assert(isInU(j));
            }
        }
    }
}

bool Graph::isPotentialTent(int v) {
    if (!isNone(v))
        return false;
    if (degree(v) != 3 || degWithDouble(v) != 3)
        return false;
    int cnt_u = 0;
    for (auto u : neighbours[v]) {
        if (isInU(u))
            cnt_u++;
    }
    return cnt_u == 2;
}

bool Graph::reduceTents() {
    vector<int> queue;
    bool ret = false;

    for (int i = 0; i < getN(); ++i)
        if (isPotentialTent(i))
            queue.push_back(i);
    while(!queue.empty()) {
        int v = queue.back();
        queue.pop_back();
        if (!isPotentialTent(v))
            continue;
        int u = -1, vid = -1;
        for (int i = 0; i < 3; ++i)
            if (isNone(neighbours[v][i])) {
                vid = i;
                u = neighbours[v][i];
            }
        assert(u >= 0 && vid >= 0);
        int uid = posAtNeigh[v][vid];
        int w = getN();
        N++;
        neighbours.push_back(vector<int>());
        status.push_back(STATUS_IN_SAFE);
        neighDouble.push_back(vector<bool>());
        posAtNeigh.push_back(vector<int>());
        doubleNum.push_back(0);
        fuRepr.push_back(-1);
        fuWeight.push_back(0);
        oldNumbers.push_back(-1);

        neighbours[w].push_back(v);
        neighDouble[w].push_back(false);
        posAtNeigh[w].push_back(vid);
        neighbours[v][vid] = w;
        posAtNeigh[v][vid] = 0;

        neighbours[w].push_back(u);
        neighDouble[w].push_back(false);
        posAtNeigh[w].push_back(uid);
        neighbours[u][uid] = w;
        posAtNeigh[u][uid] = 1;

        ret = true;

        if (isPotentialTent(u))
            queue.push_back(u);
    }
    return ret;
}

int Graph::degLowerBound() {
    int m = 0, n = 0;
    vector<int> degrees;
    for (int i = 0; i < getN(); ++i) {
        if (isNone(i) || isInU(i)) {
            m += degWithDouble(i);
            n++;
        }
        if (isNone(i))
            degrees.push_back(degWithDouble(i));
    }
    m /= 2;
    std::sort(degrees.begin(), degrees.end());
    std::reverse(degrees.begin(), degrees.end());
    for (int i = 0; i < (int)degrees.size(); ++i) {
        if (m < n)
            return i;
        m -= degrees[i];
        n--;
    }
    return (int)degrees.size();
}
