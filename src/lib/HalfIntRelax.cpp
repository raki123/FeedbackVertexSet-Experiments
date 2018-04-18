#include "HalfIntRelax.h"

#include <cassert>
#include <iostream>
#include <set>
#include <queue>
#include <cmath>

using namespace std;

bool isEq(double a, double b) {
	return abs(a - b) <= EPS;
}
 
HalfIntRelax::HalfIntRelax(Graph& g, int s_) : s(s_), curVis(0) {
    buildLocalGraph(g);
	visi.resize(N, -1);
    type[s] = vType::s;
}

bool HalfIntRelax::isType(int v, vType t) {
    return type[v] == t;
}

void HalfIntRelax::buildLocalGraph(Graph& g) {
    N = g.getN();
    vector<pair<int, int>> E = g.createAllEdgeList();
    M = E.size();
	neighbours.resize(N);
	numberE.resize(N);
    for (int i = 0; i < E.size(); i++) {
        int a = E[i].first, b = E[i].second;
        neighbours[a].push_back(b);
        neighbours[b].push_back(a);
        numberE[a].push_back(i);
        numberE[b].push_back(i);
    }

    what.resize(M, eWhat::O);
    type.resize(N, vType::O);
    prev.resize(N, NOEDGE);
}

/*  
 * We hold variable curVis denoting what is the number of current dfs
 * Thanks to it we don't have to reset whole visited vector in every dfs
 */
 
void HalfIntRelax::prepareNewVisit() {
	curVis++;
}

bool HalfIntRelax::visited(int v) {
	return visi[v] == curVis;
}

void HalfIntRelax::visit(int v) {
	visi[v] = curVis;
}

void HalfIntRelax::zeroPrev() {
	for(int i = 0; i < prev.size(); i++)
        	prev[i] = NOEDGE;
}

vector<double> HalfIntRelax::getCover() {
	return cover;
}

/* 
 * Utilities for branching and kernelization
 */

vector<int> HalfIntRelax::getX1() {
	assert(cover.size());
	vector<int> res;
	for(int i = 0; i < cover.size(); i++)
		if(cover[i] == 1.0)
			res.push_back(i);

	return res;
}

vector<int> HalfIntRelax::getReachableThrough0() {
	assert(cover.size());
	vector<int> res;
	prepareNewVisit();
	queue<int> q;
	visit(s);
	q.push(s);

	while(!q.empty()) {
		int cur = q.front();
		q.pop();

		for(int i = 0; i < neighbours[cur].size(); i++) {
			int nei = neighbours[cur][i];
			if (!visited(nei) && cover[nei] == 0.0) {;
				visit(nei);
				q.push(nei);
				res.push_back(nei);
			}
		}
	}

	return res;
}

bool HalfIntRelax::isTreeCC(int v, int edge) {
	visit(v);

	for(int i = 0; i < neighbours[v].size(); i++) {
		int nei = neighbours[v][i];
		int eNei = numberE[v][i];
		if ((visited(nei) && (edge != eNei)))  // It takes care also of hitting s
			return false; 

		if (!visited(nei) && cover[nei] < 1.0) {
			if(!isTreeCC(nei, eNei)) 
				return false;
		}
	}

	return true;
}

vector<int> HalfIntRelax::getBridgesToTreeCC() {
	assert(cover.size());
	vector<int> bridges;

	prepareNewVisit();
	visit(s);
	
	for(int i = 0; i < neighbours[s].size(); i++) {
		int nei = neighbours[s][i];
		int eNei = numberE[s][i];
		if (isTreeCC(nei, eNei))
			bridges.push_back(nei);
	}
	return bridges;
}

double HalfIntRelax::getCoverSum(){
	return coverSum;
}

/*
 * Functions for transforming to cover from packing
 */

void HalfIntRelax::findCover() {
	coverSum = 0.0;
	findPacking();
	transformToCover();
	for(int i = 0; i < cover.size(); i++)
		coverSum += cover[i];
}

void HalfIntRelax::transformToCover() {
	cover.resize(N, 0);
	prepareNewVisit();
	visit(s);

	for(int i = 0; i < neighbours[s].size(); i++) {
		int nei = neighbours[s][i];
		int eNei = numberE[s][i];

		if(what[eNei] == eWhat::I) {
			if(!visited(nei)) {
				bool was = false;
				coverDfs(nei, eNei, was);
			}
		}
		else { // + because of double edges
			if(what[eNei] != eWhat::O)
				cover[nei] += 0.5;
		}
	}
	assert(checkCover());
}

void HalfIntRelax::coverDfs(int v, int edge, bool& was) {
	visit(v);

	if (prev[v] != NOEDGE && what[prev[v]] == eWhat::O) {
		cover[v] = 1.0;
		was = true;
	}

	int isSNei = 0;

	for(int i = 0; i < neighbours[v].size(); i++) {
		int nei = neighbours[v][i];
		int eNei = numberE[v][i];

		if (what[eNei] == eWhat::I && !visited(nei)) { // s should be marked visited already
			coverDfs(nei, eNei, was);
		}

		if (what[eNei] == eWhat::I && nei == s)
			isSNei++;
	}

	// If that was pure path, not touched during last augmentation => we put 0.5 on its ends
	if (!was) 
		cover[v] = isSNei * 0.5;

}

/*
 * Main part of algorithm - finding half integral packing.
 */

void HalfIntRelax::findPacking() {
    for (int i = 0; i < N; i++) type[i] = vType::O;
    type[s] = vType::s;

    for (int i = 0; i < M; i++) what[i] = eWhat::O;

    bool stop;
    do {
        stop = true;
        zeroPrev();
        
        for (int i = 0; i < neighbours[s].size(); i++) {
            int cur = neighbours[s][i];
            int nr = numberE[s][i];

            if (what[nr] == eWhat::O) {
                if (findAugmentingPath(cur, nr)) {
                    stop = false;
                    update();
                    break;
                }
            }
        }
    } while (!stop);
}

bool HalfIntRelax::findAugmentingPath(int v, int edge) {

    if ((isType(v, vType::N) || isType(v, vType::H) || isType(v, vType::T)) ||
        (isType(v, vType::s) && what[edge] == eWhat::O) ||
        (prev[v] != -1 &&
         !(what[prev[v]] == eWhat::O && what[edge] == eWhat::O && isType(v, vType::I)))) {
    	stos.push_back({v, edge});
    	return true;
    }

    if (isType(v, vType::s) || prev[v] != NOEDGE)
    	return false;

    prev[v] = edge;

    for (int i = 0; i < neighbours[v].size(); i++) {
    	int u = neighbours[v][i];
    	int uENr = numberE[v][i];
    	if (uENr != edge &&  	   // 2. Hit path and didn't make a single slide through it							
    		 !(what[edge] == eWhat::O  && what[uENr] == eWhat::O && isType(v, vType::I))) {
    		if (findAugmentingPath(u, uENr)) {
    			stos.push_back({v, edge});
    			return true;
			}
    	}
    }

    return false;
}

void HalfIntRelax::update_hit_flower() {
	//We are changing number of flowers legs to even and than decompose flower into paths
	int v = stos[0].first;
	int eFrom = stos[0].second;
	xorStack(0, stos.size() - 1);

	if (type[v] == vType::T || type[v] == vType::N) {
		// Make path from the stack + leg we hit
		// Clear the rest of the leg 
		int legEnd = clearLeg(v, eFrom);

		assert(legEnd != -1);

		if (type[v] == vType::T) {
			type[v] = vType::I;
			what[eFrom] = eWhat::I;
		}

		// Decompose the flower into paths (end of the leg)
		decomposeFlower(legEnd, eFrom, false);
	}

	if (type[v] == vType::H) {
		// Make a leg out from the stack
		assert(reachesSThroughI(v));
		findAndChange(v, vType::s, vType::N, eWhat::I, eWhat::N);
		type[v] = vType::T;

		// Decompose the flower into paths, include the leg
		decomposeFlower(v, eFrom, true);
	}
}

void HalfIntRelax::decomposeFlower(int v, int edge, bool includeBeg) {
	int legCnt = countLegs(v);
	int curLeg = -1, nextLeg = -1, firstLeg;

	vector<pair<int, int> > cycle = getCycle(v, edge);
	assert(cycle.size());

	int op = 0;

	if (includeBeg) {
		//We hit H -> there should exist a leg ending at v now
		curLeg = 0;
	}
	else {
		// We hit T or N -> the legs ending at v were erased
		curLeg = findNextLeg(0, cycle);
		for(int i = 1; i <= curLeg; i++) 
			cycle.push_back(cycle[i]);

		if (curLeg < 0) { // Cycle doesn't have legs anymore (everything was cleared already)
			// We want to clear it
			curLeg = -1;
			op = 1;
			legCnt = 1;
		}
	}

	while(legCnt > 0) {
		nextLeg = findNextLeg(curLeg, cycle);
		if (curLeg != cycle.size() - 1 && nextLeg < 0)
			nextLeg = cycle.size() - 1;
		assert(nextLeg > 0);

		if (!op) {
			connectWithI(curLeg, nextLeg, cycle);
			findAndChange(cycle[curLeg].first, vType::s, vType::I, eWhat::N, eWhat::I);
			findAndChange(cycle[nextLeg].first, vType::s, vType::I, eWhat::N, eWhat::I);
		}
		else {
			connectWithO(curLeg, nextLeg, cycle);
		}

		curLeg = nextLeg;
		legCnt--;
		op = (op + 1) % 2;
	}
}

int HalfIntRelax::getEdgeNr(int v1, int v2) {
	for(int i = 0; i < neighbours[v1].size(); i++)
		if (neighbours[v1][i] == v2)
			return numberE[v1][i];
	return -1;
}

void HalfIntRelax::connectWithI(int curLeg, int nextLeg, vector<pair<int, int>>& cycle) {
	for(int i = curLeg; i <= nextLeg; i++) {
		type[cycle[i].first] = vType::I;
		if (i > curLeg) {
			what[cycle[i].second] = eWhat::I;
		}
	}
}

void HalfIntRelax::connectWithO(int curLeg, int nextLeg, vector<pair<int, int>>& cycle) {
	for(int i = curLeg + 1; i <= nextLeg; i++) {
		if (i < nextLeg && type[cycle[i].first] != vType::I)
			type[cycle[i].first] = vType::O;
		if(i)
			what[cycle[i].second] = eWhat::O;
	}
}

int HalfIntRelax::findNextLeg(int pos, vector<pair<int, int> >& cyc) const {
	for(int i = pos + 1; i < cyc.size(); i++)
		if(type[cyc[i].first] == vType::T)
			return i;
	return -1;
}

vector<pair<int, int> > HalfIntRelax::getCycle(int v, int edge) {
	vector<pair<int, int> > res;
	prepareNewVisit();
	visit(v);
	res.push_back({v, edge});
	getCycleDFS(v, edge, res, v);
	return res;
}

void HalfIntRelax::getCycleDFS(int v, int edge, vector<pair<int, int>>& res, int start) {
	for(int i = 0; i < neighbours[v].size(); i++) {
		int nei = neighbours[v][i];
		int eNei = numberE[v][i];
		if(what[eNei] == eWhat::H  &&  (!visited(nei) || (nei == start && edge != eNei))) {
			res.push_back({nei, eNei});
			visit(nei);
			getCycleDFS(nei, eNei, res, start);
		}
	}
}

int HalfIntRelax::clearLeg(int v, int edge) {
	vector<pair<int, int> > st;
	st.push_back({v, edge});
	prepareNewVisit();
	visit(v);

	int legEnd = -1;

	if (type[v] == vType::T) {
		legEnd = v;
	}

	dfsClearLeg(v, edge, st, legEnd, v);

	return legEnd;
}

void HalfIntRelax::dfsClearLeg(int v, int edge, vector<pair<int, int>>& st, int& legEnd, int start) {
	int cur = v;
	int cur_e = edge;

	assert(what[cur_e] == eWhat::N || cur == start);

	if (cur != start && type[cur] == vType::T) {
		legEnd = cur;
		// Delete contents of part of leg from v to cur
		for(int i = 1; i < st.size() - 1; i++) {
			int vNr = st[i].first;
			int eNr = st[i].second;
			type[vNr] = vType::O;
			what[eNr] = eWhat::O;
		}

		what[cur_e] = eWhat::O;
		type[cur] = vType::H; 
		st.pop_back();
		return;
	}

	if (type[cur] == vType::s) {
		// Change type from N to I
		for(int i = 0; i < st.size() - 1; i++) {
			int vNr = st[i].first;
			int eNr = st[i].second;
			type[vNr] = vType::I;
			what[eNr] = eWhat::I;
		}

		what[cur_e] = eWhat::I;
		st.pop_back();
		return;
	}

	for(int i = 0; i < neighbours[cur].size(); i++) {
		int nei = neighbours[cur][i];
		int eNei = numberE[cur][i];

		if(!visited(nei) && what[eNei] == eWhat::N) {
			st.push_back({nei, eNei});
			visit(nei);
			dfsClearLeg(nei, eNei, st, legEnd, start);
		}
	}	

	st.pop_back();
}

int HalfIntRelax::countLegs(int v) {
	// v is supposed to be on cycle or have neighbours on cycle
	int res = 0;
	queue<int> q;
	q.push(v);
	prepareNewVisit();
	visit(v);

	while(!q.empty()) { 
		int cur = q.front();
		q.pop();

		if (type[cur] == vType::T)
			res++;

		for(int i = 0; i < neighbours[cur].size(); i++) {
			int nei = neighbours[cur][i];
			int eNei = numberE[cur][i];

			if (!visited(nei) && what[eNei] == eWhat::H) {
				visit(nei);
				q.push(nei);
			}
		}
	}
	return res;
}

void HalfIntRelax::update_hit_s() {
	xorStack(0, stos.size() - 1);
	type[s] = vType::s;
}


void HalfIntRelax::xorStack(int beg, int end) {
	for(int i = end; i >= beg; i--) {
		if(what[stos[i].second] == eWhat::O) {
			what[stos[i].second] = eWhat::I;
		}
		else if(what[stos[i].second] == eWhat::I) {
			what[stos[i].second] = eWhat::O;
		}
	}

	for(int i = end; i > beg; i--) {
		//update type[v]
		// 2 edges to v && both == O -> type[v] = O, else type[v] = I
		int cur_v = stos[i].first;
		int cur_e = stos[i].second; // Edge we came through to the vertex
		int next_e = stos[i - 1].second; // Edge we go out through to the next v
		if (what[next_e] == eWhat::O && what[cur_e] == eWhat::O)
			type[cur_v] = vType::O;
		else
			type[cur_v] = vType::I;		
	}
}

void HalfIntRelax::update_hit_self() {
    // 1. find cycle
    int cycBeg = -1, cycEnd = 0; // cycle ends at stos[0]

    for (int i = stos.size() - 1; i >= 0; i--) {
    	if (stos[i].first == stos[0].first) {
    		cycBeg = i;
    		break;
    	}
    }

    assert(cycBeg > 0);

    // Set the first vertex and edge in the cycle
    // stos[cycBeg].second - the last edge going into the cycle
    // if it wasn't on I path -> after xoring it will become a leg
    // Another option - if it was on path, that intersects with cycle 
    // only on the first vertex

    if (what[stos[cycBeg].second] == eWhat::O || 
    	(what[stos[cycBeg].second] == eWhat::I 
    		&& what[stos[0].second] == eWhat::O 
    		&& what[stos[cycBeg-1].second] == eWhat::O))
    	type[stos[0].first] = vType::T;
    else
    	type[stos[0].first] = vType::H;

 
	// cycle -> turn into flower with last leg coming in
    for (int i = cycBeg - 1; i >= 0; i--) {
    	// That was the first or last v on intersection of I and cycle 
    	// -> it will have a leg coming in
    	if (i) {
	    	if(isLegEnd(i))
	    		type[stos[i].first] = vType::T;
	    	else
	    		type[stos[i].first] = vType::H;
    	}

    	what[stos[i].second] = eWhat::H;
    }

    // Path to cycle -> xor all of it

    xorStack(cycBeg, stos.size() - 1);

   // Change all I coming into cycle to Legs (N)
	for(int i = 0; i < cycBeg; i++) {
		if (type[stos[i].first] == vType::T) {
			if (reachesSThroughI(stos[i].first)) {
				int end = findAndChange(stos[i].first, vType::s, vType::N, eWhat::I, eWhat::N);
				assert(end == s);
			}
			else {
				type[stos[i].first] = vType::H;
				// Delete the path through I if it led back to the cycle
				int otherEnd = findAndChange(stos[i].first, vType::T, vType::O, eWhat::I, eWhat::O);
				type[otherEnd] = vType::H;
			}
		}
	}

    //Clear unreachable I's, may be required if we slide through the same path I twice during 
    // one augmentation dfs
    clearUnreachableI();
}

bool HalfIntRelax::isLegEnd(int posAtStack) {
// Works before updating the states of edges
	if(posAtStack == 0) return false;

	//We came through stos[posAtStack].second, go to next through stos[posAtStack-1].second
	// If that was the first or last v on intersection of I and cycle 
	// -> it will have a leg coming in

	return ((what[stos[posAtStack].second] == eWhat::I && what[stos[posAtStack - 1].second] == eWhat::O) ||
			(what[stos[posAtStack].second] == eWhat::O && what[stos[posAtStack - 1].second] == eWhat::I));
}


bool HalfIntRelax::reachesSThroughI(int v) {
	prepareNewVisit();
// We don't want it to hit the vertex we came from
	visit(v);


	for (int i = 0; i < neighbours[v].size(); i++) {
		int eNr = numberE[v][i];
		int actNei = neighbours[v][i];
		if (!visited(actNei) && what[eNr] == eWhat::I) {
			return reachSDfs(actNei, eNr);
		}
	}
	// we should never give a vertex without I edge in neighbour
	assert(false);
	return false;
}

bool HalfIntRelax::reachSDfs(int v, int edge) {
	visit(v);
	if (v == s)
		return true;

	// May give a little speedup, it should work correctly without it
	if (what[edge] == eWhat::I && type[v] != vType::I)
		return false;

	for (int i = 0; i < neighbours[v].size(); i++) {
		int eNr = numberE[v][i];
		int actNei = neighbours[v][i];
		if (!visited(actNei) && what[eNr] == eWhat::I)
			return reachSDfs(actNei, eNr);
	}

	return false;
}

// Returns some vertex of type target reachable through oldE
int HalfIntRelax::findAndChange(int v, vType target, vType newT, eWhat oldE, eWhat newE) {
// it is different function than dfs, because we don't want to change type in v
	int res = -1;
	prepareNewVisit();
	visit(v);

	for (int i = 0; i < neighbours[v].size(); i++) {
		int eNr = numberE[v][i];
		int actNei = neighbours[v][i];

		if (!visited(actNei) && what[eNr] == oldE) {
			int curRes = dfsTypeChange(actNei, eNr, target, newT, oldE, newE);
			if(curRes != -1)
				res = curRes;
		}
	}

	return res;
}

int HalfIntRelax::dfsTypeChange(int v, int edge, vType targetType, vType newT, eWhat oldE, eWhat newE) {
	int res = -1;
	what[edge] = newE;
	if(type[v] == targetType)
		return v;
	type[v] = newT;
	visit(v);

	for (int i = 0; i < neighbours[v].size(); i++) {
		int eNr = numberE[v][i];
		int actNei = neighbours[v][i];
		if (!visited(actNei) && what[eNr] == oldE) {
			int curRes = dfsTypeChange(actNei, eNr, targetType, newT, oldE, newE);
			if (curRes != -1)
				res = curRes;
		}
	}
	return res;
}

//Function clears old parts of paths that no longer are reachable from s
void HalfIntRelax::clearUnreachableI() {
	findAndChange(s, vType::s, vType::I, eWhat::I, eWhat::I);

	for (int i = 0; i < neighbours.size(); i++) {
		if (type[i] == vType::I && !visited(i)) {
			for(int j = 0; j < neighbours[i].size(); j++) {
				int eNr = numberE[i][j];
				what[eNr] = eWhat::O;
			}
			type[i] = vType::O;
		}
	}
}

void HalfIntRelax::update() {
	assert(!stos.empty());

	pair<int, int> p = stos[0];
	int v = p.first, e = p.second;

	if (isType(v, vType::N) || isType(v, vType::H) || isType(v, vType::T))
		update_hit_flower();
	else if (v == s) 
		update_hit_s();
	else
		update_hit_self();	

	stos.clear();
}

/* 
 *  Utilities for checking correctness of solution.
*/

vType HalfIntRelax::brutVType(int v) {
	map<eWhat, int> mapa;
	int cnt = 0;
	mapa[eWhat::O] = 0;
	mapa[eWhat::N] = 0;
	mapa[eWhat::H] = 0;
	mapa[eWhat::I] = 0;

	for(int i = 0; i < numberE[v].size(); i++) {
		mapa[what[numberE[v][i]]]++;
		if(what[numberE[v][i]] != eWhat::O)
			cnt++;
	}

	if (v == s)
		return vType::s;

	if (cnt == 2 && mapa[eWhat::I] == 2) 
		return vType::I;

	if (cnt == 2 && mapa[eWhat::N] == 2)
		return vType::N;

	if (cnt == 2 && mapa[eWhat::H] == 2)
		return vType::H;

	if (cnt == 3 && mapa[eWhat::H] == 2 && mapa[eWhat::N] == 1)
		return vType::T;

	assert(cnt == 0);
	return vType::O;
}

bool HalfIntRelax::brutCheckAllTypes() {
	for(int i = 0; i < N; i++) {
		if(type[i] != brutVType(i)) {
			return false;
		}
	}
	return true;
}

bool HalfIntRelax::checkCover() {
	double sum = 0.0;

	for(int i = 0; i < cover.size(); i++)
		sum += cover[i];

	double cnt = 0;

	for(int i = 0; i < neighbours[s].size(); i++) {
		if (what[numberE[s][i]] != eWhat::O)
			cnt += 0.5;
	}
	assert(cnt == sum);

	prepareNewVisit();
	bool check = checkCoverDfs(s, -1);
	return check;
}

bool HalfIntRelax::checkCoverDfs(int v, int edge) {
	visit(v);

	for(int i = 0; i < neighbours[v].size(); i++) {
		int nei = neighbours[v][i];
		int eNei = numberE[v][i];
		if (eNei != edge && cover[nei] < 1.0) {
			if(visited(nei) || nei == s || (isEq(cover[nei], 0.5) && v != s)) {
				return false;
			}
			if (isEq(cover[nei], 0.0) /*&& !vis[nei] */)
				if(!checkCoverDfs(nei, eNei))
				return false;
		}
	}

	return true;
}

// Debug functions
void HalfIntRelax::printStackTrace() {
	cerr << "\nStack:\n";
    for (int i = 0; i < stos.size(); i++) {
    	// v0, e
        cerr << "(" << stos[i].first << ", " << stos[i].second << ") \n";
    }
    cerr << "\n";
}

void HalfIntRelax::printVType(int i) {
	switch(type[i]) {
			// {O, I, N, H, T, s};
			case vType::O : cerr << "O "; break;
			case vType::I : cerr << "I "; break;
			case vType::N : cerr << "N "; break;
			case vType::H : cerr << "H "; break;
			case vType::T : cerr << "T "; break;
			case vType::s : cerr << "s "; break;
		}
}

void HalfIntRelax::printVTypes() {
	for(int i = 0; i < N; i++) {
		printVType(i);		
	}
	cerr << "\n";
}

void HalfIntRelax::printEType(int eNr) {
	switch(what[eNr]) {
			// {O, I, N, H};
		case eWhat::O : cerr << "O"; break;
		case eWhat::I : cerr << "I"; break;
		case eWhat::N : cerr << "N"; break;
		case eWhat::H : cerr << "H"; break;
	}
}

void HalfIntRelax::printETypes(int i) {
	cerr << i << ": ";
	for(int j = 0; j < neighbours[i].size(); j++) {
		cerr << "(" << neighbours[i][j] << ", " << numberE[i][j] << ", ";
		printEType(numberE[i][j]);
		cerr << "), ";
	}
	cerr << "\n";
}

void HalfIntRelax::printETypes() {
	for(int i = 0; i < neighbours.size(); i++) {
		printETypes(i);
	}
}

void HalfIntRelax::printAllDebug() {
	cerr << "s = " << s << "\n";
	printVTypes();
	printETypes();
	printStackTrace();
	if(cover.size())
		printCover();
}

void HalfIntRelax::printCover() {
	for(int i = 0; i < cover.size(); i++){
		cerr << cover[i] << " ";
	}
	cerr << "\n";
}

void HalfIntRelax::printPrev() {
	cerr << "Prevs = \n";
	for(int i = 0; i < neighbours.size(); i++) 
		cerr << prev[i] << " ";
	cerr << "\n";
}