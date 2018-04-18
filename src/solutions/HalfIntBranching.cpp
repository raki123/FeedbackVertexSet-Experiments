#include "HalfIntBranching.h"
#include <iostream>
using namespace std;

HalfIntBranching::HalfIntBranching() : BranchingSolution(), kernelize_switch(true) { }

HalfIntBranching::HalfIntBranching(int cap) : BranchingSolution(cap), kernelize_switch(true) {}

vector<int> HalfIntBranching::solve(Graph& g) {;
    vector<int> apx = find_approx(g);
    if (!isSub || size_cap >= apx.size()) {
        size_cap = apx.size();
        sol = apx;
        solFound = true;
    }

    // May land here from recursion on connected components
    // Then we should look for undeletable guys
    int s = -1, deg = -1;
    for (int i = 0; i < g.getN(); ++i)
      if (g.isInU(i) && g.degree(i) > deg) {
        s = i;
        deg = g.degree(i);
      }
    if (s >= 0)
        solve2(g, s);
    else
        solve1(g);

    return sol;
}

int HalfIntBranching::getLeftK(Graph& g) {
	return size_cap - g.getXSize();
}
 
ReductionResult HalfIntBranching::kernelize(Graph& g) {
    bool change = false;

    for(int s = 0; s < g.getN(); s++) {
        if(g.degree(s)) {
            HalfIntRelax hi(g, s);
            hi.findCover();

            vector<int> X1 = hi.getX1();
            vector<int> bridges = hi.getBridgesToTreeCC();
            
            if(hi.getCoverSum() > getLeftK(g)) {
               if (g.isNone(s))
                   g.putInSolution(s);
               else {
                   assert(g.isInU(s));
                   return RED_NOSOLUTION;
               }
               continue;
            }

            for(int i = 0; i < X1.size(); i++)
                change |= g.addEdge(s, X1[i], true);

            for(int i = 0; i < bridges.size(); i++) 
                g.removeEdge(s, bridges[i]);

            if(bridges.size())
                change = true;
        }
    }

    return change ? RED_CHANGED : RED_NOTHING;
}

ReductionResult HalfIntBranching::highTotalDoubleCnt(Graph &g) {
    int k = getLeftK(g);
    if (g.allDoubleCnt() > k*k)
        return RED_NOSOLUTION;
    return RED_NOTHING;
}

ReductionResult HalfIntBranching::reductions(Graph& g) {
    ReductionResult r;

    r = reduceConnectedComponents(g);
    if (r != RED_NOTHING)
        return r;

    if(g.reduceTooManyDouble(getLeftK(g)))
        return RED_CHANGED;

    r = highTotalDoubleCnt(g);
    if (r != RED_NOTHING)
        return r;

    if (kernelize_switch) {
        r = kernelize(g);
        if (r != RED_NOTHING)
            return r;
    }

    return RED_NOTHING;
}

void HalfIntBranching::solve1(Graph& g) {
    ReductionResult r = reduce(g);
    if (r == RED_NOSOLUTION)
        return;

    pair<int, int> md = g.getMaxDegV();
    int s = md.first;

    if (s == -1) 
        return;

    Graph g2(g);

    g.putInSolution(s);
    solve1(g);

    g2.putInSafe(s);
    solve2(g2, s);
}

void HalfIntBranching::solve2(Graph& g, int s) {
    kernelize_switch = false;
    ReductionResult r = reduce(g);
    kernelize_switch = true;
    if (r == RED_NOSOLUTION)
        return;

    if (!g.degree(s)) {
        solve1(g);
        return;
    }

    assert(g.isInU(s));

    HalfIntRelax hi(g, s);
    hi.findCover();
    
    if ((hi.getCoverSum() * 2) != g.degree(s)) {
        vector<int> toDel = hi.getX1();
        vector<int> toMerge = hi.getReachableThrough0();

        for(int i = 0; i < toDel.size(); i++)
            if(!g.isInX(toDel[i]))
                g.purePutInSolution(toDel[i]);

        for(int i = 0; i < toMerge.size(); i++)
            g.merge(s, toMerge[i]); 

        solve2(g, s);
    }
    else {
        pair<int, int> md = g.getMaxDegNeigh(s);
        int v = md.first;
        if (v == -1)
            return;
        
        Graph g2(g);

        g.putInSolution(v);
        solve2(g, s);

        g2.merge(s, v);
        solve2(g2, s);
    }
}
