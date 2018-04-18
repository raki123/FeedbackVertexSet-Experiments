#ifndef HALFINTBRANCHING_H
#define HALFINTBRANCHING_H

#include "../lib/Graph.h"
#include "../lib/HalfIntRelax.h"
#include "Solution.h"

class HalfIntBranching : public BranchingSolution {
protected:
    bool kernelize_switch;

public:
    HalfIntBranching();
    HalfIntBranching(int);

    Solution *clone(int cap) { return new HalfIntBranching(cap); }
    virtual vector<int> solve(Graph &g);

    void solve1(Graph &g);
    void solve2(Graph &g, int s);

    // Main kernelization function from Iwatas paper
    ReductionResult kernelize(Graph &g);
    // Checks if total number of double edges is > k^2
    ReductionResult highTotalDoubleCnt(Graph &g);
    virtual ReductionResult reductions(Graph &g);

    // Returns how many more vertices we can add to the solution
    int getLeftK(Graph &g);
};

class HalfIntBranchingFull : public HalfIntBranching {
public:
    HalfIntBranchingFull() : HalfIntBranching(){};
    HalfIntBranchingFull(int cap) : HalfIntBranching(cap){};

    Solution *clone(int cap) { return new HalfIntBranchingFull(cap); }

    virtual bool impossible(Graph &brG) {
        return brG.degLowerBound() + (int)brG.getX().size() >= size_cap;
    }
};

class HalfIntBranchingFullDeg3 : public HalfIntBranchingFull {
public:
    HalfIntBranchingFullDeg3() : HalfIntBranchingFull(){};
    HalfIntBranchingFullDeg3(int cap) : HalfIntBranchingFull(cap){};

    Solution *clone(int cap) { return new HalfIntBranchingFullDeg3(cap); }

    virtual ReductionResult reductions(Graph &g) {
        ReductionResult r = HalfIntBranching::reductions(g);
        if (r != RED_NOTHING) return r;

        return reduceMaxDeg3(g);
    }
};

#endif
