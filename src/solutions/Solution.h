#ifndef SOLUTION_H
#define SOLUTION_H

#include "../lib/Graph.h"
#include "../lib/CubicSolver.h"
#include <vector>
#include <iostream>
#include <cassert>
using std::cerr;
using std::vector;

//#define DEBUG_CHECK_DOUBLE_EDGES
#define SPLIT_CC
#define USE_CUBIC_SOLVER

class Solution {
public:
    bool solFound;
    vector<int> branching_hints;
    virtual ~Solution(){};
  	virtual vector<int> solve(Graph& g) = 0;
    virtual Solution *clone(int) = 0;
};

enum ReductionResult {
    RED_NOTHING,
    RED_CHANGED,
    RED_NOSOLUTION,
};

class BranchingSolution : public Solution {
public:

    BranchingSolution() : isSub(false) { solFound = false; depth = 0; };
    BranchingSolution(int cap) : size_cap(cap), isSub(true) { solFound = false; depth = 0; };

    virtual Solution *clone(int cap) { return new BranchingSolution(cap); }

    vector<int> solve(Graph& g) {
        vector<int> apx = find_approx(g);
        if (!isSub || size_cap >= apx.size()) {
            size_cap = apx.size();
            sol = apx;
            solFound = true;
        }
        branching_hints = apx;
        std::reverse(branching_hints.begin(), branching_hints.end());
        branching(g);
        return sol;
    }

    virtual bool impossible(Graph &brG) {
        return false;
    }

    ReductionResult simpleReductions(Graph &brG) {
        bool change = false;
        // RR1.5 Delete selfloops
        change |= brG.deleteSelfLoops();
#ifdef DEBUG_CHECK_DOUBLE_EDGES
        for(int i = 0; i < brG.getN(); i++)
            assert(brG.doubleNr(i) == brG.doubleNr2(i));
#endif

        // RR2 v with deg <= 1 -> delete
        change |= brG.deleteDeg1Vertices();
#ifdef DEBUG_CHECK_DOUBLE_EDGES
        for(int i = 0; i < brG.getN(); i++)
            assert(brG.doubleNr(i) == brG.doubleNr2(i));
#endif

        // RR 2.5 Reduce deg 2 vertices
        change |= brG.reduceDeg2Queue();
#ifdef DEBUG_CHECK_DOUBLE_EDGES
        for(int i = 0; i < brG.getN(); i++)
            assert(brG.doubleNr(i) == brG.doubleNr2(i));
#endif

        // RR3  v with >= 2 edges to the same CC in U -> X
        change |= brG.checkSameCCNeighbours();
#ifdef DEBUG_CHECK_DOUBLE_EDGES
        for(int i = 0; i < brG.getN(); i++)
            assert(brG.doubleNr(i) == brG.doubleNr2(i));
#endif

        return change ? RED_CHANGED : RED_NOTHING;
    }

    // Reduces connected components. Returns
    ReductionResult reduceConnectedComponents(Graph &brG) {
        // Split into connected components
        int largestCC = brG.recalcComponents();
        int ccNum = brG.getNumberOfCCs();

        bool changed = false;

        for(int i = 0; i < ccNum; i++) {
            if (i == largestCC)
                continue;

            if(!brG.isAnyNoneInCC(i))
                continue;

            Graph curCC(brG, i);
            int ccSize = curCC.getN();

            BranchingSolution *s = (BranchingSolution *)clone(size_cap - brG.getXSize());
            s->solve(curCC);

            if (!s->solFound) {
                delete s;
                return RED_NOSOLUTION;
            }

            // Set all vertices to their statuses from the cc solution
            changed = true;

            // If empty, then no solution (exceeded size cap)
            if(s->sol.empty()) {
                curCC.printGraph();
                cerr << s->solFound << "  " << s->sol.size() << "   " << s->size_cap << " " << s->branching_hints.size() << std::endl;
                assert(false);
            }
            for (auto v : s->sol) {
                brG.purePutInSolution(curCC.getOriginalNr(v));
            }
            for (int j = 0; j < curCC.getNeighbours().size(); ++j) {
                int oldNr = curCC.getOriginalNr(j);
                if (!brG.isInX(oldNr)) {
                    brG.pureDeleteVertex(oldNr);
                }
            }
            delete s;
        }

        return changed ? RED_CHANGED : RED_NOTHING;
    }

    ReductionResult reduceMaxDeg3(Graph &brG) {
        pair<int, int> maxDeg = brG.getMaxDegVWDouble();
        if (maxDeg.second > 3)
            return RED_NOTHING;

        CubicSolver cs(brG.getStatuses(), brG.getNeighbours());

        vector<int> ret = cs.solve();
        for (auto a : ret)
            brG.purePutInSolution(a);

        for (int i = 0; i < brG.getN(); ++i)
            if (brG.isNone(i) || brG.isInU(i))
                brG.pureDeleteVertex(i);

        return RED_CHANGED;
    }

    ReductionResult reduceTents(Graph &brG) {
        return brG.reduceTents() ? RED_CHANGED : RED_NOTHING;
    }

    virtual ReductionResult reductions(Graph &brG) {
        ReductionResult r = RED_NOTHING;
#ifdef SPLIT_CC
        r = reduceConnectedComponents(brG);
        if (r != RED_NOTHING)
            return r;
#endif
#ifdef USE_CUBIC_SOLVER
        r = reduceMaxDeg3(brG);
        if (r != RED_NOTHING)
            return r;
#endif
        return r;
    }

    ReductionResult reduce(Graph &brG) {
        while(true) {
            if (brG.getXSize() >= size_cap || this->impossible(brG))
                return RED_NOSOLUTION;

            if (brG.isEmpty()) {
                trySolUpdate(brG);
                return RED_NOSOLUTION;
            }

            ReductionResult r = simpleReductions(brG);

            if (r == RED_NOSOLUTION)
                return r;
            else if (r == RED_CHANGED)
                continue;

            r = reductions(brG);
            if (r != RED_CHANGED)
                return r;
        }
    }

    int chooseMaxDeg(Graph &brG) {
        pair<int, int> maxDeg = brG.getMaxDegV();
        return maxDeg.first;
    }
    int chooseMaxDeg2(Graph &brG) {
        pair<int, int> maxDeg = brG.getMaxDegVWDouble();
        return maxDeg.first;
    }
    int chooseMaxDegToUndeletable(Graph &brG) {
        pair<int, int> maxDeg = brG.getMaxDegToUndeletable();
        return maxDeg.first;
    }
    int chooseMaxDegPreferDouble(Graph &brG) {
        pair<int, int> maxDeg = brG.getMaxDegPreferDouble();
        return maxDeg.first;
    }
    virtual int chooseBranchingVertex(Graph &brG) {
        return chooseMaxDeg(brG);
    }

    virtual void stepIntoDeletionBranch(Graph &brG) { depth++; }
    virtual void stepOutOfDeletionBranch(Graph &brG) { depth--; }
    virtual void stepIntoUndeletableBranch(Graph &brG) { depth++; }
    virtual void stepOutOfUndeletableBranch(Graph &brG) { depth--; }

    void branching(Graph& brG) {
        if (reduce(brG) == RED_NOSOLUTION)
            return;

        int v = chooseBranchingVertex(brG);
        if (v < 0)
            return;

        Graph g2(brG);

        brG.putInSolution(v);
        stepIntoDeletionBranch(brG);
        branching(brG);
        stepOutOfDeletionBranch(brG);

        // branch 2
        //  v -> U

        g2.putInSafe(v);
        stepIntoUndeletableBranch(brG);
        branching(g2);
        stepOutOfUndeletableBranch(brG);
    }

    vector<int> find_approx(Graph& g) {
        Graph gApx(g);

        while(gApx.checkCycles()) {
            ReductionResult r = simpleReductions(gApx);
            if (r == RED_CHANGED)
                continue;
            pair<int, int> md = gApx.getMaxDegV();
            int v = md.first;
            gApx.putInSolution(v);
        }

        return gApx.getX();
    }

    vector<int> getSolution() {
        return sol;
    }

    int getCap() {
        return size_cap;
    }

    void setCap(int c) {
        size_cap = c;
    }

    void trySolUpdate(Graph& brG) {
        if (brG.getXSize() < size_cap) {
            sol = brG.getX();
            size_cap = sol.size();
            solFound = true;
        }
    }


protected:
    vector<int> sol;
    bool isSub;
    int depth;
    int size_cap;
};

class ApproximateSolution : public BranchingSolution {
public:
    ApproximateSolution() { solFound = false; };

    virtual Solution *clone(int cap) { return new ApproximateSolution(); }

    virtual vector<int> solve(Graph &g) {
        branching_hints = find_approx(g);
        solFound = true;
        return branching_hints;
    }
};

#endif
