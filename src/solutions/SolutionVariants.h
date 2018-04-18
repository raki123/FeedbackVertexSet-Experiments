#ifndef SOLUTION_VARIANTS_H
#define SOLUTION_VARIANTS_H

#include "../lib/Graph.h"
#include "Solution.h"

class Sol8 : public BranchingSolution {
public:

	Sol8() : BranchingSolution() { br2 = 0; };
	Sol8(int cap) : BranchingSolution(cap) { br2 = 0; };

    virtual Solution *clone(int cap) { return new Sol8(cap); }

    virtual void stepIntoUndeletableBranch(Graph &brG) { BranchingSolution::stepIntoUndeletableBranch(brG); br2++; }
    virtual void stepOutOfUndeletableBranch(Graph &brG) { BranchingSolution::stepOutOfUndeletableBranch(brG); br2--; }
   
   virtual bool impossible(Graph &brG) {
        return br2 >= 2 * size_cap;
    }

    virtual int chooseBranchingVertex(Graph &brG) {
        return chooseMaxDeg2(brG);
    }

protected:
    int br2;
};

class Sol8noCC : public Sol8 {
public:
    Sol8noCC() : Sol8() { };
    Sol8noCC(int cap) : Sol8(cap) { };

    virtual Solution *clone(int cap) { return new Sol8noCC(cap); }

    virtual ReductionResult reductions(Graph &brG) {
        return reduceMaxDeg3(brG);
    }

};

class ICBranchingSolution : public BranchingSolution {
public:
    ICBranchingSolution() : BranchingSolution() { };
    ICBranchingSolution(int cap) : BranchingSolution(cap) { };

    virtual Solution *clone(int cap) { return new ICBranchingSolution(cap); }

    virtual int chooseBranchingVertex(Graph &brG) {
        while (!branching_hints.empty()) {
            int v = branching_hints.back();
            branching_hints.pop_back();
            if (brG.isNone(v))
                return v;
        }
        return chooseMaxDegToUndeletable(brG);
    }
};

class Sol5 : public ICBranchingSolution {
public:
    Sol5() : ICBranchingSolution() { };
    Sol5(int cap) : ICBranchingSolution(cap) { };

    virtual Solution *clone(int cap) { return new Sol5(cap); }

    virtual ReductionResult reductions(Graph &brG) {
        return reduceConnectedComponents(brG);
    }
};

class Sol5noCC : public Sol5 {
public:
    Sol5noCC() : Sol5() { };
    Sol5noCC(int cap) : Sol5(cap) { };

    virtual Solution *clone(int cap) { return new Sol5noCC(cap); }

    virtual ReductionResult reductions(Graph &brG) {
        return RED_NOTHING;
    }
};

class Sol5Full : public Sol5 {
public:
    Sol5Full() : Sol5() { };
    Sol5Full(int cap) : Sol5(cap) { };

    virtual Solution *clone(int cap) { return new Sol5Full(cap); }

    virtual bool impossible(Graph &brG) {
        return  brG.degLowerBound() + (int)brG.getX().size() >= size_cap;
    }
};

class Sol3_62 : public ICBranchingSolution {
public:
    Sol3_62() : ICBranchingSolution() { };
    Sol3_62(int cap) : ICBranchingSolution(cap) { };

    virtual Solution *clone(int cap) { return new Sol3_62(cap); }

    virtual ReductionResult reductions(Graph &brG) {
        ReductionResult r;
        r = reduceConnectedComponents(brG);
        if (r != RED_NOTHING)
            return r;
        r = reduceTents(brG);
        if (r != RED_NOTHING)
            return r;
        r = reduceMaxDeg3(brG);
        if (r != RED_NOTHING)
            return r;
        return r;
    }
};

class Sol3_62noCC : public Sol3_62 {
public:
    Sol3_62noCC() : Sol3_62() { };
    Sol3_62noCC(int cap) : Sol3_62(cap) { };

    virtual Solution *clone(int cap) { return new Sol3_62noCC(cap); }

    virtual ReductionResult reductions(Graph &brG) {
        ReductionResult r;
        r = reduceTents(brG);
        if (r != RED_NOTHING)
            return r;
        r = reduceMaxDeg3(brG);
        if (r != RED_NOTHING)
            return r;
        return r;
    }
};

class Sol3_62Full : public Sol3_62 {
public:
    Sol3_62Full() : Sol3_62() { };
    Sol3_62Full(int cap) : Sol3_62(cap) { };

    virtual Solution *clone(int cap) { return new Sol3_62Full(cap); }

    virtual bool impossible(Graph &brG) {
        return  brG.degLowerBound() + (int)brG.getX().size() >= size_cap;
    }
};

class SolPreferDouble : public BranchingSolution {
public:

    SolPreferDouble() : BranchingSolution() { };
    SolPreferDouble(int cap) : BranchingSolution(cap) { };

    virtual Solution *clone(int cap) { return new SolPreferDouble(cap); }

    virtual int chooseBranchingVertex(Graph &brG) {
        return chooseMaxDegPreferDouble(brG);
    }

    virtual bool impossible(Graph &brG) {
        return  brG.degLowerBound() + (int)brG.getX().size() >= size_cap;
    }
};

class SolPreferUndel : public BranchingSolution {
public:

    SolPreferUndel() : BranchingSolution() { };
    SolPreferUndel(int cap) : BranchingSolution(cap) { };

    virtual Solution *clone(int cap) { return new SolPreferUndel(cap); }

    virtual int chooseBranchingVertex(Graph &brG) {
        return chooseMaxDegToUndeletable(brG);
    }

    virtual bool impossible(Graph &brG) {
        return  brG.degLowerBound() + (int)brG.getX().size() >= size_cap;
    }
};

class Sol8Full : public Sol8 {
public:
    Sol8Full() : Sol8() {};
    Sol8Full(int cap) : Sol8(cap) {};
    virtual Solution *clone(int cap) { return new HalfIntBranchingFull(cap); }

    virtual bool impossible(Graph &brG) {
        return br2 >= 2 * size_cap || brG.degLowerBound() + (int)brG.getX().size() >= size_cap;
    }
};

class Sol8FullNoDeg3 : public Sol8Full {
public:
    Sol8FullNoDeg3() : Sol8Full() { };
    Sol8FullNoDeg3(int cap) : Sol8Full(cap) { };

    virtual Solution *clone(int cap) { return new Sol8FullNoDeg3(cap); }
    
    virtual ReductionResult reductions(Graph &brG) {
        return reduceConnectedComponents(brG);
    }

    virtual bool impossible(Graph &brG) {
        return br2 >= 3 * size_cap || brG.degLowerBound() + (int)brG.getX().size() >= size_cap;
    }
};

class Sol8NoDeg3 : public Sol8 {
public:
    Sol8NoDeg3() : Sol8() { };
    Sol8NoDeg3(int cap) : Sol8(cap) { };

    virtual Solution *clone(int cap) { return new Sol8NoDeg3(cap); }

    virtual ReductionResult reductions(Graph &brG) {
        return reduceConnectedComponents(brG);
    }

    virtual bool impossible(Graph &brG) {
        return br2 >= 3 * size_cap;
    }
};


class SolPreferDoubleIC : public SolPreferDouble {
public:

    SolPreferDoubleIC() : SolPreferDouble() {};
    SolPreferDoubleIC(int cap) : SolPreferDouble(cap) {};

    virtual Solution *clone(int cap) { return new SolPreferDoubleIC(cap); }

    virtual int chooseBranchingVertex(Graph &brG) {
        while (!branching_hints.empty()) {
            int v = branching_hints.back();
            branching_hints.pop_back();
            if (brG.isNone(v))
                return v;
        }
        return chooseMaxDegPreferDouble(brG);
    }
};

class SolPreferUndelIC : public SolPreferUndel {
public:

    SolPreferUndelIC() : SolPreferUndel() {};
    SolPreferUndelIC(int cap) : SolPreferUndel(cap) {};

    virtual Solution *clone(int cap) { return new SolPreferUndelIC(cap); }

    virtual int chooseBranchingVertex(Graph &brG) {
        while (!branching_hints.empty()) {
            int v = branching_hints.back();
            branching_hints.pop_back();
            if (brG.isNone(v))
                return v;
        }
        return chooseMaxDegToUndeletable(brG);
    }
};

class Sol8FullIC : public Sol8Full {
public:
    Sol8FullIC() : Sol8Full() {};
    Sol8FullIC(int cap) : Sol8Full(cap) {};

    virtual Solution *clone(int cap) { return new Sol8FullIC(cap); }

    virtual int chooseBranchingVertex(Graph &brG) {
        while (!branching_hints.empty()) {
            int v = branching_hints.back();
            branching_hints.pop_back();
            if (brG.isNone(v))
                return v;
        }
        return chooseMaxDeg2(brG);
    }
};

#endif