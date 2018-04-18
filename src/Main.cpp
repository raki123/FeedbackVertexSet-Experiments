#include <cstring>
#include <iostream>
#include "lib/Graph.h"
#include "lib/Translator.h"
#include "solutions/HalfIntBranching.h"
#include "solutions/SolutionVariants.h"

using namespace std;

// Takes solution kind as argument
int main(int argc, char *argv[]) {
    Graph g;
    Translator t;

    t.read_input(&g);
    Solution *s;
    bool noCC = argc > 2 && !strcmp(argv[2], "noCC");
    bool full = argc > 2 && !strcmp(argv[2], "Full");
    bool IC = argc > 2 && (!strcmp(argv[2], "FullIC") || !strcmp(argv[2], "IC"));
    bool nodeg3 = argc > 2 && !strcmp(argv[2], "NoDeg3");
    bool deg3 = argc > 2 && !strcmp(argv[2], "Deg3");
    if (argc <= 1 || !strcmp(argv[1], "8")) {
        if (noCC)
            s = new Sol8noCC();
        else if (full)
            s = new Sol8Full();
        else if (IC)
            s = new Sol8FullIC();
        else if (nodeg3)
            s = new Sol8NoDeg3();
        else
            s = new Sol8();
    } else if (!strcmp(argv[1], "5")) {
        if (noCC)
            s = new Sol5noCC();
        else if (full)
            s = new Sol5Full();
        else
            s = new Sol5();
    } else if (!strcmp(argv[1], "3.62")) {
        if (noCC)
            s = new Sol3_62noCC();
        else if (full)
            s = new Sol3_62Full();
        else
            s = new Sol3_62();
    } else if (!strcmp(argv[1], "PD")) {
        if (IC)
            s = new SolPreferDoubleIC();
        else
            s = new SolPreferDouble();
    } else if (!strcmp(argv[1], "PU")) {
        if (IC)
            s = new SolPreferUndelIC();
        else
            s = new SolPreferUndel();
    } else if (!strcmp(argv[1], "Full")) {
        if (nodeg3)
            s = new Sol8FullNoDeg3();
        else
            s = new Sol8Full();
    } else if (!strcmp(argv[1], "LP")) {
        if (full)
            s = new HalfIntBranchingFull();
        else if (deg3)
            s = new HalfIntBranchingFullDeg3();
        else
            s = new HalfIntBranching();
    } else if (!strcmp(argv[1], "Appx")) {
        s = new ApproximateSolution();
    } else {
        cerr << "Wrong arguments!" << std::endl;
        return -1;
    }

    vector<int> sol = s->solve(g);
    delete s;
    cout << sol.size() << "\n";
    t.printNames(sol);

    return 0;
}
