#include "Translator.h"
#include <iostream>
#include <sstream>

using namespace std;

void Translator::read_input(Graph *g) {
    map<string, int> namesMap;
    vector<map<int, int> > neighMap;

    N = M = 0;
    string line;
    string s1, s2;

    vector<pair<int, int>> edges;

    while (getline(cin, line)) {
        if (line[0] == '#' || !line.size())
            continue;

        stringstream ss(line);

        ss >> s1;
        ss >> s2;
        int a = add_vertex(s1, namesMap);
        int b = add_vertex(s2, namesMap);
        edges.push_back(make_pair(a, b));
        M++;
    }

    N = namesMap.size();

    g->buildGraph(N, M, edges);
    g->populateDeg2Queue();
}

void Translator::printNames(vector<int> v) {
    for (int i = 0; i < v.size(); i++)
        cout << names[v[i]] << " ";
    cout << "\n";
}

string Translator::getName(int v) {
    return names[v];
}

vector<string> Translator::list_of_names(vector<int> v) {
    vector<string> res;
    for (int i = 0; i < v.size(); i++) {
        res.push_back(getName(v[i]));
    }
    return res;
}

int Translator::add_vertex(string name, map<string, int>& namesMap) {
    if (namesMap.count(name))
        return namesMap[name];
    names.push_back(name);
    int nr = names.size() - 1;
    namesMap[name] = nr;
    return nr;
}
