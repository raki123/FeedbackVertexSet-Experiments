/*
 Copyright 2017 Marcin Pilipczuk.

 This file used to be a part of fvs_pace_challenge,
 an implementation of FPT algorithm for Feedback Vertex Set,
 a submission to track B of PACE Challenge 2016.
 It has been adjusted to fit the new graph interface.
*/

#ifndef CUBICSOLVER_H
#define CUBICSOLVER_H

#include <deque>
#include <cassert>
#include <set>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>

//#define DEG3DEBUGINFO

using namespace std;

class CubicSolver{
private:
    vector<int> nei;
    vector<int> blossom, backlabel, sid, mate;
    vector<vector<int> > undels;
    vector<pair<int, int> > transform;
    vector<int> fau;
    vector<bool> taken;
    vector<int> verts_names;
    int ne, nec, next_sid, nv, nu;
    unordered_map<int, int> names2ids;
    vector<set<int> > h;
    vector<vector<int> > blossoms, tips;
    deque<int> q;
    vector<int> comp_id, comp_roots;
    vector<vector<int> > path;

#ifdef DEG3DEBUGINFO
    vector<int> revlabel;
    vector<bool> old_taken;
    vector<vector<int> > old_augments;
    vector<pair<int, int> > blossom_starts, blossom_tips;
#endif

    const int CC_LABEL_OFFSET = -4;
    const int INF = 1000*1000*1001;

    inline int VJ2E(int v, int j);
    inline int E2V(int e);
    inline int E2J(int e);
    inline int NE1(int e);
    inline int NE2(int e);
    inline int CC2E(int i);
    inline int E2TFM(int e);
    inline bool IS_REGULAR(int e);
    inline bool IS_CC(int e);
    inline bool IS_TRANSFORM(int e);
    inline bool IS_LABELLED(int e);

    /*** Building dependency graph ***/
    bool build_from_not_taken_dfs(int u, int e, int e0, int target);
    void build_from_not_taken(int e);
    void build_from_component_edges();

    void build_dependency_graph(void);

    /*** Handling search paths and labeling ***/
    int flatten_edge(int e);
    void compute_path(int e, int back, int rev);
    void label(int e, int back, int rev);

    /*** Transforms and blossoms ***/
    /* Make transform (f1, f2) */
    int create_transform(int f1, int f2);

    /* Initiate blossom */
    int setup_blossom();

    /* Add vertex to blossom. Merges blossoms if required. */
    void add_vertex_to_blossom(int x, int bid, bool is_tip=false);

    /* Checks if e and f are in the same blossom. f can be negative (returns false then). */
    bool same_blossom(int e, int f);

    /* Perform blossom step. */
    void new_blossom(int e0, int e1, int b0);

    /* Augment solution along a search path. */
    void augment(vector<int> &path);

    /* One iteration of the augmenting path algorithm: returns if it managed to augment it. */
    bool iteration();

    /* Returns computed solution */
    vector<int> give_solution();

    /* Find and union implementation to check consistency & compute current components */
    int fau_find(int x);
    bool fau_join(int x, int y);

    /* Recomputes find and union, checking if the solution is feasible. */
    void check_solution();

    /* Floods undeletable connected component; needed for constructor. */
    void dfs_undeletable(int x, int cmp_id, vector<int> &vis, const vector<int> &status, const vector<vector<int> > &nei);

public:
    /* Constructs base graph 
       The first argument is a list of flags: <0 means ignore vertex (deleted), 0 means deletable, 1 undeletable
       The second one is a list of neighbours for every vertex.
       Every deletable vertex needs to have degree 3 and there cannot be parallel edges.
       These properties are NOT checked in the implementation.
     */
    CubicSolver(vector<int> &status, vector<vector<int> > &neighbors);

    /* Finds minimum feedback vertex set */
    vector<int> solve();
};

#endif //CUBICSOLVER_H
