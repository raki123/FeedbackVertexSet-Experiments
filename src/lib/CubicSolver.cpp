/*
 Copyright 2017 Marcin Pilipczuk.

 This file used to be a part of fvs_pace_challenge,
 an implementation of FPT algorithm for Feedback Vertex Set,
 a submission to track B of PACE Challenge 2016.
 It has been adjusted to fit the new graph interface.
*/

#include <deque>
#include <cassert>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include "CubicSolver.h"

inline int CubicSolver::VJ2E(int v, int j){ return 3*v + j; }
inline int CubicSolver::E2V(int e){ assert(IS_REGULAR(e)); return e/3; }
inline int CubicSolver::E2J(int e){ assert(IS_REGULAR(e)); return e%3; }
inline int CubicSolver::NE1(int e){ assert(IS_REGULAR(e)); return VJ2E(E2V(e), (E2J(e)+1)%3); }
inline int CubicSolver::NE2(int e){ assert(IS_REGULAR(e)); return VJ2E(E2V(e), (E2J(e)+2)%3); }
inline int CubicSolver::CC2E(int i){ return i + ne; }
inline int CubicSolver::E2TFM(int e){ assert(IS_TRANSFORM(e)); return e - nec; }
inline bool CubicSolver::IS_REGULAR(int e){ return e >= 0 && e < ne; }
inline bool CubicSolver::IS_CC(int e){ return e >= ne && e < nec; }
inline bool CubicSolver::IS_TRANSFORM(int e){ return e >= nec && e < nec + (int)transform.size(); }

inline bool CubicSolver::IS_LABELLED(int e){ return backlabel[e] <= CC_LABEL_OFFSET || backlabel[e] >= 0; }


bool CubicSolver::build_from_not_taken_dfs(int u, int e, int e0, int target){
    bool res = false;
    if (u == target)
        res = true;
    else
        for (__typeof(undels[u].begin()) it = undels[u].begin(); it != undels[u].end(); ++it){
            int f = *it;
            if (!taken[E2V(f)])
                continue;
            if (E2J(f) == 0){
                if (e != NE2(f))
                    res |= build_from_not_taken_dfs(nei[NE1(f)], NE2(f), e0, target);
                if (e != NE1(f))
                    res |= build_from_not_taken_dfs(nei[NE2(f)], NE1(f), e0, target);
            } else {
                int f1 = VJ2E(E2V(f), 3-E2J(f));
                if (e != f1)
                    res |= build_from_not_taken_dfs(nei[VJ2E(E2V(f), 0)], f1, e0, target);
            }
        }
    if (res && e >= 0){
        assert(E2J(e) != 0);
        h[e].insert(e0);
        h[e0].insert(e);
    }
    return res;
}

void CubicSolver::build_from_not_taken(int e){
    int u1 = fau_find(nei[NE1(e)]);
    int u2 = fau_find(nei[NE2(e)]);
    if (u1 == u2)
        build_from_not_taken_dfs(nei[NE1(e)], -1, e, nei[NE2(e)]);
    else{
        h[e].insert(CC2E(comp_id[u1]));
        h[e].insert(CC2E(comp_id[u2]));
        build_from_not_taken_dfs(nei[NE1(e)], -1, e, u1);
        build_from_not_taken_dfs(nei[NE2(e)], -1, e, u2);
    }

}

void CubicSolver::build_from_component_edges(){
    for (int e = 0; e < ne; ++e)
        if (E2J(e) != 0 && !taken[E2V(e)]){
            int f1 = fau_find(nei[NE1(e)]);
            int f2 = fau_find(nei[NE2(e)]);
            if (f1 != f2){
                h[CC2E(comp_id[f1])].insert(e);
                h[CC2E(comp_id[f2])].insert(e);
            }
        }
}

void CubicSolver::build_dependency_graph(void){
    h.clear();
    h.resize((unsigned)nec);
    build_from_component_edges();
    for (int e = 0; e < ne; ++e)
        if (E2J(e) != 0 && !taken[E2V(e)])
            build_from_not_taken(e);
}

int CubicSolver::flatten_edge(int e){
    assert(IS_TRANSFORM(e) || IS_REGULAR(e));
    if (IS_TRANSFORM(e))
        return transform[E2TFM(e)].first;
    else
        return e;
}

void CubicSolver::compute_path(int e, int back, int rev){
    assert(path[e].size() == 0);
    if (back >= 0)
        path[e] = path[back];
    else
        path[e].push_back(back);
    if (rev >= 0){
        for (__typeof(path[rev].rbegin()) it = path[rev].rbegin(); it != path[rev].rend(); it++) {
            if (flatten_edge(e) == flatten_edge(*it))
                break;
            else
                path[e].push_back(*it);
        }
    } else
        path[e].push_back(mate[e]);
    path[e].push_back(e);
}

void CubicSolver::label(int e, int back, int rev){
    if (backlabel[e] == -1){
        backlabel[e] = back;
#ifdef DEG3DEBUGINFO
        revlabel[e] = rev;
#endif
        compute_path(e, back, rev);
        sid[e] = next_sid++;
        q.push_back(e);
    }
}

// Make transform (f1, f2)
int CubicSolver::create_transform(int f1, int f2){
    assert(f1 >= 0 && f2>= 0);
    assert(IS_REGULAR(f1) && IS_REGULAR(f2));
    assert(backlabel[f1] == -1 && backlabel[f2] == -1);
    int id = (int)h.size();
    h.push_back(h[f1]);
    for (__typeof(h[f2].begin()) it = h[f2].begin(); it != h[f2].end(); ++it) {
        if (h[id].count(*it))
            h[id].erase(*it);
        else
            h[id].insert(*it);
    }
    for (__typeof(h[id].begin()) it = h[id].begin(); it != h[id].end(); ++it)
        h[*it].insert(id);
    transform.push_back(make_pair(f1, f2));
    blossom.push_back(-1);
    mate.push_back(mate[f1]);
    path.push_back(vector<int>());
    backlabel.push_back(-1);
#ifdef DEG3DEBUGINFO
    revlabel.push_back(-1);
#endif
    sid.push_back(INF);
    return id;
}

// Initiate blossom
int CubicSolver::setup_blossom(){
    int bid = (int)blossoms.size();
    blossoms.push_back(vector<int>());
    tips.push_back(vector<int>());
    return bid;
}

void CubicSolver::add_vertex_to_blossom(int x, int bid, bool is_tip){
    if (blossom[x] >= 0) {
        if (blossom[x] != bid) {
            int bold = blossom[x];
            int label_count = 0;
            for (__typeof(blossoms[bold].begin()) it = blossoms[bold].begin(); it != blossoms[bold].end(); ++it) {
                blossom[*it] = bid;
                blossoms[bid].push_back(*it);
            }
            for (__typeof(tips[bold].begin()) it = tips[bold].begin(); it != tips[bold].end(); ++it) {
                if (is_tip) {
                    assert(backlabel[*it] == -1);
                    tips[bid].push_back(*it);
                } else {// old tip
                    if (backlabel[*it] == -1)
                        backlabel[*it] = -2;
                    else label_count++;
                }
            }
            assert(is_tip || label_count == 1);
            // blossoms[bold].clear();
        }
    } else {
        blossom[x] = bid;
        blossoms[bid].push_back(x);
        if (is_tip){
            assert(backlabel[x] == -1);
            tips[bid].push_back(x);
        }
    }
}

bool CubicSolver::iteration(){
    // Preprocess: if not taken with three different components, take and finish
    for (int v = 0; v < nv; ++v)
        if (!taken[v]) {
            int f1 = fau_find(nei[VJ2E(v, 0)]);
            int f2 = fau_find(nei[VJ2E(v, 1)]);
            int f3 = fau_find(nei[VJ2E(v, 2)]);
            if (f1 != f2 && f1 != f3 && f2 != f3){
                taken[v] = true;
                return true;
            }
        }
    // Initialization: singleton edges from dummy "root" to roots of every connected component
    comp_roots.clear();
    comp_id.clear(); comp_id.resize((unsigned)nu, -1);
    for (int u = 0; u < nu; ++u)
        if (fau_find(u) == u){
            comp_id[u] = (int)comp_roots.size();
            comp_roots.push_back(u);
        }
    // Initialization: blossoms are singletons, IDs are +infty, no backwards labels
    ne = 3*nv;
    nec = ne + (int)comp_roots.size();
    next_sid = 1;
    transform.clear(); tips.clear();
    blossom.clear(); blossom.resize((unsigned)nec, -1);
    backlabel.clear(); backlabel.resize((unsigned)nec, -1);
    sid.clear(); sid.resize((unsigned)nec, INF);
    path.clear(); path.resize((unsigned)nec);
    mate.clear(); mate.resize((unsigned)nec, -1);
    blossoms.clear();
#ifdef DEG3DEBUGINFO
    blossom_tips.clear(); blossom_starts.clear();
    revlabel.clear(); revlabel.resize((unsigned)nec, -1);
#endif
    // Initialization: make mates
    for (int v = 0; v < nv ; ++v){
        mate[VJ2E(v, 1)] = VJ2E(v, 2);
        mate[VJ2E(v, 2)] = VJ2E(v, 1);
    }
    for (int c = 0; c < (int)comp_roots.size(); ++c)
        mate[CC2E(c)] = CC2E(c);
    // Initialization: build dependency graph
    build_dependency_graph();
    // Initialization: Enqueue singular edges between root and connected components
    q.clear();
    for (int i = 0; i < (int)comp_roots.size(); ++i)
        label(CC2E(i), CC_LABEL_OFFSET-i, -1);
    // Main loop
    while(q.size() > 0){
        int e = q.front();
        q.pop_front();
        vector<pair<int, int> > enei;
        for(__typeof(h[e].begin()) it = h[e].begin(); it != h[e].end(); ++it)
            if (!same_blossom(e, *it))
                enei.push_back(make_pair(sid[*it], *it));
        sort(enei.begin(), enei.end());
        for (__typeof(enei.begin()) itf = enei.begin(); itf != enei.end(); ++itf){
            int f = itf->second;
            if (same_blossom(e, f)) continue;
            int mf = mate[f];
            assert(f >= 0 && mf >= 0);
            if (IS_LABELLED(f)){
                if (sid[f] < sid[e]){
                    set<int> f_blossoms, f_edges;
                    for (__typeof(path[f].rbegin()) it = path[f].rbegin(); it != path[f].rend(); ++it){
                        if (*it < 0)
                            break;
                        if (blossom[*it] >= 0)
                            f_blossoms.insert(blossom[*it]);
                        else
                            f_edges.insert(*it);
                    }
                    int b = -1;
                    for (__typeof(path[e].rbegin()) it = path[e].rbegin(); it != path[e].rend(); ++it)
                        if (f_edges.count(*it) > 0 || (*it >= 0 && f_blossoms.count(blossom[*it]) > 0)){
                            b = *it;
                            break;
                        }
                    if (b == -1){
                        // Augment
#ifdef DEG3DEBUGINFO
                        old_augments.clear();
                        old_taken = taken;
#endif
                        augment(path[e]);
                        augment(path[f]);
                        return true;
                    } else {
                        // Blossom
                        assert(b >= 0);
#ifdef DEG3DEBUGINFO
                        blossom_starts.push_back(make_pair(e, f));
#endif
                        new_blossom(e, f, b);
                    }
                }
            } else if (!IS_LABELLED(mf) && blossom[f] == -1){
                assert(backlabel[f] == -1 && backlabel[mf] == -1);
                if (h[e].count(mf) > 0){
                    // Degenerate blossom
#ifdef DEG3DEBUGINFO
                    blossom_starts.push_back(make_pair(f, mf));
                    blossom_tips.push_back(make_pair(f, mf));
#endif
                    int t = create_transform(f, mf);
                    int bid = setup_blossom();
                    add_vertex_to_blossom(f, bid, true);
                    add_vertex_to_blossom(mf, bid, true);
                    add_vertex_to_blossom(t, bid);
                    label(t, e, -1);
                } else {
                    // Grow step
                    label(mf, e, -1);
                }
            }
        }
    }
    return false;
}

bool CubicSolver::same_blossom(int e, int f){
    if (e == f)
        return true;
    if (e < 0 || f < 0 || blossom[e] == -1 || blossom[f] == -1)
        return false;
    return blossom[e] == blossom[f];
}

void CubicSolver::new_blossom(int e0, int e1, int b0){
    int e[2] = {e0, e1};
    int b[2] = {-1, -1};
    int t[2] = {e1, e0};
    for (int i = 0; i < 2; ++i){
        for (__typeof(path[e[i]].rbegin()) it = path[e[i]].rbegin(); it != path[e[i]].rend(); ++it)
            if (same_blossom(b0, *it)) {
                b[i] = *it;
                break;
            } else
                t[i] = *it;
    }
    assert(b[0] >= 0 && b[1] >=0 && t[0] >= 0 && t[1] >= 0);
    if (b[0] != b[1])
        t[0] = t[1] = -1;
    int bid = setup_blossom();
    vector<pair<int, pair<int, int> > > update_queue;
#ifdef DEG3DEBUGINFO
    blossom_tips.push_back(make_pair(t[0], t[1]));
#endif
    int last_mate = -1;
    for (int i = 0; i < 2; ++i)
        for (__typeof(path[e[i]].rbegin()) it = path[e[i]].rbegin(); it != path[e[i]].rend(); ++it) {
            if (*it == b[i])
                break;
            if (*it != t[i] && backlabel[*it] == -1) {
                assert(mate[last_mate] == *it);
                update_queue.push_back(make_pair(-(int) sid[last_mate], make_pair(*it, i)));
            } else
                last_mate = *it;
        }
    sort(update_queue.begin(), update_queue.end());
    for (__typeof(update_queue.begin()) it = update_queue.begin(); it != update_queue.end(); ++it)
        label(it->second.first, e[1-it->second.second], e[it->second.second]);
    for (int i = 0; i < 2; ++i)
        for (__typeof(path[e[i]].rbegin()) it = path[e[i]].rbegin(); it != path[e[i]].rend(); ++it) {
            if (*it == b[i])
                break;
            add_vertex_to_blossom(*it, bid, same_blossom(*it, t[0]) || same_blossom(*it, t[1]));
        }
    if (t[0] >= 0){
        int tf = create_transform(t[0], t[1]);
        add_vertex_to_blossom(tf, bid);
        label(tf, e[1], e[0]);
    } else {
        // Merge blossoms
        int bold = blossom[b[0]];
        assert(bold >= 0);
        add_vertex_to_blossom(b[0], bid, true);
    }
}

void CubicSolver::augment(vector<int> &path){
#ifdef DEG3DEBUGINFO
    old_augments.push_back(path);
#endif
    bool parity = true;
    for (__typeof(path.rbegin()) it = path.rbegin(); it != path.rend(); ++it) {
        if (parity) {
            if (IS_TRANSFORM(*it))
                taken[E2V(transform[E2TFM(*it)].first)] = !taken[E2V(transform[E2TFM(*it)].first)];
            else if (IS_REGULAR(*it))
                taken[E2V(*it)] = !taken[E2V(*it)];
        }
        parity = !parity;
    }
}

vector<int> CubicSolver::give_solution(){
    vector<int> res;
    for (int i = 0; i < nv; ++i)
        if (!taken[i])
            res.push_back(verts_names[i]);
    return res;
}

int CubicSolver::fau_find(int x){
    return fau[x] < 0 ? x : (fau[x] = fau_find(fau[x]));
}

bool CubicSolver::fau_join(int x, int y){
    x = fau_find(x);
    y = fau_find(y);
    if (x == y) return false;
    if (fau[x] > fau[y]) swap(x, y);
    else if (fau[x] == fau[y]) fau[x]--;
    fau[y] = x;
    return true;
}

/* Debug function. Recomputes fau, checking if the solution is feasible. */
void CubicSolver::check_solution(){
    fau.clear(); fau.resize(nu, -1);
    for (int i = 0; i < nv; ++i)
        if (taken[i])
            if (!fau_join(nei[VJ2E(i, 0)], nei[VJ2E(i, 1)]) || !fau_join(nei[VJ2E(i, 0)], nei[VJ2E(i, 2)]))
                throw new string("SOLUTION NOT FEASIBLE!");
}

void CubicSolver::dfs_undeletable(int x, int cmp_id, vector<int> &vis, const vector<int> &status, const vector<vector<int> > &neighbors) {
    vis[x] = cmp_id;
    for (__typeof(neighbors[x].begin())  it = neighbors[x].begin(); it != neighbors[x].end(); it++)
        if (status[*it] > 0 && vis[*it] == -1)
            dfs_undeletable(*it, cmp_id, vis, status, neighbors);
}

CubicSolver::CubicSolver(vector<int> &status, vector<vector<int> > &neighbors){
    nu = nv = 0;
    int n = (int)status.size();
    vector<int> vis(n, -1);
    for (int v = 0; v < n; ++v) {
        if (status[v] > 0) {
            if (vis[v] == -1)
                dfs_undeletable(v, nu++, vis, status, neighbors);
            names2ids[v] = vis[v];
        } else if (status[v] == 0){
            int w = nv++;
            verts_names.push_back(v);
            names2ids[v] = w;
        }
    }
    fau.resize((unsigned)nu, -1);
    nei.resize((unsigned)(3*nv));
    taken.resize((unsigned)nv, false);
    undels.resize((unsigned)nu);
    vector<int> pos((unsigned)nv, 0);
    for (int i = 0; i < nv; ++i){
        int v = verts_names[i];
        for (__typeof(neighbors[v].begin())  it = neighbors[v].begin(); it != neighbors[v].end(); it++){
            int w = *it;
            int j = names2ids[w];
            if (status[w] > 0){
                int u = VJ2E(i, pos[i]++);
                undels[j].push_back(u);
                nei[u] = j;
            } else if (status[w] == 0 && i < j){ // subdivide an edge between two deletables
                int k = nu++;
                undels.push_back(vector<int>());
                fau.push_back(-1);
                int ui = VJ2E(i, pos[i]++);
                int uj = VJ2E(j, pos[j]++);
                nei[ui] = k;
                nei[uj] = k;
                undels[k].push_back(ui);
                undels[k].push_back(uj);
            }
        }
    }
}

vector<int> CubicSolver::solve(){
    while(iteration())
        check_solution();
    check_solution();
    return give_solution();
}

