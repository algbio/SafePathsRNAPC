#include <algorithms/safe_paths.h>

#include <algorithms/greedy_approx.h>

#include <lemon/network_simplex.h>
#include <lemon/edmonds_karp.h>
#include <lemon/dfs.h>
#include <lemon/bfs.h>


using namespace lemon;



std::vector<ListDigraph::Arc> path_through_direct_edges(ListDigraph& red, ListDigraph::ArcMap<ListDigraph::Arc>& direct, ListDigraph::Node s, ListDigraph::Node t) {

    std::vector<std::pair<ListDigraph::Arc , ListDigraph::Node>> restorage_list;
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        ListDigraph::Arc d_e = direct[e];
        if (e != d_e) { // If it is a reverse edge
            restorage_list.push_back({e, red.target(e)});
            red.changeTarget(e, red.source(e));
        }
    }
    std::vector<ListDigraph::Arc> path;
    Bfs<ListDigraph> bfs_to_t(red);
    bfs_to_t.run(s, t);
    ListDigraph::Node temp_v = t;
    ListDigraph::Arc temp_e(INVALID);
    while ((temp_e = bfs_to_t.predArc(temp_v)) != INVALID) {
        path.push_back(temp_e);
        temp_v = red.source(temp_e);
    }


    for (auto& pair : restorage_list) {
        red.changeTarget(pair.first, pair.second);
    }

    return path;
}



std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_MPC(ListDigraph& g) {

    // Build the Min-Flow network reduction
    ListDigraph red;

    ListDigraph::NodeMap<ListDigraph::Node> v_in(g);
    ListDigraph::NodeMap<ListDigraph::Node> v_out(g);
    ListDigraph::NodeMap<ListDigraph::Node> original(red);
    ListDigraph::ArcMap<int> cost(red, 0);
    ListDigraph::ArcMap<int> demand(red, 0);
    ListDigraph::NodeMap<int> supply(red, 0);

    ListDigraph::Node s = red.addNode();
    ListDigraph::Node t = red.addNode();
    ListDigraph::Arc st = red.addArc(s, t);
    supply[s] = countNodes(g);
    supply[t] = -countNodes(g);


    // Set split vertices
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        v_in[v] = red.addNode();
        v_out[v] = red.addNode();
        original[v_in[v]] = v;
        original[v_out[v]] = v;

        ListDigraph::Arc split = red.addArc(v_in[v], v_out[v]);
        demand[split] = 1;

        ListDigraph::Arc sv = red.addArc(s, v_in[v]);
        cost[sv] = 1;

        red.addArc(v_out[v], t);
    }

    // Set edges connecting split vertices
    for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        red.addArc(v_out[u], v_in[v]);
    }



    // Use NetworkSimplex for solving the min-flow
    NetworkSimplex<ListDigraph> ns(red);
    ns.lowerMap(demand).costMap(cost).supplyMap(supply).run();
    int64_t width = ns.totalCost();


    // Obtain Flow solution
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ns.flowMap(flowMap);

    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;
    std::vector<std::vector<ListDigraph::Arc>> path_cover_edges_red;

    // Remove 0 flow edges and st

    // Stores the modified edges with the corresponding target
    // (we use the strategy to move the target to the source instead of removing)
    std::vector<std::pair<ListDigraph::Arc , ListDigraph::Node>> restorage_list = {{st, t}};
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        if (flowMap[e] == 0) {
            restorage_list.push_back({e, red.target(e)});
        }
    }
    for (auto& pair : restorage_list) {
        red.changeTarget(pair.first, red.source(pair.first));
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        std::vector<ListDigraph::Arc> path_edges_red;

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t) {
                if (path.empty() || path.back() != original[v]) {
                    path.push_back(original[v]);
                } else if (red.source(e) != s && (!path.empty() && path.back() == original[v])) {
                    path_edges_red.push_back(e);
                }
            }
            flowMap[e]--;
            if (flowMap[e] == 0) {
                restorage_list.push_back({e, red.target(e)});
                red.changeTarget(e, red.source(e));
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());
        std::reverse(path_edges_red.begin(), path_edges_red.end());

        path_cover.push_back(path);
        path_cover_edges_red.push_back(path_edges_red);


        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }



    // Compute Safe Paths

    // Restore removed edges in the process of computing the MPC
    for (auto& pair: restorage_list) {
        red.changeTarget(pair.first, pair.second);
    }

    std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_per_path;


    for (int i = 0; i < path_cover.size(); ++i) {
        std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>> path_maximal_safe_paths_pair;
        std::vector<std::vector<ListDigraph::Node>> path_maximal_safe_paths;

        std::vector<ListDigraph::Node>& path = path_cover[i];
        std::vector<ListDigraph::Arc>& path_edges_red = path_cover_edges_red[i];

        int x = 0, y = 0;
        bool fail_to_expand = false;
        while (y+1 < path.size()) {
            std::vector<ListDigraph::Arc> transitive_edges;
            ListDigraph::Node x_p = path[y+1];
            ListDigraph::Arc e = path_edges_red[y];


            // Compute reduction
            red.changeTarget(e, red.source(e));
            for (int z = x+1; z <= y; ++z) {
                ListDigraph::Node v = path[z];
                for (ListDigraph::InArcIt to_v(g, v); to_v != INVALID; ++to_v) {
                    ListDigraph::Node u = g.source(to_v);
                    if (u != path[z-1]) {
                        transitive_edges.push_back(red.addArc(v_out[u], v_in[x_p]));
                    }
                }
            }

            // Compute new width
            ns.reset();
            ns.lowerMap(demand).costMap(cost).supplyMap(supply).run();
            int64_t new_width = ns.totalCost(); // Maybe we have to check something else here for the RNA path cover case

            if (width == new_width) { // It is not safe
                // Report the path between x and y, move x to the right, and (if necessary) y to the right
                if (x != y && !fail_to_expand) {
                    std::vector<ListDigraph::Node> maximal_safe_path;
                    for (int z = x; z <= y; ++z) {
                        maximal_safe_path.push_back(path[z]);
                    }
                    path_maximal_safe_paths.push_back(maximal_safe_path);
                    fail_to_expand = true;
                }
                ++x;
                if (x > y) {
                    ++y;
                }
            } else { //Path is safe
                // Move y to the right
                ++y;
                fail_to_expand = false;
            }

            // Remove transitive edges and add e
            red.changeTarget(e, v_in[x_p]);
            for (ListDigraph::Arc e : transitive_edges) {
                red.erase(e);
            }
        }

        // (possibly) report the last path
        if (x != y) {
            std::vector<ListDigraph::Node> maximal_safe_path;
            for (int z = x; z <= y; ++z) {
                maximal_safe_path.push_back(path[z]);
            }
            path_maximal_safe_paths.push_back(maximal_safe_path);
        }

        path_maximal_safe_paths_pair.first = path;
        path_maximal_safe_paths_pair.second = path_maximal_safe_paths;
        path_maximal_safe_paths_per_path.push_back(path_maximal_safe_paths_pair);
    }

    return path_maximal_safe_paths_per_path;
}



std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> greedy_path_maximal_safe_paths_MPC(ListDigraph& g) {

    // Build the Min-Flow network reduction
    std::vector<std::vector<ListDigraph::Arc>> paths = greedy_approximation_MPC_edges(g);

    // Compute mu values according to the current path cover
    ListDigraph::ArcMap<int64_t> mu(g, 0); // Number of paths using this edge
    ListDigraph::NodeMap<int64_t> mu_v(g, 0); // Number of paths using this vertex
    ListDigraph::NodeMap<int64_t> starting_at(g, 0); // Number of paths starting at this vertex
    ListDigraph::NodeMap<int64_t> ending_at(g, 0); // Number of paths ending at this vertex
    for (auto path: paths) {
        for (int i = 0; i < path.size(); ++i) {
            auto edge = path[i];
            mu[edge]++;
            mu_v[g.source(edge)]++;
            if (i == path.size()-1) {
                mu_v[g.target(edge)]++;
                ending_at[g.target(edge)]++;
            }
            if (i == 0) {
                starting_at[g.source(edge)]++;
            }
        }
    }
    // Put flow in the corresponding vertices covered by a path of length one (does not appear in the answer)
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        if (mu_v[v] == 0) {
            mu_v[v]++;
            starting_at[v]++;
            ending_at[v]++;
        }
    }

    ListDigraph red;

    ListDigraph::NodeMap<ListDigraph::Arc> split_edges(g);
    ListDigraph::NodeMap<ListDigraph::Arc> from_s(g);
    ListDigraph::NodeMap<ListDigraph::Arc> to_t(g);
    ListDigraph::NodeMap<ListDigraph::Node> original(red);


    // Capacities of the MaxFlow reduction
    ListDigraph::ArcMap<int64_t> feasible_flow(red, 0); // From the approximation
    ListDigraph::ArcMap<int64_t> capacities(red, 0); // For the Max-Flow reduction, it is flow-demand


    ListDigraph::Node s = red.addNode();
    ListDigraph::Node t = red.addNode();


    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        ListDigraph::Node v_in = red.addNode();
        ListDigraph::Node v_out = red.addNode();
        original[v_in] = v;
        original[v_out] = v;

        ListDigraph::Arc split = red.addArc(v_in, v_out);
        split_edges[v] = split;
        feasible_flow[split] = mu_v[v];
        capacities[split] = mu_v[v] - 1;

        ListDigraph::Arc sv = red.addArc(s, v_in);
        ListDigraph::Arc vt = red.addArc(v_out, t);

        from_s[v] = sv;
        to_t[v] = vt;

        feasible_flow[sv] = starting_at[v];
        feasible_flow[vt] = ending_at[v];

        capacities[sv] = starting_at[v];
        capacities[vt] = ending_at[v];
    }

    for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        ListDigraph::Arc red_e = red.addArc(red.target(split_edges[u]), red.source(split_edges[v]));
        feasible_flow[red_e] = mu[e];
        capacities[red_e] = mu[e];
    }

    ListDigraph::ArcMap<ListDigraph::Arc> direct(red);
    ListDigraph::ArcMap<ListDigraph::Arc> reverse(red);
    std::vector<ListDigraph::Arc> edges;
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        edges.push_back(e);
    }

    for (auto e : edges) {
        ListDigraph::Arc rev_e = red.addArc(red.target(e), red.source(e));
        capacities[rev_e] = countNodes(g);
        direct[e] = e;
        direct[rev_e] = e;
        reverse[e] = rev_e;
        reverse[rev_e] = rev_e;
    }


    // Run Max-Flow algorithm
    EdmondsKarp<ListDigraph, ListDigraph::ArcMap<int64_t>> ek(red, capacities, s, t);

    // Set the flowMap to store the result in run
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ek.flowMap(flowMap);
    ek.run();


    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;
    std::vector<std::vector<ListDigraph::Arc>> path_cover_edges_red;

    // Remove 0 flow edges and st

    // Stores the modified edges with the corresponding target
    // (we use the strategy to move the target to the source instead of removing)
    std::vector<std::pair<ListDigraph::Arc , ListDigraph::Node>> restorage_list;
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        ListDigraph::Arc d_e = direct[e];
        if (e != d_e) { // If it is a reverse edge
            // The flow on that edge is computed as the previous flow minus the one discounted by that edge,
            // plus the flow in the reverse direction (discounted in the Max-flow, therefore pushed in the Min-flow)
            flowMap[d_e] = feasible_flow[d_e] + flowMap[e] - flowMap[d_e];
            flowMap[e] = 0;

            restorage_list.push_back({e, red.target(e)});
            if (flowMap[d_e] == 0) {
                restorage_list.push_back({d_e, red.target(d_e)});
            }

            capacities[e] = countNodes(g);
            capacities[d_e] = flowMap[d_e];
        }
    }
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        ListDigraph::Arc split = split_edges[v];
        capacities[split]--;
    }
    for (auto& pair : restorage_list) {
        red.changeTarget(pair.first, red.source(pair.first));
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        std::vector<ListDigraph::Arc> path_edges_red;

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t) {
                if (path.empty() || path.back() != original[v]) {
                    path.push_back(original[v]);
                } else if (red.source(e) != s && (!path.empty() && path.back() == original[v])) {
                    path_edges_red.push_back(e);
                }
            }
            flowMap[e]--;
            if (flowMap[e] == 0) {
                restorage_list.push_back({e, red.target(e)});
                red.changeTarget(e, red.source(e));
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());
        std::reverse(path_edges_red.begin(), path_edges_red.end());

        path_cover.push_back(path);
        path_cover_edges_red.push_back(path_edges_red);

        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }

    int64_t width = path_cover.size();


    // Compute Safe Paths

    // Restore removed edges in the process of computing the MPC
    for (auto& pair: restorage_list) {
        red.changeTarget(pair.first, pair.second);
    }

    std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_per_path;


    for (int i = 0; i < path_cover.size(); ++i) {
        std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>> path_maximal_safe_paths_pair;
        std::vector<std::vector<ListDigraph::Node>> path_maximal_safe_paths;

        std::vector<ListDigraph::Node>& path = path_cover[i];
        std::vector<ListDigraph::Arc>& path_edges_red = path_cover_edges_red[i];

        int x = 0, y = 0;
        bool fail_to_expand = false;
        while (y+1 < path.size()) {
            std::vector<ListDigraph::Arc> transitive_edges;
            ListDigraph::Node x_p_1 = path[y];
            ListDigraph::Node x_p = path[y+1];
            ListDigraph::Arc e = path_edges_red[y];
            int64_t mu_e = capacities[e];
            ListDigraph::Arc rev_e = reverse[e];


            // Compute reduction
            red.changeTarget(e, red.source(e));
            red.changeTarget(rev_e, red.source(rev_e));
            capacities[from_s[x_p]] += mu_e;
            capacities[to_t[x_p_1]] += mu_e;
            for (int z = x+1; z <= y; ++z) {
                ListDigraph::Node v = path[z];
                for (ListDigraph::InArcIt to_v(g, v); to_v != INVALID; ++to_v) {
                    ListDigraph::Node u = g.source(to_v);
                    if (u != path[z-1]) {
                        ListDigraph::Arc tran_e = red.addArc(red.target(split_edges[u]), red.source(split_edges[x_p]));
                        capacities[tran_e] = 0;
                        ListDigraph::Arc rev_tran_e = red.addArc(red.target(tran_e), red.source(tran_e));
                        capacities[rev_tran_e] = countNodes(g);

                        transitive_edges.push_back(tran_e);
                        transitive_edges.push_back(rev_tran_e);
                    }
                }
            }

            // Compute new width
            EdmondsKarp<ListDigraph, ListDigraph::ArcMap<int64_t>> ek(red, capacities, s, t);

            // Set the flowMap to store the result in run
            ListDigraph::ArcMap<int64_t> flowMap(red);
            ek.flowMap(flowMap);
            ek.run();

            int64_t new_width = width + mu_e - ek.flowValue();

            if (width == new_width) { // It is not safe
                // Report the path between x and y, move x to the right, and (if necessary) y to the right
                if (x != y && !fail_to_expand) {
                    std::vector<ListDigraph::Node> maximal_safe_path;
                    for (int z = x; z <= y; ++z) {
                        maximal_safe_path.push_back(path[z]);
                    }
                    path_maximal_safe_paths.push_back(maximal_safe_path);
                    fail_to_expand = true;
                }
                ++x;
                if (x > y) {
                    ++y;
                }
            } else { //Path is safe
                // Move y to the right
                ++y;
                fail_to_expand = false;
            }

            // Remove transitive edges and add e
            capacities[from_s[x_p]] -= mu_e;
            capacities[to_t[x_p_1]] -= mu_e;
            red.changeTarget(e, red.source(rev_e));
            red.changeTarget(rev_e, red.source(e));
            for (ListDigraph::Arc e : transitive_edges) {
                red.erase(e);
            }
        }

        // (possibly) report the last path
        if (x != y) {
            std::vector<ListDigraph::Node> maximal_safe_path;
            for (int z = x; z <= y; ++z) {
                maximal_safe_path.push_back(path[z]);
            }
            path_maximal_safe_paths.push_back(maximal_safe_path);
        }

        path_maximal_safe_paths_pair.first = path;
        path_maximal_safe_paths_pair.second = path_maximal_safe_paths;
        path_maximal_safe_paths_per_path.push_back(path_maximal_safe_paths_pair);
    }

    return path_maximal_safe_paths_per_path;
}



std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_MPC(ListDigraph& g, std::vector<ListDigraph::Node>& S, std::vector<ListDigraph::Node>& T) {

    ListDigraph::NodeMap<bool> in_S(g, false);
    for (ListDigraph::Node v : S) {
        in_S[v] = true;
    }

    // Build the Min-Flow network reduction
    ListDigraph red;

    ListDigraph::NodeMap<ListDigraph::Node> v_in(g);
    ListDigraph::NodeMap<ListDigraph::Node> v_out(g);
    ListDigraph::NodeMap<ListDigraph::Node> original(red);
    ListDigraph::ArcMap<int> cost(red, 0);
    ListDigraph::ArcMap<int> demand(red, 0);
    ListDigraph::NodeMap<int> supply(red, 0);

    ListDigraph::Node s = red.addNode();
    ListDigraph::Node t = red.addNode();
    ListDigraph::Arc st = red.addArc(s, t);
    supply[s] = countNodes(g);
    supply[t] = -countNodes(g);


    // Set split vertices
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        v_in[v] = red.addNode();
        v_out[v] = red.addNode();
        original[v_in[v]] = v;
        original[v_out[v]] = v;

        ListDigraph::Arc split = red.addArc(v_in[v], v_out[v]);
        demand[split] = 1;
    }
    for (ListDigraph::Node v : S) {
        ListDigraph::Arc sv = red.addArc(s, v_in[v]);
        cost[sv] = 1;
    }
    for (ListDigraph::Node v : T) {
        red.addArc(v_out[v], t);
    }

    // Set edges connecting split vertices
    for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        red.addArc(v_out[u], v_in[v]);
    }



    // Use NetworkSimplex for solving the min-flow
    NetworkSimplex<ListDigraph> ns(red);
    ns.lowerMap(demand).costMap(cost).supplyMap(supply).run();
    int64_t width = ns.totalCost();


    // Obtain Flow solution
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ns.flowMap(flowMap);

    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;
    std::vector<std::vector<ListDigraph::Arc>> path_cover_edges_red;

    // Remove 0 flow edges and st

    // Stores the modified edges with the corresponding target
    // (we use the strategy to move the target to the source instead of removing)
    std::vector<std::pair<ListDigraph::Arc , ListDigraph::Node>> restorage_list = {{st, t}};
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        if (flowMap[e] == 0) {
            restorage_list.push_back({e, red.target(e)});
        }
    }
    for (auto& pair : restorage_list) {
        red.changeTarget(pair.first, red.source(pair.first));
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        std::vector<ListDigraph::Arc> path_edges_red;

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t) {
                if (path.empty() || path.back() != original[v]) {
                    path.push_back(original[v]);
                } else if (red.source(e) != s && (!path.empty() && path.back() == original[v])) {
                    path_edges_red.push_back(e);
                }
            }
            flowMap[e]--;
            if (flowMap[e] == 0) {
                restorage_list.push_back({e, red.target(e)});
                red.changeTarget(e, red.source(e));
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());
        std::reverse(path_edges_red.begin(), path_edges_red.end());

        path_cover.push_back(path);
        path_cover_edges_red.push_back(path_edges_red);


        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }



    // Compute Safe Paths

    // Restore removed edges
    for (auto& pair: restorage_list) {
        red.changeTarget(pair.first, pair.second);
    }

    std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_per_path;


    for (int i = 0; i < path_cover.size(); ++i) {
        std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>> path_maximal_safe_paths_pair;
        std::vector<std::vector<ListDigraph::Node>> path_maximal_safe_paths;

        std::vector<ListDigraph::Node>& path = path_cover[i];
        std::vector<ListDigraph::Arc>& path_edges_red = path_cover_edges_red[i];

        int x = 0, y = 0;
        bool fail_to_expand = false;
        while (y+1 < path.size()) {
            std::vector<ListDigraph::Arc> transitive_edges;
            ListDigraph::Node x_p = path[y+1];
            ListDigraph::Arc e = path_edges_red[y];


            // Compute reduction
            red.changeTarget(e, red.source(e));
            for (int z = x+1; z <= y; ++z) {
                ListDigraph::Node v = path[z];
                for (ListDigraph::InArcIt to_v(g, v); to_v != INVALID; ++to_v) {
                    ListDigraph::Node u = g.source(to_v);
                    if (u != path[z-1]) {
                        transitive_edges.push_back(red.addArc(v_out[u], v_in[x_p]));
                    }
                } if (in_S[v]) {
                    transitive_edges.push_back(red.addArc(s, v_in[x_p]));
                }
            }

            // Compute new width
            ns.reset();
            NetworkSimplex<ListDigraph>::ProblemType result = ns.lowerMap(demand).costMap(cost).supplyMap(supply).run();
            int64_t new_width = ns.totalCost(); // Maybe we have to check something else here for the RNA path cover case

            if (width == new_width && result == NetworkSimplex<ListDigraph>::OPTIMAL) { // It is not safe
                // Report the path between x and y, move x to the right, and (if necessary) y to the right
                if (x != y && !fail_to_expand) {
                    std::vector<ListDigraph::Node> maximal_safe_path;
                    for (int z = x; z <= y; ++z) {
                        maximal_safe_path.push_back(path[z]);
                    }
                    path_maximal_safe_paths.push_back(maximal_safe_path);
                    fail_to_expand = true;
                }
                ++x;
                if (x > y) {
                    ++y;
                }
            } else { //Path is safe
                // Move y to the right
                ++y;
                fail_to_expand = false;
            }

            // Remove transitive edges and add e
            red.changeTarget(e, v_in[x_p]);
            for (ListDigraph::Arc e : transitive_edges) {
                red.erase(e);
            }
        }

        // (possibly) report the last path
        if (x != y) {
            std::vector<ListDigraph::Node> maximal_safe_path;
            for (int z = x; z <= y; ++z) {
                maximal_safe_path.push_back(path[z]);
            }
            path_maximal_safe_paths.push_back(maximal_safe_path);
        }

        path_maximal_safe_paths_pair.first = path;
        path_maximal_safe_paths_pair.second = path_maximal_safe_paths;
        path_maximal_safe_paths_per_path.push_back(path_maximal_safe_paths_pair);
    }

    return path_maximal_safe_paths_per_path;
}



std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> greedy_path_maximal_safe_paths_MPC(ListDigraph& g, std::vector<ListDigraph::Node>& S, std::vector<ListDigraph::Node>& T) {

    ListDigraph::NodeMap<bool> in_S(g, false);
    for (ListDigraph::Node v : S) {
        in_S[v] = true;
    }

    // Build the Min-Flow network reduction
    std::vector<std::vector<ListDigraph::Arc>> paths = greedy_approximation_MPC_edges(g, S, T);

    // Compute mu values according to the current path cover
    ListDigraph::ArcMap<int64_t> mu(g, 0); // Number of paths using this edge
    ListDigraph::NodeMap<int64_t> mu_v(g, 0); // Number of paths using this vertex
    ListDigraph::NodeMap<int64_t> starting_at(g, 0); // Number of paths starting at this vertex
    ListDigraph::NodeMap<int64_t> ending_at(g, 0); // Number of paths ending at this vertex
    for (auto path: paths) {
        for (int i = 0; i < path.size(); ++i) {
            auto edge = path[i];
            mu[edge]++;
            mu_v[g.source(edge)]++;
            if (i == path.size()-1) {
                mu_v[g.target(edge)]++;
                ending_at[g.target(edge)]++;
            }
            if (i == 0) {
                starting_at[g.source(edge)]++;
            }
        }
    }
    // Put flow in the corresponding vertices covered by a path of length one (does not appear in the answer)
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        if (mu_v[v] == 0) {
            mu_v[v]++;
            starting_at[v]++;
            ending_at[v]++;
        }
    }

    ListDigraph red;

    ListDigraph::NodeMap<ListDigraph::Arc> split_edges(g);
    ListDigraph::NodeMap<ListDigraph::Node> original(red);


    // Capacities of the MaxFlow reduction
    ListDigraph::ArcMap<int64_t> feasible_flow(red, 0); // From the approximation
    ListDigraph::ArcMap<int64_t> capacities(red, 0); // For the Max-Flow reduction, it is flow-demand


    ListDigraph::Node s = red.addNode();
    ListDigraph::Node t = red.addNode();


    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        ListDigraph::Node v_in = red.addNode();
        ListDigraph::Node v_out = red.addNode();
        original[v_in] = v;
        original[v_out] = v;

        ListDigraph::Arc split = red.addArc(v_in, v_out);
        split_edges[v] = split;
        feasible_flow[split] = mu_v[v];
        capacities[split] = mu_v[v] - 1;
    }
    for (ListDigraph::Node v : S) {
        ListDigraph::Arc sv = red.addArc(s, red.source(split_edges[v]));
        feasible_flow[sv] = starting_at[v];
        capacities[sv] = starting_at[v];
    }
    for (ListDigraph::Node v : T) {
        ListDigraph::Arc vt = red.addArc(red.target(split_edges[v]), t);
        feasible_flow[vt] = ending_at[v];
        capacities[vt] = ending_at[v];
    }

    for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        ListDigraph::Arc red_e = red.addArc(red.target(split_edges[u]), red.source(split_edges[v]));
        feasible_flow[red_e] = mu[e];
        capacities[red_e] = mu[e];
    }

    ListDigraph::ArcMap<ListDigraph::Arc> direct(red);
    ListDigraph::ArcMap<ListDigraph::Arc> reverse(red);
    std::vector<ListDigraph::Arc> edges;
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        edges.push_back(e);
    }

    for (auto e : edges) {
        ListDigraph::Arc rev_e = red.addArc(red.target(e), red.source(e));
        capacities[rev_e] = countNodes(g);
        direct[e] = e;
        direct[rev_e] = e;
        reverse[e] = rev_e;
        reverse[rev_e] = rev_e;
    }


    // Run Max-Flow algorithm
    EdmondsKarp<ListDigraph, ListDigraph::ArcMap<int64_t>> ek(red, capacities, s, t);

    // Set the flowMap to store the result in run
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ek.flowMap(flowMap);
    ek.run();


    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;
    std::vector<std::vector<ListDigraph::Arc>> path_cover_edges_red;

    // Remove 0 flow edges and st

    // Stores the modified edges with the corresponding target
    // (we use the strategy to move the target to the source instead of removing)
    std::vector<std::pair<ListDigraph::Arc , ListDigraph::Node>> restorage_list;
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        ListDigraph::Arc d_e = direct[e];
        if (e != d_e) { // If it is a reverse edge
            // The flow on that edge is computed as the previous flow minus the one discounted by that edge,
            // plus the flow in the reverse direction (discounted in the Max-flow, therefore pushed in the Min-flow)
            flowMap[d_e] = feasible_flow[d_e] + flowMap[e] - flowMap[d_e];
            flowMap[e] = 0;

            restorage_list.push_back({e, red.target(e)});
            if (flowMap[d_e] == 0) {
                restorage_list.push_back({d_e, red.target(d_e)});
            }

            capacities[e] = countNodes(g);
            capacities[d_e] = flowMap[d_e];
        }
    }
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        ListDigraph::Arc split = split_edges[v];
        capacities[split]--;
    }
    for (auto& pair : restorage_list) {
        red.changeTarget(pair.first, red.source(pair.first));
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        std::vector<ListDigraph::Arc> path_edges_red;

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t) {
                if (path.empty() || path.back() != original[v]) {
                    path.push_back(original[v]);
                } else if (red.source(e) != s && (!path.empty() && path.back() == original[v])) {
                    path_edges_red.push_back(e);
                }
            }
            flowMap[e]--;
            if (flowMap[e] == 0) {
                restorage_list.push_back({e, red.target(e)});
                red.changeTarget(e, red.source(e));
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());
        std::reverse(path_edges_red.begin(), path_edges_red.end());

        path_cover.push_back(path);
        path_cover_edges_red.push_back(path_edges_red);

        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }

    int64_t width = path_cover.size();


    // Compute Safe Paths

    // Restore removed edges in the process of computing the MPC
    for (auto& pair: restorage_list) {
        red.changeTarget(pair.first, pair.second);
    }

    std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_per_path;


    for (int i = 0; i < path_cover.size(); ++i) {
        std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>> path_maximal_safe_paths_pair;
        std::vector<std::vector<ListDigraph::Node>> path_maximal_safe_paths;

        std::vector<ListDigraph::Node>& path = path_cover[i];
        std::vector<ListDigraph::Arc>& path_edges_red = path_cover_edges_red[i];

        int x = 0, y = 0;
        bool fail_to_expand = false;
        while (y+1 < path.size()) {
            std::vector<ListDigraph::Arc> transitive_edges;
            ListDigraph::Node x_p_1 = path[y];
            ListDigraph::Node x_p = path[y+1];
            ListDigraph::Arc e = path_edges_red[y];
            int64_t mu_e = capacities[e];
            ListDigraph::Arc rev_e = reverse[e];


            // Compute reduction

            red.changeTarget(e, red.source(e));
            red.changeTarget(rev_e, red.source(rev_e));
            for (int z = x+1; z <= y; ++z) {
                ListDigraph::Node v = path[z];
                for (ListDigraph::InArcIt to_v(g, v); to_v != INVALID; ++to_v) {
                    ListDigraph::Node u = g.source(to_v);
                    if (u != path[z-1]) {
                        ListDigraph::Arc tran_e = red.addArc(red.target(split_edges[u]), red.source(split_edges[x_p]));
                        direct[tran_e] = tran_e;
                        capacities[tran_e] = 0;

                        ListDigraph::Arc rev_tran_e = red.addArc(red.target(tran_e), red.source(tran_e));
                        direct[rev_tran_e] = tran_e;
                        capacities[rev_tran_e] = countNodes(g);

                        transitive_edges.push_back(tran_e);
                        transitive_edges.push_back(rev_tran_e);
                    }
                }
                if (in_S[v]) {
                    ListDigraph::Arc tran_e = red.addArc(s, red.source(split_edges[x_p]));
                    direct[tran_e] = tran_e;
                    capacities[tran_e] = 0;

                    ListDigraph::Arc rev_tran_e = red.addArc(red.target(tran_e), red.source(tran_e));
                    direct[rev_tran_e] = tran_e;
                    capacities[rev_tran_e] = countNodes(g);

                    transitive_edges.push_back(tran_e);
                    transitive_edges.push_back(rev_tran_e);
                }
            }

            // Find path from s to x_p in red and push mu_e flow through it
            std::vector<ListDigraph::Arc> s_x_p_path = path_through_direct_edges(red, direct, s, red.source(rev_e));
            for (ListDigraph::Arc temp_e: s_x_p_path) {
                capacities[temp_e] += mu_e;
            }
            // Find path from x_p_1 to t in red and push mu_e flow through it
            std::vector<ListDigraph::Arc> x_p_1_t_path = path_through_direct_edges(red, direct, red.source(e), t);
            for (ListDigraph::Arc temp_e: x_p_1_t_path) {
                capacities[temp_e] += mu_e;
            }


            if (s_x_p_path.size() == 0 || x_p_1_t_path.size() == 0) {
                //Path is safe
                // Move y to the right
                ++y;
                fail_to_expand = false;
            } else {
                // Compute new width
                EdmondsKarp<ListDigraph, ListDigraph::ArcMap<int64_t>> ek(red, capacities, s, t);

                // Set the flowMap to store the result in run
                ListDigraph::ArcMap<int64_t> flowMap(red);
                ek.flowMap(flowMap);
                ek.run();

                int64_t new_width = width + mu_e - ek.flowValue();

                if (width == new_width) { // It is not safe
                    // Report the path between x and y, move x to the right, and (if necessary) y to the right
                    if (x != y && !fail_to_expand) {
                        std::vector<ListDigraph::Node> maximal_safe_path;
                        for (int z = x; z <= y; ++z) {
                            maximal_safe_path.push_back(path[z]);
                        }
                        path_maximal_safe_paths.push_back(maximal_safe_path);
                        fail_to_expand = true;
                    }
                    ++x;
                    if (x > y) {
                        ++y;
                    }
                } else { //Path is safe
                    // Move y to the right
                    ++y;
                    fail_to_expand = false;
                }
            }


            // Remove transitive edges and add e
            for (ListDigraph::Arc temp_e : s_x_p_path) {
                capacities[temp_e] -= mu_e;
            }
            for (ListDigraph::Arc temp_e : x_p_1_t_path) {
                capacities[temp_e] -= mu_e;
            }
            red.changeTarget(e, red.source(rev_e));
            red.changeTarget(rev_e, red.source(e));
            for (ListDigraph::Arc e : transitive_edges) {
                red.erase(e);
            }
        }

        // (possibly) report the last path
        if (x != y) {
            std::vector<ListDigraph::Node> maximal_safe_path;
            for (int z = x; z <= y; ++z) {
                maximal_safe_path.push_back(path[z]);
            }
            path_maximal_safe_paths.push_back(maximal_safe_path);
        }

        path_maximal_safe_paths_pair.first = path;
        path_maximal_safe_paths_pair.second = path_maximal_safe_paths;
        path_maximal_safe_paths_per_path.push_back(path_maximal_safe_paths_pair);
    }

    return path_maximal_safe_paths_per_path;
}



std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_U_MPC(ListDigraph& g, std::vector<lemon::ListDigraph::Node>& U) {

    // Compute in_U
    ListDigraph::NodeMap<bool> in_U(g, false);
    for (ListDigraph::Node v : U) {
        in_U[v] = true;
    }

    // Build the Min-Flow network reduction
    ListDigraph red;

    ListDigraph::NodeMap<ListDigraph::Node> v_in(g);
    ListDigraph::NodeMap<ListDigraph::Node> v_out(g);
    ListDigraph::NodeMap<ListDigraph::Node> original(red);
    ListDigraph::ArcMap<int> cost(red, 0);
    ListDigraph::ArcMap<int> demand(red, 0);
    ListDigraph::NodeMap<int> supply(red, 0);

    ListDigraph::Node s = red.addNode();
    ListDigraph::Node t = red.addNode();
    ListDigraph::Arc st = red.addArc(s, t);
    supply[s] = countNodes(g);
    supply[t] = -countNodes(g);


    // Set split vertices
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        v_in[v] = red.addNode();
        v_out[v] = red.addNode();
        original[v_in[v]] = v;
        original[v_out[v]] = v;

        ListDigraph::Arc split = red.addArc(v_in[v], v_out[v]);
        if (in_U[v]) {
            demand[split] = 1;
        }

        ListDigraph::Arc sv = red.addArc(s, v_in[v]);
        cost[sv] = 1;

        red.addArc(v_out[v], t);
    }

    // Set edges connecting split vertices
    for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        red.addArc(v_out[u], v_in[v]);
    }



    // Use NetworkSimplex for solving the min-flow
    NetworkSimplex<ListDigraph> ns(red);
    ns.lowerMap(demand).costMap(cost).supplyMap(supply).run();
    int64_t width = ns.totalCost();


    // Obtain Flow solution
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ns.flowMap(flowMap);

    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;
    std::vector<std::vector<ListDigraph::Arc>> path_cover_edges_red;

    // Remove 0 flow edges and st

    // Stores the modified edges with the corresponding target
    // (we use the strategy to move the target to the source instead of removing)
    std::vector<std::pair<ListDigraph::Arc , ListDigraph::Node>> restorage_list = {{st, t}};
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        if (flowMap[e] == 0) {
            restorage_list.push_back({e, red.target(e)});
        }
    }
    for (auto& pair : restorage_list) {
        red.changeTarget(pair.first, red.source(pair.first));
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        std::vector<ListDigraph::Arc> path_edges_red;

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t) {
                if (path.empty() || path.back() != original[v]) {
                    path.push_back(original[v]);
                } else if (red.source(e) != s && (!path.empty() && path.back() == original[v])) {
                    path_edges_red.push_back(e);
                }
            }
            flowMap[e]--;
            if (flowMap[e] == 0) {
                restorage_list.push_back({e, red.target(e)});
                red.changeTarget(e, red.source(e));
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());
        std::reverse(path_edges_red.begin(), path_edges_red.end());

        path_cover.push_back(path);
        path_cover_edges_red.push_back(path_edges_red);


        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }



    // Compute Safe Paths

    // Restore removed edges in the process of computing the MPC
    for (auto& pair: restorage_list) {
        red.changeTarget(pair.first, pair.second);
    }

    std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_per_path;


    for (int i = 0; i < path_cover.size(); ++i) {
        std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>> path_maximal_safe_paths_pair;
        std::vector<std::vector<ListDigraph::Node>> path_maximal_safe_paths;

        std::vector<ListDigraph::Node>& path = path_cover[i];
        std::vector<ListDigraph::Arc>& path_edges_red = path_cover_edges_red[i];

        int x = 0, y = 0;
        bool fail_to_expand = false;
        while (y+1 < path.size()) {
            std::vector<ListDigraph::Arc> transitive_edges;
            ListDigraph::Node x_p = path[y+1];
            ListDigraph::Arc e = path_edges_red[y];


            // Compute reduction
            red.changeTarget(e, red.source(e));
            for (int z = x+1; z <= y; ++z) {
                ListDigraph::Node v = path[z];
                for (ListDigraph::InArcIt to_v(g, v); to_v != INVALID; ++to_v) {
                    ListDigraph::Node u = g.source(to_v);
                    if (u != path[z-1]) {
                        transitive_edges.push_back(red.addArc(v_out[u], v_in[x_p]));
                    }
                }
            }

            // Compute new width
            ns.reset();
            ns.lowerMap(demand).costMap(cost).supplyMap(supply).run();
            int64_t new_width = ns.totalCost(); // Maybe we have to check something else here for the RNA path cover case

            if (width == new_width) { // It is not safe
                // Report the path between x and y, move x to the right, and (if necessary) y to the right
                if (x != y && !fail_to_expand) {
                    std::vector<ListDigraph::Node> maximal_safe_path;
                    for (int z = x; z <= y; ++z) {
                        maximal_safe_path.push_back(path[z]);
                    }
                    path_maximal_safe_paths.push_back(maximal_safe_path);
                    fail_to_expand = true;
                }
                ++x;
                if (x > y) {
                    ++y;
                }
            } else { //Path is safe
                // Move y to the right
                ++y;
                fail_to_expand = false;
            }

            // Remove transitive edges and add e
            red.changeTarget(e, v_in[x_p]);
            for (ListDigraph::Arc e : transitive_edges) {
                red.erase(e);
            }
        }

        // (possibly) report the last path
        if (x != y) {
            std::vector<ListDigraph::Node> maximal_safe_path;
            for (int z = x; z <= y; ++z) {
                maximal_safe_path.push_back(path[z]);
            }
            path_maximal_safe_paths.push_back(maximal_safe_path);
        }

        path_maximal_safe_paths_pair.first = path;
        path_maximal_safe_paths_pair.second = path_maximal_safe_paths;
        path_maximal_safe_paths_per_path.push_back(path_maximal_safe_paths_pair);
    }

    return path_maximal_safe_paths_per_path;
}



std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> greedy_path_maximal_safe_paths_U_MPC(ListDigraph& g, std::vector<lemon::ListDigraph::Node>& U) {

    // Compute in_U
    ListDigraph::NodeMap<bool> in_U(g, false);
    for (ListDigraph::Node v : U) {
        in_U[v] = true;
    }

    // Build the Min-Flow network reduction
    std::vector<std::vector<ListDigraph::Arc>> paths = greedy_approximation_U_MPC_edges(g, U);

    // Compute mu values according to the current path cover
    ListDigraph::ArcMap<int64_t> mu(g, 0); // Number of paths using this edge
    ListDigraph::NodeMap<int64_t> mu_v(g, 0); // Number of paths using this vertex
    ListDigraph::NodeMap<int64_t> starting_at(g, 0); // Number of paths starting at this vertex
    ListDigraph::NodeMap<int64_t> ending_at(g, 0); // Number of paths ending at this vertex
    for (auto path: paths) {
        for (int i = 0; i < path.size(); ++i) {
            auto edge = path[i];
            mu[edge]++;
            mu_v[g.source(edge)]++;
            if (i == path.size()-1) {
                mu_v[g.target(edge)]++;
                ending_at[g.target(edge)]++;
            }
            if (i == 0) {
                starting_at[g.source(edge)]++;
            }
        }
    }
    // Put flow in the corresponding vertices covered by a path of length one (does not appear in the answer)
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        if (mu_v[v] == 0 && in_U[v]) {
            mu_v[v]++;
            starting_at[v]++;
            ending_at[v]++;
        }
    }

    ListDigraph red;

    ListDigraph::NodeMap<ListDigraph::Arc> split_edges(g);
    ListDigraph::NodeMap<ListDigraph::Arc> from_s(g);
    ListDigraph::NodeMap<ListDigraph::Arc> to_t(g);
    ListDigraph::NodeMap<ListDigraph::Node> original(red);


    // Capacities of the MaxFlow reduction
    ListDigraph::ArcMap<int64_t> feasible_flow(red, 0); // From the approximation
    ListDigraph::ArcMap<int64_t> capacities(red, 0); // For the Max-Flow reduction, it is flow-demand


    ListDigraph::Node s = red.addNode();
    ListDigraph::Node t = red.addNode();


    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        ListDigraph::Node v_in = red.addNode();
        ListDigraph::Node v_out = red.addNode();
        original[v_in] = v;
        original[v_out] = v;

        ListDigraph::Arc split = red.addArc(v_in, v_out);
        split_edges[v] = split;
        feasible_flow[split] = mu_v[v];
        capacities[split] = mu_v[v];
        if (in_U[v]) {
            capacities[split]--;
        }

        ListDigraph::Arc sv = red.addArc(s, v_in);
        ListDigraph::Arc vt = red.addArc(v_out, t);

        from_s[v] = sv;
        to_t[v] = vt;

        feasible_flow[sv] = starting_at[v];
        feasible_flow[vt] = ending_at[v];

        capacities[sv] = starting_at[v];
        capacities[vt] = ending_at[v];
    }

    for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        ListDigraph::Arc red_e = red.addArc(red.target(split_edges[u]), red.source(split_edges[v]));
        feasible_flow[red_e] = mu[e];
        capacities[red_e] = mu[e];
    }

    ListDigraph::ArcMap<ListDigraph::Arc> direct(red);
    ListDigraph::ArcMap<ListDigraph::Arc> reverse(red);
    std::vector<ListDigraph::Arc> edges;
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        edges.push_back(e);
    }

    for (auto e : edges) {
        ListDigraph::Arc rev_e = red.addArc(red.target(e), red.source(e));
        capacities[rev_e] = countNodes(g);
        direct[e] = e;
        direct[rev_e] = e;
        reverse[e] = rev_e;
        reverse[rev_e] = rev_e;
    }


    // Run Max-Flow algorithm
    EdmondsKarp<ListDigraph, ListDigraph::ArcMap<int64_t>> ek(red, capacities, s, t);

    // Set the flowMap to store the result in run
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ek.flowMap(flowMap);
    ek.run();


    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;
    std::vector<std::vector<ListDigraph::Arc>> path_cover_edges_red;

    // Remove 0 flow edges and st

    // Stores the modified edges with the corresponding target
    // (we use the strategy to move the target to the source instead of removing)
    std::vector<std::pair<ListDigraph::Arc , ListDigraph::Node>> restorage_list;
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        ListDigraph::Arc d_e = direct[e];
        if (e != d_e) { // If it is a reverse edge
            // The flow on that edge is computed as the previous flow minus the one discounted by that edge,
            // plus the flow in the reverse direction (discounted in the Max-flow, therefore pushed in the Min-flow)
            flowMap[d_e] = feasible_flow[d_e] + flowMap[e] - flowMap[d_e];
            flowMap[e] = 0;

            restorage_list.push_back({e, red.target(e)});
            if (flowMap[d_e] == 0) {
                restorage_list.push_back({d_e, red.target(d_e)});
            }

            capacities[e] = countNodes(g);
            capacities[d_e] = flowMap[d_e];
        }
    }
    for (ListDigraph::Node v : U) {
        ListDigraph::Arc split = split_edges[v];
        capacities[split]--;
    }
    for (auto& pair : restorage_list) {
        red.changeTarget(pair.first, red.source(pair.first));
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        std::vector<ListDigraph::Arc> path_edges_red;

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t) {
                if (path.empty() || path.back() != original[v]) {
                    path.push_back(original[v]);
                } else if (red.source(e) != s && (!path.empty() && path.back() == original[v])) {
                    path_edges_red.push_back(e);
                }
            }
            flowMap[e]--;
            if (flowMap[e] == 0) {
                restorage_list.push_back({e, red.target(e)});
                red.changeTarget(e, red.source(e));
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());
        std::reverse(path_edges_red.begin(), path_edges_red.end());

        path_cover.push_back(path);
        path_cover_edges_red.push_back(path_edges_red);

        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }

    int64_t width = path_cover.size();


    // Compute Safe Paths

    // Restore removed edges in the process of computing the MPC
    for (auto& pair: restorage_list) {
        red.changeTarget(pair.first, pair.second);
    }

    std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_per_path;


    for (int i = 0; i < path_cover.size(); ++i) {
        std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>> path_maximal_safe_paths_pair;
        std::vector<std::vector<ListDigraph::Node>> path_maximal_safe_paths;

        std::vector<ListDigraph::Node>& path = path_cover[i];
        std::vector<ListDigraph::Arc>& path_edges_red = path_cover_edges_red[i];

        int x = 0, y = 0;
        bool fail_to_expand = false;
        while (y+1 < path.size()) {
            std::vector<ListDigraph::Arc> transitive_edges;
            ListDigraph::Node x_p_1 = path[y];
            ListDigraph::Node x_p = path[y+1];
            ListDigraph::Arc e = path_edges_red[y];
            int64_t mu_e = capacities[e];
            ListDigraph::Arc rev_e = reverse[e];


            // Compute reduction
            red.changeTarget(e, red.source(e));
            red.changeTarget(rev_e, red.source(rev_e));
            capacities[from_s[x_p]] += mu_e;
            capacities[to_t[x_p_1]] += mu_e;
            for (int z = x+1; z <= y; ++z) {
                ListDigraph::Node v = path[z];
                for (ListDigraph::InArcIt to_v(g, v); to_v != INVALID; ++to_v) {
                    ListDigraph::Node u = g.source(to_v);
                    if (u != path[z-1]) {
                        ListDigraph::Arc tran_e = red.addArc(red.target(split_edges[u]), red.source(split_edges[x_p]));
                        capacities[tran_e] = 0;
                        ListDigraph::Arc rev_tran_e = red.addArc(red.target(tran_e), red.source(tran_e));
                        capacities[rev_tran_e] = countNodes(g);

                        transitive_edges.push_back(tran_e);
                        transitive_edges.push_back(rev_tran_e);
                    }
                }
            }

            // Compute new width
            EdmondsKarp<ListDigraph, ListDigraph::ArcMap<int64_t>> ek(red, capacities, s, t);

            // Set the flowMap to store the result in run
            ListDigraph::ArcMap<int64_t> flowMap(red);
            ek.flowMap(flowMap);
            ek.run();

            int64_t new_width = width + mu_e - ek.flowValue();

            if (width == new_width) { // It is not safe
                // Report the path between x and y, move x to the right, and (if necessary) y to the right
                if (x != y && !fail_to_expand) {
                    std::vector<ListDigraph::Node> maximal_safe_path;
                    for (int z = x; z <= y; ++z) {
                        maximal_safe_path.push_back(path[z]);
                    }
                    path_maximal_safe_paths.push_back(maximal_safe_path);
                    fail_to_expand = true;
                }
                ++x;
                if (x > y) {
                    ++y;
                }
            } else { //Path is safe
                // Move y to the right
                ++y;
                fail_to_expand = false;
            }

            // Remove transitive edges and add e
            capacities[from_s[x_p]] -= mu_e;
            capacities[to_t[x_p_1]] -= mu_e;
            red.changeTarget(e, red.source(rev_e));
            red.changeTarget(rev_e, red.source(e));
            for (ListDigraph::Arc e : transitive_edges) {
                red.erase(e);
            }
        }

        // (possibly) report the last path
        if (x != y) {
            std::vector<ListDigraph::Node> maximal_safe_path;
            for (int z = x; z <= y; ++z) {
                maximal_safe_path.push_back(path[z]);
            }
            path_maximal_safe_paths.push_back(maximal_safe_path);
        }

        path_maximal_safe_paths_pair.first = path;
        path_maximal_safe_paths_pair.second = path_maximal_safe_paths;
        path_maximal_safe_paths_per_path.push_back(path_maximal_safe_paths_pair);
    }

    return path_maximal_safe_paths_per_path;
}



std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_U_MPC(ListDigraph& g, std::vector<ListDigraph::Node>& S, std::vector<ListDigraph::Node>& T, std::vector<lemon::ListDigraph::Node>& U) {

    // Compute in_U
    ListDigraph::NodeMap<bool> in_U(g, false);
    ListDigraph::NodeMap<bool> in_S(g, false);
    for (ListDigraph::Node v : U) {
        in_U[v] = true;
    }
    for (ListDigraph::Node v : S) {
        in_S[v] = true;
    }

    // Build the Min-Flow network reduction
    ListDigraph red;

    ListDigraph::NodeMap<ListDigraph::Node> v_in(g);
    ListDigraph::NodeMap<ListDigraph::Node> v_out(g);
    ListDigraph::NodeMap<ListDigraph::Node> original(red);
    ListDigraph::ArcMap<int> cost(red, 0);
    ListDigraph::ArcMap<int> demand(red, 0);
    ListDigraph::NodeMap<int> supply(red, 0);

    ListDigraph::Node s = red.addNode();
    ListDigraph::Node t = red.addNode();
    ListDigraph::Arc st = red.addArc(s, t);
    supply[s] = countNodes(g);
    supply[t] = -countNodes(g);


    // Set split vertices
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        v_in[v] = red.addNode();
        v_out[v] = red.addNode();
        original[v_in[v]] = v;
        original[v_out[v]] = v;

        ListDigraph::Arc split = red.addArc(v_in[v], v_out[v]);
        if (in_U[v]) {
            demand[split] = 1;
        }
    }
    for (ListDigraph::Node v : S) {
        ListDigraph::Arc sv = red.addArc(s, v_in[v]);
        cost[sv] = 1;
    }
    for (ListDigraph::Node v : T) {
        red.addArc(v_out[v], t);
    }

    // Set edges connecting split vertices
    for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        red.addArc(v_out[u], v_in[v]);
    }



    // Use NetworkSimplex for solving the min-flow
    NetworkSimplex<ListDigraph> ns(red);
    ns.lowerMap(demand).costMap(cost).supplyMap(supply).run();
    int64_t width = ns.totalCost();


    // Obtain Flow solution
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ns.flowMap(flowMap);

    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;
    std::vector<std::vector<ListDigraph::Arc>> path_cover_edges_red;

    // Remove 0 flow edges and st

    // Stores the modified edges with the corresponding target
    // (we use the strategy to move the target to the source instead of removing)
    std::vector<std::pair<ListDigraph::Arc , ListDigraph::Node>> restorage_list = {{st, t}};
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        if (flowMap[e] == 0) {
            restorage_list.push_back({e, red.target(e)});
        }
    }
    for (auto& pair : restorage_list) {
        red.changeTarget(pair.first, red.source(pair.first));
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        std::vector<ListDigraph::Arc> path_edges_red;

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t) {
                if (path.empty() || path.back() != original[v]) {
                    path.push_back(original[v]);
                } else if (red.source(e) != s && (!path.empty() && path.back() == original[v])) {
                    path_edges_red.push_back(e);
                }
            }
            flowMap[e]--;
            if (flowMap[e] == 0) {
                restorage_list.push_back({e, red.target(e)});
                red.changeTarget(e, red.source(e));
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());
        std::reverse(path_edges_red.begin(), path_edges_red.end());

        path_cover.push_back(path);
        path_cover_edges_red.push_back(path_edges_red);


        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }



    // Compute Safe Paths

    // Restore removed edges
    for (auto& pair: restorage_list) {
        red.changeTarget(pair.first, pair.second);
    }

    std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_per_path;


    for (int i = 0; i < path_cover.size(); ++i) {
        std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>> path_maximal_safe_paths_pair;
        std::vector<std::vector<ListDigraph::Node>> path_maximal_safe_paths;

        std::vector<ListDigraph::Node>& path = path_cover[i];
        std::vector<ListDigraph::Arc>& path_edges_red = path_cover_edges_red[i];

        int x = 0, y = 0;
        bool fail_to_expand = false;
        while (y+1 < path.size()) {
            std::vector<ListDigraph::Arc> transitive_edges;
            ListDigraph::Node x_p = path[y+1];
            ListDigraph::Arc e = path_edges_red[y];


            // Compute reduction
            red.changeTarget(e, red.source(e));
            for (int z = x+1; z <= y; ++z) {
                ListDigraph::Node v = path[z];
                for (ListDigraph::InArcIt to_v(g, v); to_v != INVALID; ++to_v) {
                    ListDigraph::Node u = g.source(to_v);
                    if (u != path[z-1]) {
                        transitive_edges.push_back(red.addArc(v_out[u], v_in[x_p]));
                    }
                } if (in_S[v]) {
                    transitive_edges.push_back(red.addArc(s, v_in[x_p]));
                }
            }

            // Compute new width
            ns.reset();
            NetworkSimplex<ListDigraph>::ProblemType result = ns.lowerMap(demand).costMap(cost).supplyMap(supply).run();
            int64_t new_width = ns.totalCost(); // Maybe we have to check something else here for the RNA path cover case

            if (width == new_width && result == NetworkSimplex<ListDigraph>::OPTIMAL) { // It is not safe
                // Report the path between x and y, move x to the right, and (if necessary) y to the right
                if (x != y && !fail_to_expand) {
                    std::vector<ListDigraph::Node> maximal_safe_path;
                    for (int z = x; z <= y; ++z) {
                        maximal_safe_path.push_back(path[z]);
                    }
                    path_maximal_safe_paths.push_back(maximal_safe_path);
                    fail_to_expand = true;
                }
                ++x;
                if (x > y) {
                    ++y;
                }
            } else { //Path is safe
                // Move y to the right
                ++y;
                fail_to_expand = false;
            }

            // Remove transitive edges and add e
            red.changeTarget(e, v_in[x_p]);
            for (ListDigraph::Arc e : transitive_edges) {
                red.erase(e);
            }
        }

        // (possibly) report the last path
        if (x != y) {
            std::vector<ListDigraph::Node> maximal_safe_path;
            for (int z = x; z <= y; ++z) {
                maximal_safe_path.push_back(path[z]);
            }
            path_maximal_safe_paths.push_back(maximal_safe_path);
        }

        path_maximal_safe_paths_pair.first = path;
        path_maximal_safe_paths_pair.second = path_maximal_safe_paths;
        path_maximal_safe_paths_per_path.push_back(path_maximal_safe_paths_pair);
    }

    return path_maximal_safe_paths_per_path;
}



std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> greedy_path_maximal_safe_paths_U_MPC(ListDigraph& g, std::vector<ListDigraph::Node>& S, std::vector<ListDigraph::Node>& T, std::vector<lemon::ListDigraph::Node>& U) {

    // Compute in_U
    ListDigraph::NodeMap<bool> in_U(g, false);
    ListDigraph::NodeMap<bool> in_S(g, false);
    for (ListDigraph::Node v : U) {
        in_U[v] = true;
    }
    for (ListDigraph::Node v : S) {
        in_S[v] = true;
    }

    // Build the Min-Flow network reduction
    std::vector<std::vector<ListDigraph::Arc>> paths = greedy_approximation_U_MPC_edges(g, S, T, U);

    // Compute mu values according to the current path cover
    ListDigraph::ArcMap<int64_t> mu(g, 0); // Number of paths using this edge
    ListDigraph::NodeMap<int64_t> mu_v(g, 0); // Number of paths using this vertex
    ListDigraph::NodeMap<int64_t> starting_at(g, 0); // Number of paths starting at this vertex
    ListDigraph::NodeMap<int64_t> ending_at(g, 0); // Number of paths ending at this vertex
    for (auto path: paths) {
        for (int i = 0; i < path.size(); ++i) {
            auto edge = path[i];
            mu[edge]++;
            mu_v[g.source(edge)]++;
            if (i == path.size()-1) {
                mu_v[g.target(edge)]++;
                ending_at[g.target(edge)]++;
            }
            if (i == 0) {
                starting_at[g.source(edge)]++;
            }
        }
    }
    // Put flow in the corresponding vertices covered by a path of length one (does not appear in the answer)
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        if (mu_v[v] == 0 && in_U[v]) {
            mu_v[v]++;
            starting_at[v]++;
            ending_at[v]++;
        }
    }

    ListDigraph red;

    ListDigraph::NodeMap<ListDigraph::Arc> split_edges(g);
    ListDigraph::NodeMap<ListDigraph::Node> original(red);


    // Capacities of the MaxFlow reduction
    ListDigraph::ArcMap<int64_t> feasible_flow(red, 0); // From the approximation
    ListDigraph::ArcMap<int64_t> capacities(red, 0); // For the Max-Flow reduction, it is flow-demand


    ListDigraph::Node s = red.addNode();
    ListDigraph::Node t = red.addNode();


    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        ListDigraph::Node v_in = red.addNode();
        ListDigraph::Node v_out = red.addNode();
        original[v_in] = v;
        original[v_out] = v;

        ListDigraph::Arc split = red.addArc(v_in, v_out);
        split_edges[v] = split;
        feasible_flow[split] = mu_v[v];
        capacities[split] = mu_v[v];
        if (in_U[v]) {
            capacities[split]--;
        }
    }
    for (ListDigraph::Node v : S) {
        ListDigraph::Arc sv = red.addArc(s, red.source(split_edges[v]));
        feasible_flow[sv] = starting_at[v];
        capacities[sv] = starting_at[v];
    }
    for (ListDigraph::Node v : T) {
        ListDigraph::Arc vt = red.addArc(red.target(split_edges[v]), t);
        feasible_flow[vt] = ending_at[v];
        capacities[vt] = ending_at[v];
    }

    for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        ListDigraph::Arc red_e = red.addArc(red.target(split_edges[u]), red.source(split_edges[v]));
        feasible_flow[red_e] = mu[e];
        capacities[red_e] = mu[e];
    }

    ListDigraph::ArcMap<ListDigraph::Arc> direct(red);
    ListDigraph::ArcMap<ListDigraph::Arc> reverse(red);
    std::vector<ListDigraph::Arc> edges;
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        edges.push_back(e);
    }

    for (auto e : edges) {
        ListDigraph::Arc rev_e = red.addArc(red.target(e), red.source(e));
        capacities[rev_e] = countNodes(g);
        direct[e] = e;
        direct[rev_e] = e;
        reverse[e] = rev_e;
        reverse[rev_e] = rev_e;
    }


    // Run Max-Flow algorithm
    EdmondsKarp<ListDigraph, ListDigraph::ArcMap<int64_t>> ek(red, capacities, s, t);

    // Set the flowMap to store the result in run
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ek.flowMap(flowMap);
    ek.run();


    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;
    std::vector<std::vector<ListDigraph::Arc>> path_cover_edges_red;

    // Remove 0 flow edges and st

    // Stores the modified edges with the corresponding target
    // (we use the strategy to move the target to the source instead of removing)
    std::vector<std::pair<ListDigraph::Arc , ListDigraph::Node>> restorage_list;
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        ListDigraph::Arc d_e = direct[e];
        if (e != d_e) { // If it is a reverse edge
            // The flow on that edge is computed as the previous flow minus the one discounted by that edge,
            // plus the flow in the reverse direction (discounted in the Max-flow, therefore pushed in the Min-flow)
            flowMap[d_e] = feasible_flow[d_e] + flowMap[e] - flowMap[d_e];
            flowMap[e] = 0;

            restorage_list.push_back({e, red.target(e)});
            if (flowMap[d_e] == 0) {
                restorage_list.push_back({d_e, red.target(d_e)});
            }

            capacities[e] = countNodes(g);
            capacities[d_e] = flowMap[d_e];
        }
    }
    for (ListDigraph::Node v : U) {
        ListDigraph::Arc split = split_edges[v];
        capacities[split]--;
    }
    for (auto& pair : restorage_list) {
        red.changeTarget(pair.first, red.source(pair.first));
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        std::vector<ListDigraph::Arc> path_edges_red;

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t) {
                if (path.empty() || path.back() != original[v]) {
                    path.push_back(original[v]);
                } else if (!path.empty() && path.back() == original[v]) {
                    path_edges_red.push_back(e);
                }
            } else {
                path_edges_red.push_back(e);
            }
            flowMap[e]--;
            if (flowMap[e] == 0) {
                restorage_list.push_back({e, red.target(e)});
                red.changeTarget(e, red.source(e));
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());
        std::reverse(path_edges_red.begin(), path_edges_red.end());

        path_cover.push_back(path);
        path_cover_edges_red.push_back(path_edges_red);

        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }

    int64_t width = path_cover.size();

    // Compute the paths that goes through every edge
    // These are stored as a vector of pairs {i, j}
    // where i is the position of the path in path_cover_edges,
    // and j is the position of the edge in that path
    ListDigraph::ArcMap<std::vector<std::pair<int64_t, int64_t>>> paths_through(red, {});
    for (int i = 0; i < path_cover_edges_red.size(); ++i) {
        std::vector<ListDigraph::Arc>& path_edges_red = path_cover_edges_red[i];
        for (int j = 0; j < path_edges_red.size(); ++j) {
            ListDigraph::Arc e = path_edges_red[j];
            paths_through[e].push_back({i, j});
        }
    }

    // Compute Safe Paths

    // Restore removed edges in the process of computing the MPC
    for (auto& pair: restorage_list) {
        red.changeTarget(pair.first, pair.second);
    }

    std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_per_path;


    for (int i = 0; i < path_cover.size(); ++i) {
        std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>> path_maximal_safe_paths_pair;
        std::vector<std::vector<ListDigraph::Node>> path_maximal_safe_paths;

        std::vector<ListDigraph::Node>& path = path_cover[i];
        std::vector<ListDigraph::Arc>& path_edges_red = path_cover_edges_red[i];

        int x = 0, y = 0;
        bool fail_to_expand = false;
        while (y+1 < path.size()) {
            std::vector<ListDigraph::Arc> transitive_edges;
            ListDigraph::Node x_p_1 = path[y];
            ListDigraph::Node x_p = path[y+1];
            ListDigraph::Arc e = path_edges_red[y+1];
            int64_t mu_e = capacities[e];
            ListDigraph::Arc rev_e = reverse[e];


            // Compute reduction

            red.changeTarget(e, red.source(e));
            red.changeTarget(rev_e, red.source(rev_e));
            for (int z = x+1; z <= y; ++z) {
                ListDigraph::Node v = path[z];
                for (ListDigraph::InArcIt to_v(g, v); to_v != INVALID; ++to_v) {
                    ListDigraph::Node u = g.source(to_v);
                    if (u != path[z-1]) {
                        ListDigraph::Arc tran_e = red.addArc(red.target(split_edges[u]), red.source(split_edges[x_p]));
                        direct[tran_e] = tran_e;
                        capacities[tran_e] = 0;

                        ListDigraph::Arc rev_tran_e = red.addArc(red.target(tran_e), red.source(tran_e));
                        direct[rev_tran_e] = tran_e;
                        capacities[rev_tran_e] = countNodes(g);

                        transitive_edges.push_back(tran_e);
                        transitive_edges.push_back(rev_tran_e);
                    }
                }
                if (in_S[v]) {
                    ListDigraph::Arc tran_e = red.addArc(s, red.source(split_edges[x_p]));
                    direct[tran_e] = tran_e;
                    capacities[tran_e] = 0;

                    ListDigraph::Arc rev_tran_e = red.addArc(red.target(tran_e), red.source(tran_e));
                    direct[rev_tran_e] = tran_e;
                    capacities[rev_tran_e] = countNodes(g);

                    transitive_edges.push_back(tran_e);
                    transitive_edges.push_back(rev_tran_e);
                }
            }

            // Redistribute the flow (if possible)

            // Remove reverse edges
            std::vector<std::pair<ListDigraph::Arc , ListDigraph::Node>> reverse_edges;
            for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
                ListDigraph::Arc d_e = direct[e];
                if (e != d_e) { // If it is a reverse edge
                    reverse_edges.push_back({e, red.target(e)});
                    red.changeTarget(e, red.source(e));
                }
            }


            Bfs<ListDigraph> bfs_from_s(red);
            bfs_from_s.run(s);

            // Reverse direct edges
            std::vector<ListDigraph::Arc> direct_edges;
            for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
                ListDigraph::Arc d_e = direct[e];
                if (e == d_e) { // If it is a direct edge
                    direct_edges.push_back(e);
                }
            }
            for (ListDigraph::Arc e : direct_edges) {
                ListDigraph::Node source = red.source(e);
                ListDigraph::Node target = red.target(e);
                red.changeSource(e, target);
                red.changeTarget(e, source);
            }

            Bfs<ListDigraph> bfs_to_t(red);
            bfs_to_t.run(t);

            // Reverse reversed direct edges
            for (ListDigraph::Arc e : direct_edges) {
                ListDigraph::Node source = red.source(e);
                ListDigraph::Node target = red.target(e);
                red.changeSource(e, target);
                red.changeTarget(e, source);
            }
            // Restore reverse edges
            for (auto& pair : reverse_edges) {
                red.changeTarget(pair.first, pair.second);
            }

            // For every path through e find the corresponding redistribution of flow
            bool infinite_width = bfs_from_s.predArc(t) == INVALID;

            ListDigraph::ArcMap<int64_t> flow_modification(red, 0);

            if (!infinite_width) {
                for (auto pair : paths_through[e]) {
                    int64_t i = pair.first;
                    int64_t j = pair.second;
                    auto& path = path_cover[i];
                    auto& path_edges_red = path_cover_edges_red[i];

                    // First check whether s reaches the first vertex in U in path[j...path.size()-1]
                    int64_t index_first_reached_by_s = path.size();
                    for (int64_t k = j; k < path.size(); ++k) {
                        ListDigraph::Node current_vertex = path[k];
                        if (bfs_from_s.predArc(red.source(split_edges[current_vertex])) != INVALID) {
                            index_first_reached_by_s = k;
                            break;
                        } else if (in_U[current_vertex]) {
                            infinite_width = true;
                            break;
                        }
                    }

                    if (infinite_width) break;

                    int64_t index_last_reaching_t = -1;
                    for (int64_t k = j-1; k >= 0; --k) {
                        ListDigraph::Node current_vertex = path[k];
                        if (bfs_to_t.predArc(red.target(split_edges[current_vertex])) != INVALID) {
                            index_last_reaching_t = k;
                            break;
                        } else if (in_U[current_vertex]) {
                            infinite_width = true;
                            break;
                        }
                    }

                    if (infinite_width) break;

                    // Remove 1 unit of flow from path[j...index_first_reached_by_s-1]
                    for (int64_t t = j; t < index_first_reached_by_s; ++t) {
                        capacities[split_edges[path[t]]]--;
                        capacities[path_edges_red[t+1]]--;

                        flow_modification[split_edges[path[t]]]--;
                        flow_modification[path_edges_red[t+1]]--;
                    }

                    // Push 1 unit of flow in the path from s to index_first_reached_by_s
                    ListDigraph::Node first_reached_by_s = t;
                    if (index_first_reached_by_s != path.size()) {
                        first_reached_by_s = red.source(split_edges[path[index_first_reached_by_s]]);
                    }
                    ListDigraph::Arc temp_e(INVALID);
                    ListDigraph::Node temp_v = first_reached_by_s;
                    while ((temp_e = bfs_from_s.predArc(temp_v)) != INVALID) {
                        capacities[temp_e]++;
                        flow_modification[temp_e]++;
                        temp_v = red.source(temp_e);
                    }


                    // Remove 1 unit of flow from path[index_last_reaching_t+1...j-1]
                    for (int64_t t = index_last_reaching_t; t < j-1; ++t) {
                        capacities[path_edges_red[t+1]]--;
                        capacities[split_edges[path[t+1]]]--;

                        flow_modification[path_edges_red[t+1]]--;
                        flow_modification[split_edges[path[t+1]]]--;
                    }

                    // Push 1 unit of flow in the path from index_last_reaching_t to t
                    ListDigraph::Node last_reaching_t = s;
                    if (index_last_reaching_t != -1) {
                        last_reaching_t = red.target(split_edges[path[index_last_reaching_t]]);
                    }
                    temp_v = last_reaching_t;
                    while ((temp_e = bfs_to_t.predArc(temp_v)) != INVALID) {
                        capacities[temp_e]++;
                        flow_modification[temp_e]++;
                        temp_v = red.target(temp_e);
                    }

                }
            }


            if (infinite_width) {
                //Path is safe
                // Move y to the right
                ++y;
                fail_to_expand = false;
            } else {
                // Compute new width
                EdmondsKarp<ListDigraph, ListDigraph::ArcMap<int64_t>> ek(red, capacities, s, t);

                // Set the flowMap to store the result in run
                ListDigraph::ArcMap<int64_t> flowMap(red);
                ek.flowMap(flowMap);
                ek.run();

                int64_t new_width = width + mu_e - ek.flowValue();

                if (width == new_width) { // It is not safe
                    // Report the path between x and y, move x to the right, and (if necessary) y to the right
                    if (x != y && !fail_to_expand) {
                        std::vector<ListDigraph::Node> maximal_safe_path;
                        for (int z = x; z <= y; ++z) {
                            maximal_safe_path.push_back(path[z]);
                        }
                        path_maximal_safe_paths.push_back(maximal_safe_path);
                        fail_to_expand = true;
                    }
                    ++x;
                    if (x > y) {
                        ++y;
                    }
                } else { //Path is safe
                    // Move y to the right
                    ++y;
                    fail_to_expand = false;
                }
            }

            // Here put the flow back
            for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
                capacities[e] -= flow_modification[e];
            }

            // Remove transitive edges and add e
            red.changeTarget(e, red.source(rev_e));
            red.changeTarget(rev_e, red.source(e));
            for (ListDigraph::Arc e : transitive_edges) {
                red.erase(e);
            }
        }

        // (possibly) report the last path
        if (x != y) {
            std::vector<ListDigraph::Node> maximal_safe_path;
            for (int z = x; z <= y; ++z) {
                maximal_safe_path.push_back(path[z]);
            }
            path_maximal_safe_paths.push_back(maximal_safe_path);
        }

        path_maximal_safe_paths_pair.first = path;
        path_maximal_safe_paths_pair.second = path_maximal_safe_paths;
        path_maximal_safe_paths_per_path.push_back(path_maximal_safe_paths_pair);
    }

    return path_maximal_safe_paths_per_path;
}



std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_PC(ListDigraph& g, int64_t l) {

    // Build the Min-Flow network reduction
    ListDigraph red;

    ListDigraph::NodeMap<ListDigraph::Node> v_in(g);
    ListDigraph::NodeMap<ListDigraph::Node> v_out(g);
    ListDigraph::NodeMap<ListDigraph::Node> original(red);
    ListDigraph::ArcMap<int> cost(red, 0);
    ListDigraph::ArcMap<int> demand(red, 0);
    ListDigraph::NodeMap<int> supply(red, 0);

    ListDigraph::Node s = red.addNode();
    ListDigraph::Node t = red.addNode();
    ListDigraph::Arc st = red.addArc(s, t);
    supply[s] = countNodes(g);
    supply[t] = -countNodes(g);


    // Set split vertices
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        v_in[v] = red.addNode();
        v_out[v] = red.addNode();
        original[v_in[v]] = v;
        original[v_out[v]] = v;

        ListDigraph::Arc split = red.addArc(v_in[v], v_out[v]);
        demand[split] = 1;

        ListDigraph::Arc sv = red.addArc(s, v_in[v]);
        cost[sv] = 1;

        red.addArc(v_out[v], t);
    }

    // Set edges connecting split vetices
    for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        red.addArc(v_out[u], v_in[v]);
    }



    // Use NetworkSimplex for solving the min-flow
    NetworkSimplex<ListDigraph> ns(red);
    ns.lowerMap(demand).costMap(cost).supplyMap(supply).run();
    int64_t width = ns.totalCost();

    if (width > l) { // Case where there are not safe edges at all
        return {};
    }

    // Obtain Flow solution
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ns.flowMap(flowMap);

    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;
    std::vector<std::vector<ListDigraph::Arc>> path_cover_edges_red;

    // Remove 0 flow edges and st

    // Stores the modified edges with the corresponding target
    // (we use the strategy to move the target to the source instead of removing)
    std::vector<std::pair<ListDigraph::Arc , ListDigraph::Node>> restorage_list = {{st, t}};
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        if (flowMap[e] == 0) {
            restorage_list.push_back({e, red.target(e)});
        }
    }
    for (auto& pair : restorage_list) {
        red.changeTarget(pair.first, red.source(pair.first));
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        std::vector<ListDigraph::Arc> path_edges_red;

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t) {
                if (path.empty() || path.back() != original[v]) {
                    path.push_back(original[v]);
                } else if (red.source(e) != s && (!path.empty() && path.back() == original[v])) {
                    path_edges_red.push_back(e);
                }
            }
            flowMap[e]--;
            if (flowMap[e] == 0) {
                restorage_list.push_back({e, red.target(e)});
                red.changeTarget(e, red.source(e));
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());
        std::reverse(path_edges_red.begin(), path_edges_red.end());

        path_cover.push_back(path);
        path_cover_edges_red.push_back(path_edges_red);


        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }



    // Compute Safe Paths

    // Restore removed edges
    for (auto& pair: restorage_list) {
        red.changeTarget(pair.first, pair.second);
    }

    std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_per_path;


    for (int i = 0; i < path_cover.size(); ++i) {
        std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>> path_maximal_safe_paths_pair;
        std::vector<std::vector<ListDigraph::Node>> path_maximal_safe_paths;

        std::vector<ListDigraph::Node>& path = path_cover[i];
        std::vector<ListDigraph::Arc>& path_edges_red = path_cover_edges_red[i];

        int x = 0, y = 0;
        bool fail_to_expand = false;
        while (y+1 < path.size()) {
            std::vector<ListDigraph::Arc> transitive_edges;
            ListDigraph::Node x_p = path[y+1];
            ListDigraph::Arc e = path_edges_red[y];


            // Compute reduction
            red.changeTarget(e, red.source(e));
            for (int z = x+1; z <= y; ++z) {
                ListDigraph::Node v = path[z];
                for (ListDigraph::InArcIt to_v(g, v); to_v != INVALID; ++to_v) {
                    ListDigraph::Node u = g.source(to_v);
                    if (u != path[z-1]) {
                        transitive_edges.push_back(red.addArc(v_out[u], v_in[x_p]));
                    }
                }
            }

            // Compute new width
            ns.reset();
            ns.lowerMap(demand).costMap(cost).supplyMap(supply).run();
            int64_t new_width = ns.totalCost(); // Maybe we have to check something else here for the RNA path cover case

            if (new_width  <= l) { // It is not safe
                // Report the path between x and y, move x to the right, and (if necessary) y to the right
                if (x != y && !fail_to_expand) {
                    std::vector<ListDigraph::Node> maximal_safe_path;
                    for (int z = x; z <= y; ++z) {
                        maximal_safe_path.push_back(path[z]);
                    }
                    path_maximal_safe_paths.push_back(maximal_safe_path);
                    fail_to_expand = true;
                }
                ++x;
                if (x > y) {
                    ++y;
                }
            } else { //Path is safe
                // Move y to the right
                ++y;
                fail_to_expand = false;
            }

            // Remove transitive edges and add e
            red.changeTarget(e, v_in[x_p]);
            for (ListDigraph::Arc e : transitive_edges) {
                red.erase(e);
            }
        }

        // (possibly) report the last path
        if (x != y) {
            std::vector<ListDigraph::Node> maximal_safe_path;
            for (int z = x; z <= y; ++z) {
                maximal_safe_path.push_back(path[z]);
            }
            path_maximal_safe_paths.push_back(maximal_safe_path);
        }

        path_maximal_safe_paths_pair.first = path;
        path_maximal_safe_paths_pair.second = path_maximal_safe_paths;
        path_maximal_safe_paths_per_path.push_back(path_maximal_safe_paths_pair);
    }

    return path_maximal_safe_paths_per_path;
}



std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> greedy_path_maximal_safe_paths_PC(ListDigraph& g, int64_t l) {

    // Build the Min-Flow network reduction
    std::vector<std::vector<ListDigraph::Arc>> paths = greedy_approximation_MPC_edges(g);

    // Compute mu values according to the current path cover
    ListDigraph::ArcMap<int64_t> mu(g, 0); // Number of paths using this edge
    ListDigraph::NodeMap<int64_t> mu_v(g, 0); // Number of paths using this vertex
    ListDigraph::NodeMap<int64_t> starting_at(g, 0); // Number of paths starting at this vertex
    ListDigraph::NodeMap<int64_t> ending_at(g, 0); // Number of paths ending at this vertex
    for (auto path: paths) {
        for (int i = 0; i < path.size(); ++i) {
            auto edge = path[i];
            mu[edge]++;
            mu_v[g.source(edge)]++;
            if (i == path.size()-1) {
                mu_v[g.target(edge)]++;
                ending_at[g.target(edge)]++;
            }
            if (i == 0) {
                starting_at[g.source(edge)]++;
            }
        }
    }
    // Put flow in the corresponding vertices covered by a path of length one (does not appear in the answer)
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        if (mu_v[v] == 0) {
            mu_v[v]++;
            starting_at[v]++;
            ending_at[v]++;
        }
    }

    ListDigraph red;

    ListDigraph::NodeMap<ListDigraph::Arc> split_edges(g);
    ListDigraph::NodeMap<ListDigraph::Arc> from_s(g);
    ListDigraph::NodeMap<ListDigraph::Arc> to_t(g);
    ListDigraph::NodeMap<ListDigraph::Node> original(red);


    // Capacities of the MaxFlow reduction
    ListDigraph::ArcMap<int64_t> feasible_flow(red, 0); // From the approximation
    ListDigraph::ArcMap<int64_t> capacities(red, 0); // For the Max-Flow reduction, it is flow-demand


    ListDigraph::Node s = red.addNode();
    ListDigraph::Node t = red.addNode();


    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        ListDigraph::Node v_in = red.addNode();
        ListDigraph::Node v_out = red.addNode();
        original[v_in] = v;
        original[v_out] = v;

        ListDigraph::Arc split = red.addArc(v_in, v_out);
        split_edges[v] = split;
        feasible_flow[split] = mu_v[v];
        capacities[split] = mu_v[v] - 1;

        ListDigraph::Arc sv = red.addArc(s, v_in);
        ListDigraph::Arc vt = red.addArc(v_out, t);

        from_s[v] = sv;
        to_t[v] = vt;

        feasible_flow[sv] = starting_at[v];
        feasible_flow[vt] = ending_at[v];

        capacities[sv] = starting_at[v];
        capacities[vt] = ending_at[v];
    }

    for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        ListDigraph::Arc red_e = red.addArc(red.target(split_edges[u]), red.source(split_edges[v]));
        feasible_flow[red_e] = mu[e];
        capacities[red_e] = mu[e];
    }

    ListDigraph::ArcMap<ListDigraph::Arc> direct(red);
    ListDigraph::ArcMap<ListDigraph::Arc> reverse(red);
    std::vector<ListDigraph::Arc> edges;
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        edges.push_back(e);
    }

    for (auto e : edges) {
        ListDigraph::Arc rev_e = red.addArc(red.target(e), red.source(e));
        capacities[rev_e] = countNodes(g);
        direct[e] = e;
        direct[rev_e] = e;
        reverse[e] = rev_e;
        reverse[rev_e] = rev_e;
    }


    // Run Max-Flow algorithm
    EdmondsKarp<ListDigraph, ListDigraph::ArcMap<int64_t>> ek(red, capacities, s, t);

    // Set the flowMap to store the result in run
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ek.flowMap(flowMap);
    ek.run();


    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;
    std::vector<std::vector<ListDigraph::Arc>> path_cover_edges_red;

    // Remove 0 flow edges and st

    // Stores the modified edges with the corresponding target
    // (we use the strategy to move the target to the source instead of removing)
    std::vector<std::pair<ListDigraph::Arc , ListDigraph::Node>> restorage_list;
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        ListDigraph::Arc d_e = direct[e];
        if (e != d_e) { // If it is a reverse edge
            // The flow on that edge is computed as the previous flow minus the one discounted by that edge,
            // plus the flow in the reverse direction (discounted in the Max-flow, therefore pushed in the Min-flow)
            flowMap[d_e] = feasible_flow[d_e] + flowMap[e] - flowMap[d_e];
            flowMap[e] = 0;

            restorage_list.push_back({e, red.target(e)});
            if (flowMap[d_e] == 0) {
                restorage_list.push_back({d_e, red.target(d_e)});
            }

            capacities[e] = countNodes(g);
            capacities[d_e] = flowMap[d_e];
        }
    }
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        ListDigraph::Arc split = split_edges[v];
        capacities[split]--;
    }
    for (auto& pair : restorage_list) {
        red.changeTarget(pair.first, red.source(pair.first));
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        std::vector<ListDigraph::Arc> path_edges_red;

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t) {
                if (path.empty() || path.back() != original[v]) {
                    path.push_back(original[v]);
                } else if (red.source(e) != s && (!path.empty() && path.back() == original[v])) {
                    path_edges_red.push_back(e);
                }
            }
            flowMap[e]--;
            if (flowMap[e] == 0) {
                restorage_list.push_back({e, red.target(e)});
                red.changeTarget(e, red.source(e));
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());
        std::reverse(path_edges_red.begin(), path_edges_red.end());

        path_cover.push_back(path);
        path_cover_edges_red.push_back(path_edges_red);

        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }

    int64_t width = path_cover.size();
    if (width > l) { // Case where there are not safe edges at all
        return {};
    }

    // Compute Safe Paths

    // Restore removed edges in the process of computing the MPC
    for (auto& pair: restorage_list) {
        red.changeTarget(pair.first, pair.second);
    }

    std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_per_path;


    for (int i = 0; i < path_cover.size(); ++i) {
        std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>> path_maximal_safe_paths_pair;
        std::vector<std::vector<ListDigraph::Node>> path_maximal_safe_paths;

        std::vector<ListDigraph::Node>& path = path_cover[i];
        std::vector<ListDigraph::Arc>& path_edges_red = path_cover_edges_red[i];

        int x = 0, y = 0;
        bool fail_to_expand = false;
        while (y+1 < path.size()) {
            std::vector<ListDigraph::Arc> transitive_edges;
            ListDigraph::Node x_p_1 = path[y];
            ListDigraph::Node x_p = path[y+1];
            ListDigraph::Arc e = path_edges_red[y];
            int64_t mu_e = capacities[e];
            ListDigraph::Arc rev_e = reverse[e];


            // Compute reduction
            red.changeTarget(e, red.source(e));
            red.changeTarget(rev_e, red.source(rev_e));
            capacities[from_s[x_p]] += mu_e;
            capacities[to_t[x_p_1]] += mu_e;
            for (int z = x+1; z <= y; ++z) {
                ListDigraph::Node v = path[z];
                for (ListDigraph::InArcIt to_v(g, v); to_v != INVALID; ++to_v) {
                    ListDigraph::Node u = g.source(to_v);
                    if (u != path[z-1]) {
                        ListDigraph::Arc tran_e = red.addArc(red.target(split_edges[u]), red.source(split_edges[x_p]));
                        capacities[tran_e] = 0;
                        ListDigraph::Arc rev_tran_e = red.addArc(red.target(tran_e), red.source(tran_e));
                        capacities[rev_tran_e] = countNodes(g);

                        transitive_edges.push_back(tran_e);
                        transitive_edges.push_back(rev_tran_e);
                    }
                }
            }

            // Compute new width
            EdmondsKarp<ListDigraph, ListDigraph::ArcMap<int64_t>> ek(red, capacities, s, t);

            // Set the flowMap to store the result in run
            ListDigraph::ArcMap<int64_t> flowMap(red);
            ek.flowMap(flowMap);
            ek.run();

            int64_t new_width = width + mu_e - ek.flowValue();

            if (new_width <= l) { // It is not safe
                // Report the path between x and y, move x to the right, and (if necessary) y to the right
                if (x != y && !fail_to_expand) {
                    std::vector<ListDigraph::Node> maximal_safe_path;
                    for (int z = x; z <= y; ++z) {
                        maximal_safe_path.push_back(path[z]);
                    }
                    path_maximal_safe_paths.push_back(maximal_safe_path);
                    fail_to_expand = true;
                }
                ++x;
                if (x > y) {
                    ++y;
                }
            } else { //Path is safe
                // Move y to the right
                ++y;
                fail_to_expand = false;
            }

            // Remove transitive edges and add e
            capacities[from_s[x_p]] -= mu_e;
            capacities[to_t[x_p_1]] -= mu_e;
            red.changeTarget(e, red.source(rev_e));
            red.changeTarget(rev_e, red.source(e));
            for (ListDigraph::Arc e : transitive_edges) {
                red.erase(e);
            }
        }

        // (possibly) report the last path
        if (x != y) {
            std::vector<ListDigraph::Node> maximal_safe_path;
            for (int z = x; z <= y; ++z) {
                maximal_safe_path.push_back(path[z]);
            }
            path_maximal_safe_paths.push_back(maximal_safe_path);
        }

        path_maximal_safe_paths_pair.first = path;
        path_maximal_safe_paths_pair.second = path_maximal_safe_paths;
        path_maximal_safe_paths_per_path.push_back(path_maximal_safe_paths_pair);
    }

    return path_maximal_safe_paths_per_path;
}




std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_PC(ListDigraph& g, std::vector<ListDigraph::Node>& S, std::vector<ListDigraph::Node>& T, int64_t l) {

    ListDigraph::NodeMap<bool> in_S(g, false);
    for (ListDigraph::Node v : S) {
        in_S[v] = true;
    }
    // Build the Min-Flow network reduction
    ListDigraph red;

    ListDigraph::NodeMap<ListDigraph::Node> v_in(g);
    ListDigraph::NodeMap<ListDigraph::Node> v_out(g);
    ListDigraph::NodeMap<ListDigraph::Node> original(red);
    ListDigraph::ArcMap<int> cost(red, 0);
    ListDigraph::ArcMap<int> demand(red, 0);
    ListDigraph::NodeMap<int> supply(red, 0);

    ListDigraph::Node s = red.addNode();
    ListDigraph::Node t = red.addNode();
    ListDigraph::Arc st = red.addArc(s, t);
    supply[s] = countNodes(g);
    supply[t] = -countNodes(g);


    // Set split vertices
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        v_in[v] = red.addNode();
        v_out[v] = red.addNode();
        original[v_in[v]] = v;
        original[v_out[v]] = v;

        ListDigraph::Arc split = red.addArc(v_in[v], v_out[v]);
        demand[split] = 1;
    }
    for (ListDigraph::Node v : S) {
        ListDigraph::Arc sv = red.addArc(s, v_in[v]);
        cost[sv] = 1;
    }
    for (ListDigraph::Node v : T) {
        red.addArc(v_out[v], t);
    }

    // Set edges connecting split vertices
    for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        red.addArc(v_out[u], v_in[v]);
    }



    // Use NetworkSimplex for solving the min-flow
    NetworkSimplex<ListDigraph> ns(red);
    ns.lowerMap(demand).costMap(cost).supplyMap(supply).run();
    int64_t width = ns.totalCost();
    if (width > l) { // Case where there are not safe edges at all
        return {};
    }

    // Obtain Flow solution
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ns.flowMap(flowMap);

    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;
    std::vector<std::vector<ListDigraph::Arc>> path_cover_edges_red;

    // Remove 0 flow edges and st

    // Stores the modified edges with the corresponding target
    // (we use the strategy to move the target to the source instead of removing)
    std::vector<std::pair<ListDigraph::Arc , ListDigraph::Node>> restorage_list = {{st, t}};
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        if (flowMap[e] == 0) {
            restorage_list.push_back({e, red.target(e)});
        }
    }
    for (auto& pair : restorage_list) {
        red.changeTarget(pair.first, red.source(pair.first));
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        std::vector<ListDigraph::Arc> path_edges_red;

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t) {
                if (path.empty() || path.back() != original[v]) {
                    path.push_back(original[v]);
                } else if (red.source(e) != s && (!path.empty() && path.back() == original[v])) {
                    path_edges_red.push_back(e);
                }
            }
            flowMap[e]--;
            if (flowMap[e] == 0) {
                restorage_list.push_back({e, red.target(e)});
                red.changeTarget(e, red.source(e));
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());
        std::reverse(path_edges_red.begin(), path_edges_red.end());

        path_cover.push_back(path);
        path_cover_edges_red.push_back(path_edges_red);


        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }



    // Compute Safe Paths

    // Restore removed edges
    for (auto& pair: restorage_list) {
        red.changeTarget(pair.first, pair.second);
    }

    std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_per_path;


    for (int i = 0; i < path_cover.size(); ++i) {
        std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>> path_maximal_safe_paths_pair;
        std::vector<std::vector<ListDigraph::Node>> path_maximal_safe_paths;

        std::vector<ListDigraph::Node>& path = path_cover[i];
        std::vector<ListDigraph::Arc>& path_edges_red = path_cover_edges_red[i];

        int x = 0, y = 0;
        bool fail_to_expand = false;
        while (y+1 < path.size()) {
            std::vector<ListDigraph::Arc> transitive_edges;
            ListDigraph::Node x_p = path[y+1];
            ListDigraph::Arc e = path_edges_red[y];


            // Compute reduction
            red.changeTarget(e, red.source(e));
            for (int z = x+1; z <= y; ++z) {
                ListDigraph::Node v = path[z];
                for (ListDigraph::InArcIt to_v(g, v); to_v != INVALID; ++to_v) {
                    ListDigraph::Node u = g.source(to_v);
                    if (u != path[z-1]) {
                        transitive_edges.push_back(red.addArc(v_out[u], v_in[x_p]));
                    }
                } if (in_S[v]) {
                    transitive_edges.push_back(red.addArc(s, v_in[x_p]));
                }
            }

            // Compute new width
            ns.reset();
            NetworkSimplex<ListDigraph>::ProblemType result = ns.lowerMap(demand).costMap(cost).supplyMap(supply).run();
            int64_t new_width = ns.totalCost(); // Maybe we have to check something else here for the RNA path cover case

            if (new_width <= l && result == NetworkSimplex<ListDigraph>::OPTIMAL) { // It is not safe
                // Report the path between x and y, move x to the right, and (if necessary) y to the right
                if (x != y && !fail_to_expand) {
                    std::vector<ListDigraph::Node> maximal_safe_path;
                    for (int z = x; z <= y; ++z) {
                        maximal_safe_path.push_back(path[z]);
                    }
                    path_maximal_safe_paths.push_back(maximal_safe_path);
                    fail_to_expand = true;
                }
                ++x;
                if (x > y) {
                    ++y;
                }
            } else { //Path is safe
                // Move y to the right
                ++y;
                fail_to_expand = false;
            }

            // Remove transitive edges and add e
            red.changeTarget(e, v_in[x_p]);
            for (ListDigraph::Arc e : transitive_edges) {
                red.erase(e);
            }
        }

        // (possibly) report the last path
        if (x != y) {
            std::vector<ListDigraph::Node> maximal_safe_path;
            for (int z = x; z <= y; ++z) {
                maximal_safe_path.push_back(path[z]);
            }
            path_maximal_safe_paths.push_back(maximal_safe_path);
        }

        path_maximal_safe_paths_pair.first = path;
        path_maximal_safe_paths_pair.second = path_maximal_safe_paths;
        path_maximal_safe_paths_per_path.push_back(path_maximal_safe_paths_pair);
    }

    return path_maximal_safe_paths_per_path;
}



std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> greedy_path_maximal_safe_paths_PC(ListDigraph& g, std::vector<ListDigraph::Node>& S, std::vector<ListDigraph::Node>& T, int64_t l) {

    ListDigraph::NodeMap<bool> in_S(g, false);
    for (ListDigraph::Node v : S) {
        in_S[v] = true;
    }

    // Build the Min-Flow network reduction
    std::vector<std::vector<ListDigraph::Arc>> paths = greedy_approximation_MPC_edges(g, S, T);

    // Compute mu values according to the current path cover
    ListDigraph::ArcMap<int64_t> mu(g, 0); // Number of paths using this edge
    ListDigraph::NodeMap<int64_t> mu_v(g, 0); // Number of paths using this vertex
    ListDigraph::NodeMap<int64_t> starting_at(g, 0); // Number of paths starting at this vertex
    ListDigraph::NodeMap<int64_t> ending_at(g, 0); // Number of paths ending at this vertex
    for (auto path: paths) {
        for (int i = 0; i < path.size(); ++i) {
            auto edge = path[i];
            mu[edge]++;
            mu_v[g.source(edge)]++;
            if (i == path.size()-1) {
                mu_v[g.target(edge)]++;
                ending_at[g.target(edge)]++;
            }
            if (i == 0) {
                starting_at[g.source(edge)]++;
            }
        }
    }
    // Put flow in the corresponding vertices covered by a path of length one (does not appear in the answer)
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        if (mu_v[v] == 0) {
            mu_v[v]++;
            starting_at[v]++;
            ending_at[v]++;
        }
    }

    ListDigraph red;

    ListDigraph::NodeMap<ListDigraph::Arc> split_edges(g);
    ListDigraph::NodeMap<ListDigraph::Node> original(red);


    // Capacities of the MaxFlow reduction
    ListDigraph::ArcMap<int64_t> feasible_flow(red, 0); // From the approximation
    ListDigraph::ArcMap<int64_t> capacities(red, 0); // For the Max-Flow reduction, it is flow-demand


    ListDigraph::Node s = red.addNode();
    ListDigraph::Node t = red.addNode();


    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        ListDigraph::Node v_in = red.addNode();
        ListDigraph::Node v_out = red.addNode();
        original[v_in] = v;
        original[v_out] = v;

        ListDigraph::Arc split = red.addArc(v_in, v_out);
        split_edges[v] = split;
        feasible_flow[split] = mu_v[v];
        capacities[split] = mu_v[v] - 1;
    }
    for (ListDigraph::Node v : S) {
        ListDigraph::Arc sv = red.addArc(s, red.source(split_edges[v]));
        feasible_flow[sv] = starting_at[v];
        capacities[sv] = starting_at[v];
    }
    for (ListDigraph::Node v : T) {
        ListDigraph::Arc vt = red.addArc(red.target(split_edges[v]), t);
        feasible_flow[vt] = ending_at[v];
        capacities[vt] = ending_at[v];
    }

    for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        ListDigraph::Arc red_e = red.addArc(red.target(split_edges[u]), red.source(split_edges[v]));
        feasible_flow[red_e] = mu[e];
        capacities[red_e] = mu[e];
    }

    ListDigraph::ArcMap<ListDigraph::Arc> direct(red);
    ListDigraph::ArcMap<ListDigraph::Arc> reverse(red);
    std::vector<ListDigraph::Arc> edges;
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        edges.push_back(e);
    }

    for (auto e : edges) {
        ListDigraph::Arc rev_e = red.addArc(red.target(e), red.source(e));
        capacities[rev_e] = countNodes(g);
        direct[e] = e;
        direct[rev_e] = e;
        reverse[e] = rev_e;
        reverse[rev_e] = rev_e;
    }


    // Run Max-Flow algorithm
    EdmondsKarp<ListDigraph, ListDigraph::ArcMap<int64_t>> ek(red, capacities, s, t);

    // Set the flowMap to store the result in run
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ek.flowMap(flowMap);
    ek.run();


    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;
    std::vector<std::vector<ListDigraph::Arc>> path_cover_edges_red;

    // Remove 0 flow edges and st

    // Stores the modified edges with the corresponding target
    // (we use the strategy to move the target to the source instead of removing)
    std::vector<std::pair<ListDigraph::Arc , ListDigraph::Node>> restorage_list;
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        ListDigraph::Arc d_e = direct[e];
        if (e != d_e) { // If it is a reverse edge
            // The flow on that edge is computed as the previous flow minus the one discounted by that edge,
            // plus the flow in the reverse direction (discounted in the Max-flow, therefore pushed in the Min-flow)
            flowMap[d_e] = feasible_flow[d_e] + flowMap[e] - flowMap[d_e];
            flowMap[e] = 0;

            restorage_list.push_back({e, red.target(e)});
            if (flowMap[d_e] == 0) {
                restorage_list.push_back({d_e, red.target(d_e)});
            }

            capacities[e] = countNodes(g);
            capacities[d_e] = flowMap[d_e];
        }
    }
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        ListDigraph::Arc split = split_edges[v];
        capacities[split]--;
    }
    for (auto& pair : restorage_list) {
        red.changeTarget(pair.first, red.source(pair.first));
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        std::vector<ListDigraph::Arc> path_edges_red;

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t) {
                if (path.empty() || path.back() != original[v]) {
                    path.push_back(original[v]);
                } else if (red.source(e) != s && (!path.empty() && path.back() == original[v])) {
                    path_edges_red.push_back(e);
                }
            }
            flowMap[e]--;
            if (flowMap[e] == 0) {
                restorage_list.push_back({e, red.target(e)});
                red.changeTarget(e, red.source(e));
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());
        std::reverse(path_edges_red.begin(), path_edges_red.end());

        path_cover.push_back(path);
        path_cover_edges_red.push_back(path_edges_red);

        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }

    int64_t width = path_cover.size();
    if (width > l) { // Case where there are not safe edges at all
        return {};
    }

    // Compute Safe Paths

    // Restore removed edges in the process of computing the MPC
    for (auto& pair: restorage_list) {
        red.changeTarget(pair.first, pair.second);
    }

    std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_per_path;


    for (int i = 0; i < path_cover.size(); ++i) {
        std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>> path_maximal_safe_paths_pair;
        std::vector<std::vector<ListDigraph::Node>> path_maximal_safe_paths;

        std::vector<ListDigraph::Node>& path = path_cover[i];
        std::vector<ListDigraph::Arc>& path_edges_red = path_cover_edges_red[i];

        int x = 0, y = 0;
        bool fail_to_expand = false;
        while (y+1 < path.size()) {
            std::vector<ListDigraph::Arc> transitive_edges;
            ListDigraph::Node x_p_1 = path[y];
            ListDigraph::Node x_p = path[y+1];
            ListDigraph::Arc e = path_edges_red[y];
            int64_t mu_e = capacities[e];
            ListDigraph::Arc rev_e = reverse[e];


            // Compute reduction

            red.changeTarget(e, red.source(e));
            red.changeTarget(rev_e, red.source(rev_e));
            for (int z = x+1; z <= y; ++z) {
                ListDigraph::Node v = path[z];
                for (ListDigraph::InArcIt to_v(g, v); to_v != INVALID; ++to_v) {
                    ListDigraph::Node u = g.source(to_v);
                    if (u != path[z-1]) {
                        ListDigraph::Arc tran_e = red.addArc(red.target(split_edges[u]), red.source(split_edges[x_p]));
                        direct[tran_e] = tran_e;
                        capacities[tran_e] = 0;

                        ListDigraph::Arc rev_tran_e = red.addArc(red.target(tran_e), red.source(tran_e));
                        direct[rev_tran_e] = tran_e;
                        capacities[rev_tran_e] = countNodes(g);

                        transitive_edges.push_back(tran_e);
                        transitive_edges.push_back(rev_tran_e);
                    }
                }
                if (in_S[v]) {
                    ListDigraph::Arc tran_e = red.addArc(s, red.source(split_edges[x_p]));
                    direct[tran_e] = tran_e;
                    capacities[tran_e] = 0;

                    ListDigraph::Arc rev_tran_e = red.addArc(red.target(tran_e), red.source(tran_e));
                    direct[rev_tran_e] = tran_e;
                    capacities[rev_tran_e] = countNodes(g);

                    transitive_edges.push_back(tran_e);
                    transitive_edges.push_back(rev_tran_e);
                }
            }

            // Find path from s to x_p in red and push mu_e flow through it
            std::vector<ListDigraph::Arc> s_x_p_path = path_through_direct_edges(red, direct, s, red.source(rev_e));
            for (ListDigraph::Arc temp_e: s_x_p_path) {
                capacities[temp_e] += mu_e;
            }
            // Find path from x_p_1 to t in red and push mu_e flow through it
            std::vector<ListDigraph::Arc> x_p_1_t_path = path_through_direct_edges(red, direct, red.source(e), t);
            for (ListDigraph::Arc temp_e: x_p_1_t_path) {
                capacities[temp_e] += mu_e;
            }


            if (s_x_p_path.size() == 0 || x_p_1_t_path.size() == 0) {
                //Path is safe
                // Move y to the right
                ++y;
                fail_to_expand = false;
            } else {
                // Compute new width
                EdmondsKarp<ListDigraph, ListDigraph::ArcMap<int64_t>> ek(red, capacities, s, t);

                // Set the flowMap to store the result in run
                ListDigraph::ArcMap<int64_t> flowMap(red);
                ek.flowMap(flowMap);
                ek.run();

                int64_t new_width = width + mu_e - ek.flowValue();

                if (new_width <= l ) { // It is not safe
                    // Report the path between x and y, move x to the right, and (if necessary) y to the right
                    if (x != y && !fail_to_expand) {
                        std::vector<ListDigraph::Node> maximal_safe_path;
                        for (int z = x; z <= y; ++z) {
                            maximal_safe_path.push_back(path[z]);
                        }
                        path_maximal_safe_paths.push_back(maximal_safe_path);
                        fail_to_expand = true;
                    }
                    ++x;
                    if (x > y) {
                        ++y;
                    }
                } else { //Path is safe
                    // Move y to the right
                    ++y;
                    fail_to_expand = false;
                }
            }


            // Remove transitive edges and add e
            for (ListDigraph::Arc temp_e : s_x_p_path) {
                capacities[temp_e] -= mu_e;
            }
            for (ListDigraph::Arc temp_e : x_p_1_t_path) {
                capacities[temp_e] -= mu_e;
            }
            red.changeTarget(e, red.source(rev_e));
            red.changeTarget(rev_e, red.source(e));
            for (ListDigraph::Arc e : transitive_edges) {
                red.erase(e);
            }
        }

        // (possibly) report the last path
        if (x != y) {
            std::vector<ListDigraph::Node> maximal_safe_path;
            for (int z = x; z <= y; ++z) {
                maximal_safe_path.push_back(path[z]);
            }
            path_maximal_safe_paths.push_back(maximal_safe_path);
        }

        path_maximal_safe_paths_pair.first = path;
        path_maximal_safe_paths_pair.second = path_maximal_safe_paths;
        path_maximal_safe_paths_per_path.push_back(path_maximal_safe_paths_pair);
    }

    return path_maximal_safe_paths_per_path;
}



std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_U_PC(ListDigraph& g, std::vector<lemon::ListDigraph::Node>& U, int64_t l) {

    // Compute in_U
    ListDigraph::NodeMap<bool> in_U(g, false);
    for (ListDigraph::Node v : U) {
        in_U[v] = true;
    }

    // Build the Min-Flow network reduction
    ListDigraph red;

    ListDigraph::NodeMap<ListDigraph::Node> v_in(g);
    ListDigraph::NodeMap<ListDigraph::Node> v_out(g);
    ListDigraph::NodeMap<ListDigraph::Node> original(red);
    ListDigraph::ArcMap<int> cost(red, 0);
    ListDigraph::ArcMap<int> demand(red, 0);
    ListDigraph::NodeMap<int> supply(red, 0);

    ListDigraph::Node s = red.addNode();
    ListDigraph::Node t = red.addNode();
    ListDigraph::Arc st = red.addArc(s, t);
    supply[s] = countNodes(g);
    supply[t] = -countNodes(g);


    // Set split vertices
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        v_in[v] = red.addNode();
        v_out[v] = red.addNode();
        original[v_in[v]] = v;
        original[v_out[v]] = v;

        ListDigraph::Arc split = red.addArc(v_in[v], v_out[v]);
        if (in_U[v]) {
            demand[split] = 1;
        }

        ListDigraph::Arc sv = red.addArc(s, v_in[v]);
        cost[sv] = 1;

        red.addArc(v_out[v], t);
    }

    // Set edges connecting split vertices
    for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        red.addArc(v_out[u], v_in[v]);
    }



    // Use NetworkSimplex for solving the min-flow
    NetworkSimplex<ListDigraph> ns(red);
    ns.lowerMap(demand).costMap(cost).supplyMap(supply).run();
    int64_t width = ns.totalCost();
    if (width > l) { // Case where there are not safe edges at all
        return {};
    }

    // Obtain Flow solution
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ns.flowMap(flowMap);

    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;
    std::vector<std::vector<ListDigraph::Arc>> path_cover_edges_red;

    // Remove 0 flow edges and st

    // Stores the modified edges with the corresponding target
    // (we use the strategy to move the target to the source instead of removing)
    std::vector<std::pair<ListDigraph::Arc , ListDigraph::Node>> restorage_list = {{st, t}};
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        if (flowMap[e] == 0) {
            restorage_list.push_back({e, red.target(e)});
        }
    }
    for (auto& pair : restorage_list) {
        red.changeTarget(pair.first, red.source(pair.first));
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        std::vector<ListDigraph::Arc> path_edges_red;

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t) {
                if (path.empty() || path.back() != original[v]) {
                    path.push_back(original[v]);
                } else if (red.source(e) != s && (!path.empty() && path.back() == original[v])) {
                    path_edges_red.push_back(e);
                }
            }
            flowMap[e]--;
            if (flowMap[e] == 0) {
                restorage_list.push_back({e, red.target(e)});
                red.changeTarget(e, red.source(e));
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());
        std::reverse(path_edges_red.begin(), path_edges_red.end());

        path_cover.push_back(path);
        path_cover_edges_red.push_back(path_edges_red);


        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }



    // Compute Safe Paths

    // Restore removed edges in the process of computing the MPC
    for (auto& pair: restorage_list) {
        red.changeTarget(pair.first, pair.second);
    }

    std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_per_path;


    for (int i = 0; i < path_cover.size(); ++i) {
        std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>> path_maximal_safe_paths_pair;
        std::vector<std::vector<ListDigraph::Node>> path_maximal_safe_paths;

        std::vector<ListDigraph::Node>& path = path_cover[i];
        std::vector<ListDigraph::Arc>& path_edges_red = path_cover_edges_red[i];

        int x = 0, y = 0;
        bool fail_to_expand = false;
        while (y+1 < path.size()) {
            std::vector<ListDigraph::Arc> transitive_edges;
            ListDigraph::Node x_p = path[y+1];
            ListDigraph::Arc e = path_edges_red[y];


            // Compute reduction
            red.changeTarget(e, red.source(e));
            for (int z = x+1; z <= y; ++z) {
                ListDigraph::Node v = path[z];
                for (ListDigraph::InArcIt to_v(g, v); to_v != INVALID; ++to_v) {
                    ListDigraph::Node u = g.source(to_v);
                    if (u != path[z-1]) {
                        transitive_edges.push_back(red.addArc(v_out[u], v_in[x_p]));
                    }
                }
            }

            // Compute new width
            ns.reset();
            ns.lowerMap(demand).costMap(cost).supplyMap(supply).run();
            int64_t new_width = ns.totalCost(); // Maybe we have to check something else here for the RNA path cover case

            if (new_width <= l) { // It is not safe
                // Report the path between x and y, move x to the right, and (if necessary) y to the right
                if (x != y && !fail_to_expand) {
                    std::vector<ListDigraph::Node> maximal_safe_path;
                    for (int z = x; z <= y; ++z) {
                        maximal_safe_path.push_back(path[z]);
                    }
                    path_maximal_safe_paths.push_back(maximal_safe_path);
                    fail_to_expand = true;
                }
                ++x;
                if (x > y) {
                    ++y;
                }
            } else { //Path is safe
                // Move y to the right
                ++y;
                fail_to_expand = false;
            }

            // Remove transitive edges and add e
            red.changeTarget(e, v_in[x_p]);
            for (ListDigraph::Arc e : transitive_edges) {
                red.erase(e);
            }
        }

        // (possibly) report the last path
        if (x != y) {
            std::vector<ListDigraph::Node> maximal_safe_path;
            for (int z = x; z <= y; ++z) {
                maximal_safe_path.push_back(path[z]);
            }
            path_maximal_safe_paths.push_back(maximal_safe_path);
        }

        path_maximal_safe_paths_pair.first = path;
        path_maximal_safe_paths_pair.second = path_maximal_safe_paths;
        path_maximal_safe_paths_per_path.push_back(path_maximal_safe_paths_pair);
    }

    return path_maximal_safe_paths_per_path;
}



std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> greedy_path_maximal_safe_paths_U_PC(ListDigraph& g, std::vector<lemon::ListDigraph::Node>& U, int64_t l) {

    // Compute in_U
    ListDigraph::NodeMap<bool> in_U(g, false);
    for (ListDigraph::Node v : U) {
        in_U[v] = true;
    }

    // Build the Min-Flow network reduction
    std::vector<std::vector<ListDigraph::Arc>> paths = greedy_approximation_U_MPC_edges(g, U);

    // Compute mu values according to the current path cover
    ListDigraph::ArcMap<int64_t> mu(g, 0); // Number of paths using this edge
    ListDigraph::NodeMap<int64_t> mu_v(g, 0); // Number of paths using this vertex
    ListDigraph::NodeMap<int64_t> starting_at(g, 0); // Number of paths starting at this vertex
    ListDigraph::NodeMap<int64_t> ending_at(g, 0); // Number of paths ending at this vertex
    for (auto path: paths) {
        for (int i = 0; i < path.size(); ++i) {
            auto edge = path[i];
            mu[edge]++;
            mu_v[g.source(edge)]++;
            if (i == path.size()-1) {
                mu_v[g.target(edge)]++;
                ending_at[g.target(edge)]++;
            }
            if (i == 0) {
                starting_at[g.source(edge)]++;
            }
        }
    }
    // Put flow in the corresponding vertices covered by a path of length one (does not appear in the answer)
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        if (mu_v[v] == 0 && in_U[v]) {
            mu_v[v]++;
            starting_at[v]++;
            ending_at[v]++;
        }
    }

    ListDigraph red;

    ListDigraph::NodeMap<ListDigraph::Arc> split_edges(g);
    ListDigraph::NodeMap<ListDigraph::Arc> from_s(g);
    ListDigraph::NodeMap<ListDigraph::Arc> to_t(g);
    ListDigraph::NodeMap<ListDigraph::Node> original(red);


    // Capacities of the MaxFlow reduction
    ListDigraph::ArcMap<int64_t> feasible_flow(red, 0); // From the approximation
    ListDigraph::ArcMap<int64_t> capacities(red, 0); // For the Max-Flow reduction, it is flow-demand


    ListDigraph::Node s = red.addNode();
    ListDigraph::Node t = red.addNode();


    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        ListDigraph::Node v_in = red.addNode();
        ListDigraph::Node v_out = red.addNode();
        original[v_in] = v;
        original[v_out] = v;

        ListDigraph::Arc split = red.addArc(v_in, v_out);
        split_edges[v] = split;
        feasible_flow[split] = mu_v[v];
        capacities[split] = mu_v[v];
        if (in_U[v]) {
            capacities[split]--;
        }

        ListDigraph::Arc sv = red.addArc(s, v_in);
        ListDigraph::Arc vt = red.addArc(v_out, t);

        from_s[v] = sv;
        to_t[v] = vt;

        feasible_flow[sv] = starting_at[v];
        feasible_flow[vt] = ending_at[v];

        capacities[sv] = starting_at[v];
        capacities[vt] = ending_at[v];
    }

    for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        ListDigraph::Arc red_e = red.addArc(red.target(split_edges[u]), red.source(split_edges[v]));
        feasible_flow[red_e] = mu[e];
        capacities[red_e] = mu[e];
    }

    ListDigraph::ArcMap<ListDigraph::Arc> direct(red);
    ListDigraph::ArcMap<ListDigraph::Arc> reverse(red);
    std::vector<ListDigraph::Arc> edges;
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        edges.push_back(e);
    }

    for (auto e : edges) {
        ListDigraph::Arc rev_e = red.addArc(red.target(e), red.source(e));
        capacities[rev_e] = countNodes(g);
        direct[e] = e;
        direct[rev_e] = e;
        reverse[e] = rev_e;
        reverse[rev_e] = rev_e;
    }


    // Run Max-Flow algorithm
    EdmondsKarp<ListDigraph, ListDigraph::ArcMap<int64_t>> ek(red, capacities, s, t);

    // Set the flowMap to store the result in run
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ek.flowMap(flowMap);
    ek.run();


    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;
    std::vector<std::vector<ListDigraph::Arc>> path_cover_edges_red;

    // Remove 0 flow edges and st

    // Stores the modified edges with the corresponding target
    // (we use the strategy to move the target to the source instead of removing)
    std::vector<std::pair<ListDigraph::Arc , ListDigraph::Node>> restorage_list;
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        ListDigraph::Arc d_e = direct[e];
        if (e != d_e) { // If it is a reverse edge
            // The flow on that edge is computed as the previous flow minus the one discounted by that edge,
            // plus the flow in the reverse direction (discounted in the Max-flow, therefore pushed in the Min-flow)
            flowMap[d_e] = feasible_flow[d_e] + flowMap[e] - flowMap[d_e];
            flowMap[e] = 0;

            restorage_list.push_back({e, red.target(e)});
            if (flowMap[d_e] == 0) {
                restorage_list.push_back({d_e, red.target(d_e)});
            }

            capacities[e] = countNodes(g);
            capacities[d_e] = flowMap[d_e];
        }
    }
    for (ListDigraph::Node v : U) {
        ListDigraph::Arc split = split_edges[v];
        capacities[split]--;
    }
    for (auto& pair : restorage_list) {
        red.changeTarget(pair.first, red.source(pair.first));
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        std::vector<ListDigraph::Arc> path_edges_red;

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t) {
                if (path.empty() || path.back() != original[v]) {
                    path.push_back(original[v]);
                } else if (red.source(e) != s && (!path.empty() && path.back() == original[v])) {
                    path_edges_red.push_back(e);
                }
            }
            flowMap[e]--;
            if (flowMap[e] == 0) {
                restorage_list.push_back({e, red.target(e)});
                red.changeTarget(e, red.source(e));
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());
        std::reverse(path_edges_red.begin(), path_edges_red.end());

        path_cover.push_back(path);
        path_cover_edges_red.push_back(path_edges_red);

        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }

    int64_t width = path_cover.size();
    if (width > l) { // Case where there are not safe edges at all
        return {};
    }

    // Compute Safe Paths

    // Restore removed edges in the process of computing the MPC
    for (auto& pair: restorage_list) {
        red.changeTarget(pair.first, pair.second);
    }

    std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_per_path;


    for (int i = 0; i < path_cover.size(); ++i) {
        std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>> path_maximal_safe_paths_pair;
        std::vector<std::vector<ListDigraph::Node>> path_maximal_safe_paths;

        std::vector<ListDigraph::Node>& path = path_cover[i];
        std::vector<ListDigraph::Arc>& path_edges_red = path_cover_edges_red[i];

        int x = 0, y = 0;
        bool fail_to_expand = false;
        while (y+1 < path.size()) {
            std::vector<ListDigraph::Arc> transitive_edges;
            ListDigraph::Node x_p_1 = path[y];
            ListDigraph::Node x_p = path[y+1];
            ListDigraph::Arc e = path_edges_red[y];
            int64_t mu_e = capacities[e];
            ListDigraph::Arc rev_e = reverse[e];


            // Compute reduction
            red.changeTarget(e, red.source(e));
            red.changeTarget(rev_e, red.source(rev_e));
            capacities[from_s[x_p]] += mu_e;
            capacities[to_t[x_p_1]] += mu_e;
            for (int z = x+1; z <= y; ++z) {
                ListDigraph::Node v = path[z];
                for (ListDigraph::InArcIt to_v(g, v); to_v != INVALID; ++to_v) {
                    ListDigraph::Node u = g.source(to_v);
                    if (u != path[z-1]) {
                        ListDigraph::Arc tran_e = red.addArc(red.target(split_edges[u]), red.source(split_edges[x_p]));
                        capacities[tran_e] = 0;
                        ListDigraph::Arc rev_tran_e = red.addArc(red.target(tran_e), red.source(tran_e));
                        capacities[rev_tran_e] = countNodes(g);

                        transitive_edges.push_back(tran_e);
                        transitive_edges.push_back(rev_tran_e);
                    }
                }
            }

            // Compute new width
            EdmondsKarp<ListDigraph, ListDigraph::ArcMap<int64_t>> ek(red, capacities, s, t);

            // Set the flowMap to store the result in run
            ListDigraph::ArcMap<int64_t> flowMap(red);
            ek.flowMap(flowMap);
            ek.run();

            int64_t new_width = width + mu_e - ek.flowValue();

            if (new_width <= l) { // It is not safe
                // Report the path between x and y, move x to the right, and (if necessary) y to the right
                if (x != y && !fail_to_expand) {
                    std::vector<ListDigraph::Node> maximal_safe_path;
                    for (int z = x; z <= y; ++z) {
                        maximal_safe_path.push_back(path[z]);
                    }
                    path_maximal_safe_paths.push_back(maximal_safe_path);
                    fail_to_expand = true;
                }
                ++x;
                if (x > y) {
                    ++y;
                }
            } else { //Path is safe
                // Move y to the right
                ++y;
                fail_to_expand = false;
            }

            // Remove transitive edges and add e
            capacities[from_s[x_p]] -= mu_e;
            capacities[to_t[x_p_1]] -= mu_e;
            red.changeTarget(e, red.source(rev_e));
            red.changeTarget(rev_e, red.source(e));
            for (ListDigraph::Arc e : transitive_edges) {
                red.erase(e);
            }
        }

        // (possibly) report the last path
        if (x != y) {
            std::vector<ListDigraph::Node> maximal_safe_path;
            for (int z = x; z <= y; ++z) {
                maximal_safe_path.push_back(path[z]);
            }
            path_maximal_safe_paths.push_back(maximal_safe_path);
        }

        path_maximal_safe_paths_pair.first = path;
        path_maximal_safe_paths_pair.second = path_maximal_safe_paths;
        path_maximal_safe_paths_per_path.push_back(path_maximal_safe_paths_pair);
    }

    return path_maximal_safe_paths_per_path;
}



std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_U_PC(ListDigraph& g, std::vector<ListDigraph::Node>& S, std::vector<ListDigraph::Node>& T, std::vector<lemon::ListDigraph::Node>& U, int64_t l) {

    // Compute in_U
    ListDigraph::NodeMap<bool> in_U(g, false);
    ListDigraph::NodeMap<bool> in_S(g, false);
    for (ListDigraph::Node v : U) {
        in_U[v] = true;
    }
    for (ListDigraph::Node v : S) {
        in_S[v] = true;
    }

    // Build the Min-Flow network reduction
    ListDigraph red;

    ListDigraph::NodeMap<ListDigraph::Node> v_in(g);
    ListDigraph::NodeMap<ListDigraph::Node> v_out(g);
    ListDigraph::NodeMap<ListDigraph::Node> original(red);
    ListDigraph::ArcMap<int> cost(red, 0);
    ListDigraph::ArcMap<int> demand(red, 0);
    ListDigraph::NodeMap<int> supply(red, 0);

    ListDigraph::Node s = red.addNode();
    ListDigraph::Node t = red.addNode();
    ListDigraph::Arc st = red.addArc(s, t);
    supply[s] = countNodes(g);
    supply[t] = -countNodes(g);


    // Set split vertices
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        v_in[v] = red.addNode();
        v_out[v] = red.addNode();
        original[v_in[v]] = v;
        original[v_out[v]] = v;

        ListDigraph::Arc split = red.addArc(v_in[v], v_out[v]);
        if (in_U[v]) {
            demand[split] = 1;
        }
    }
    for (ListDigraph::Node v : S) {
        ListDigraph::Arc sv = red.addArc(s, v_in[v]);
        cost[sv] = 1;
    }
    for (ListDigraph::Node v : T) {
        red.addArc(v_out[v], t);
    }

    // Set edges connecting split vertices
    for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        red.addArc(v_out[u], v_in[v]);
    }



    // Use NetworkSimplex for solving the min-flow
    NetworkSimplex<ListDigraph> ns(red);
    ns.lowerMap(demand).costMap(cost).supplyMap(supply).run();
    int64_t width = ns.totalCost();
    if (width > l) { // Case where there are not safe edges at all
        return {};
    }

    // Obtain Flow solution
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ns.flowMap(flowMap);

    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;
    std::vector<std::vector<ListDigraph::Arc>> path_cover_edges_red;

    // Remove 0 flow edges and st

    // Stores the modified edges with the corresponding target
    // (we use the strategy to move the target to the source instead of removing)
    std::vector<std::pair<ListDigraph::Arc , ListDigraph::Node>> restorage_list = {{st, t}};
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        if (flowMap[e] == 0) {
            restorage_list.push_back({e, red.target(e)});
        }
    }
    for (auto& pair : restorage_list) {
        red.changeTarget(pair.first, red.source(pair.first));
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        std::vector<ListDigraph::Arc> path_edges_red;

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t) {
                if (path.empty() || path.back() != original[v]) {
                    path.push_back(original[v]);
                } else if (red.source(e) != s && (!path.empty() && path.back() == original[v])) {
                    path_edges_red.push_back(e);
                }
            }
            flowMap[e]--;
            if (flowMap[e] == 0) {
                restorage_list.push_back({e, red.target(e)});
                red.changeTarget(e, red.source(e));
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());
        std::reverse(path_edges_red.begin(), path_edges_red.end());

        path_cover.push_back(path);
        path_cover_edges_red.push_back(path_edges_red);


        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }



    // Compute Safe Paths

    // Restore removed edges
    for (auto& pair: restorage_list) {
        red.changeTarget(pair.first, pair.second);
    }

    std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_per_path;


    for (int i = 0; i < path_cover.size(); ++i) {
        std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>> path_maximal_safe_paths_pair;
        std::vector<std::vector<ListDigraph::Node>> path_maximal_safe_paths;

        std::vector<ListDigraph::Node>& path = path_cover[i];
        std::vector<ListDigraph::Arc>& path_edges_red = path_cover_edges_red[i];

        int x = 0, y = 0;
        bool fail_to_expand = false;
        while (y+1 < path.size()) {
            std::vector<ListDigraph::Arc> transitive_edges;
            ListDigraph::Node x_p = path[y+1];
            ListDigraph::Arc e = path_edges_red[y];


            // Compute reduction
            red.changeTarget(e, red.source(e));
            for (int z = x+1; z <= y; ++z) {
                ListDigraph::Node v = path[z];
                for (ListDigraph::InArcIt to_v(g, v); to_v != INVALID; ++to_v) {
                    ListDigraph::Node u = g.source(to_v);
                    if (u != path[z-1]) {
                        transitive_edges.push_back(red.addArc(v_out[u], v_in[x_p]));
                    }
                } if (in_S[v]) {
                    transitive_edges.push_back(red.addArc(s, v_in[x_p]));
                }
            }

            // Compute new width
            ns.reset();
            NetworkSimplex<ListDigraph>::ProblemType result = ns.lowerMap(demand).costMap(cost).supplyMap(supply).run();
            int64_t new_width = ns.totalCost(); // Maybe we have to check something else here for the RNA path cover case

            if (new_width <= l && result == NetworkSimplex<ListDigraph>::OPTIMAL) { // It is not safe
                // Report the path between x and y, move x to the right, and (if necessary) y to the right
                if (x != y && !fail_to_expand) {
                    std::vector<ListDigraph::Node> maximal_safe_path;
                    for (int z = x; z <= y; ++z) {
                        maximal_safe_path.push_back(path[z]);
                    }
                    path_maximal_safe_paths.push_back(maximal_safe_path);
                    fail_to_expand = true;
                }
                ++x;
                if (x > y) {
                    ++y;
                }
            } else { //Path is safe
                // Move y to the right
                ++y;
                fail_to_expand = false;
            }

            // Remove transitive edges and add e
            red.changeTarget(e, v_in[x_p]);
            for (ListDigraph::Arc e : transitive_edges) {
                red.erase(e);
            }
        }

        // (possibly) report the last path
        if (x != y) {
            std::vector<ListDigraph::Node> maximal_safe_path;
            for (int z = x; z <= y; ++z) {
                maximal_safe_path.push_back(path[z]);
            }
            path_maximal_safe_paths.push_back(maximal_safe_path);
        }

        path_maximal_safe_paths_pair.first = path;
        path_maximal_safe_paths_pair.second = path_maximal_safe_paths;
        path_maximal_safe_paths_per_path.push_back(path_maximal_safe_paths_pair);
    }

    return path_maximal_safe_paths_per_path;
}



std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> greedy_path_maximal_safe_paths_U_PC(ListDigraph& g, std::vector<ListDigraph::Node>& S, std::vector<ListDigraph::Node>& T, std::vector<lemon::ListDigraph::Node>& U, int64_t l) {

    // Compute in_U
    ListDigraph::NodeMap<bool> in_U(g, false);
    ListDigraph::NodeMap<bool> in_S(g, false);
    for (ListDigraph::Node v : U) {
        in_U[v] = true;
    }
    for (ListDigraph::Node v : S) {
        in_S[v] = true;
    }

    // Build the Min-Flow network reduction
    std::vector<std::vector<ListDigraph::Arc>> paths = greedy_approximation_U_MPC_edges(g, S, T, U);

    // Compute mu values according to the current path cover
    ListDigraph::ArcMap<int64_t> mu(g, 0); // Number of paths using this edge
    ListDigraph::NodeMap<int64_t> mu_v(g, 0); // Number of paths using this vertex
    ListDigraph::NodeMap<int64_t> starting_at(g, 0); // Number of paths starting at this vertex
    ListDigraph::NodeMap<int64_t> ending_at(g, 0); // Number of paths ending at this vertex
    for (auto path: paths) {
        for (int i = 0; i < path.size(); ++i) {
            auto edge = path[i];
            mu[edge]++;
            mu_v[g.source(edge)]++;
            if (i == path.size()-1) {
                mu_v[g.target(edge)]++;
                ending_at[g.target(edge)]++;
            }
            if (i == 0) {
                starting_at[g.source(edge)]++;
            }
        }
    }
    // Put flow in the corresponding vertices covered by a path of length one (does not appear in the answer)
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        if (mu_v[v] == 0 && in_U[v]) {
            mu_v[v]++;
            starting_at[v]++;
            ending_at[v]++;
        }
    }

    ListDigraph red;

    ListDigraph::NodeMap<ListDigraph::Arc> split_edges(g);
    ListDigraph::NodeMap<ListDigraph::Node> original(red);


    // Capacities of the MaxFlow reduction
    ListDigraph::ArcMap<int64_t> feasible_flow(red, 0); // From the approximation
    ListDigraph::ArcMap<int64_t> capacities(red, 0); // For the Max-Flow reduction, it is flow-demand


    ListDigraph::Node s = red.addNode();
    ListDigraph::Node t = red.addNode();


    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        ListDigraph::Node v_in = red.addNode();
        ListDigraph::Node v_out = red.addNode();
        original[v_in] = v;
        original[v_out] = v;

        ListDigraph::Arc split = red.addArc(v_in, v_out);
        split_edges[v] = split;
        feasible_flow[split] = mu_v[v];
        capacities[split] = mu_v[v];
        if (in_U[v]) {
            capacities[split]--;
        }
    }
    for (ListDigraph::Node v : S) {
        ListDigraph::Arc sv = red.addArc(s, red.source(split_edges[v]));
        feasible_flow[sv] = starting_at[v];
        capacities[sv] = starting_at[v];
    }
    for (ListDigraph::Node v : T) {
        ListDigraph::Arc vt = red.addArc(red.target(split_edges[v]), t);
        feasible_flow[vt] = ending_at[v];
        capacities[vt] = ending_at[v];
    }

    for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        ListDigraph::Arc red_e = red.addArc(red.target(split_edges[u]), red.source(split_edges[v]));
        feasible_flow[red_e] = mu[e];
        capacities[red_e] = mu[e];
    }

    ListDigraph::ArcMap<ListDigraph::Arc> direct(red);
    ListDigraph::ArcMap<ListDigraph::Arc> reverse(red);
    std::vector<ListDigraph::Arc> edges;
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        edges.push_back(e);
    }

    for (auto e : edges) {
        ListDigraph::Arc rev_e = red.addArc(red.target(e), red.source(e));
        capacities[rev_e] = countNodes(g);
        direct[e] = e;
        direct[rev_e] = e;
        reverse[e] = rev_e;
        reverse[rev_e] = rev_e;
    }


    // Run Max-Flow algorithm
    EdmondsKarp<ListDigraph, ListDigraph::ArcMap<int64_t>> ek(red, capacities, s, t);

    // Set the flowMap to store the result in run
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ek.flowMap(flowMap);
    ek.run();


    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;
    std::vector<std::vector<ListDigraph::Arc>> path_cover_edges_red;

    // Remove 0 flow edges and st

    // Stores the modified edges with the corresponding target
    // (we use the strategy to move the target to the source instead of removing)
    std::vector<std::pair<ListDigraph::Arc , ListDigraph::Node>> restorage_list;
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        ListDigraph::Arc d_e = direct[e];
        if (e != d_e) { // If it is a reverse edge
            // The flow on that edge is computed as the previous flow minus the one discounted by that edge,
            // plus the flow in the reverse direction (discounted in the Max-flow, therefore pushed in the Min-flow)
            flowMap[d_e] = feasible_flow[d_e] + flowMap[e] - flowMap[d_e];
            flowMap[e] = 0;

            restorage_list.push_back({e, red.target(e)});
            if (flowMap[d_e] == 0) {
                restorage_list.push_back({d_e, red.target(d_e)});
            }

            capacities[e] = countNodes(g);
            capacities[d_e] = flowMap[d_e];
        }
    }
    for (ListDigraph::Node v : U) {
        ListDigraph::Arc split = split_edges[v];
        capacities[split]--;
    }
    for (auto& pair : restorage_list) {
        red.changeTarget(pair.first, red.source(pair.first));
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        std::vector<ListDigraph::Arc> path_edges_red;

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t) {
                if (path.empty() || path.back() != original[v]) {
                    path.push_back(original[v]);
                } else if (!path.empty() && path.back() == original[v]) {
                    path_edges_red.push_back(e);
                }
            } else {
                path_edges_red.push_back(e);
            }
            flowMap[e]--;
            if (flowMap[e] == 0) {
                restorage_list.push_back({e, red.target(e)});
                red.changeTarget(e, red.source(e));
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());
        std::reverse(path_edges_red.begin(), path_edges_red.end());

        path_cover.push_back(path);
        path_cover_edges_red.push_back(path_edges_red);

        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }

    int64_t width = path_cover.size();
    if (width > l) { // Case where there are not safe edges at all
        return {};
    }

    // Compute the paths that goes through every edge
    // These are stored as a vector of pairs {i, j}
    // where i is the position of the path in path_cover_edges,
    // and j is the position of the edge in that path
    ListDigraph::ArcMap<std::vector<std::pair<int64_t, int64_t>>> paths_through(red, {});
    for (int i = 0; i < path_cover_edges_red.size(); ++i) {
        std::vector<ListDigraph::Arc>& path_edges_red = path_cover_edges_red[i];
        for (int j = 0; j < path_edges_red.size(); ++j) {
            ListDigraph::Arc e = path_edges_red[j];
            paths_through[e].push_back({i, j});
        }
    }

    // Compute Safe Paths

    // Restore removed edges in the process of computing the MPC
    for (auto& pair: restorage_list) {
        red.changeTarget(pair.first, pair.second);
    }

    std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_per_path;

    for (int i = 0; i < path_cover.size(); ++i) {
        std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>> path_maximal_safe_paths_pair;
        std::vector<std::vector<ListDigraph::Node>> path_maximal_safe_paths;

        std::vector<ListDigraph::Node>& path = path_cover[i];
        std::vector<ListDigraph::Arc>& path_edges_red = path_cover_edges_red[i];

        int x = 0, y = 0;
        bool fail_to_expand = false;
        while (y+1 < path.size()) {
            std::vector<ListDigraph::Arc> transitive_edges;
            ListDigraph::Node x_p_1 = path[y];
            ListDigraph::Node x_p = path[y+1];
            ListDigraph::Arc e = path_edges_red[y+1];
            int64_t mu_e = capacities[e];
            ListDigraph::Arc rev_e = reverse[e];


            // Compute reduction

            red.changeTarget(e, red.source(e));
            red.changeTarget(rev_e, red.source(rev_e));
            for (int z = x+1; z <= y; ++z) {
                ListDigraph::Node v = path[z];
                for (ListDigraph::InArcIt to_v(g, v); to_v != INVALID; ++to_v) {
                    ListDigraph::Node u = g.source(to_v);
                    if (u != path[z-1]) {
                        ListDigraph::Arc tran_e = red.addArc(red.target(split_edges[u]), red.source(split_edges[x_p]));
                        direct[tran_e] = tran_e;
                        capacities[tran_e] = 0;

                        ListDigraph::Arc rev_tran_e = red.addArc(red.target(tran_e), red.source(tran_e));
                        direct[rev_tran_e] = tran_e;
                        capacities[rev_tran_e] = countNodes(g);

                        transitive_edges.push_back(tran_e);
                        transitive_edges.push_back(rev_tran_e);
                    }
                }
                if (in_S[v]) {
                    ListDigraph::Arc tran_e = red.addArc(s, red.source(split_edges[x_p]));
                    direct[tran_e] = tran_e;
                    capacities[tran_e] = 0;

                    ListDigraph::Arc rev_tran_e = red.addArc(red.target(tran_e), red.source(tran_e));
                    direct[rev_tran_e] = tran_e;
                    capacities[rev_tran_e] = countNodes(g);

                    transitive_edges.push_back(tran_e);
                    transitive_edges.push_back(rev_tran_e);
                }
            }

            // Redistribute the flow (if possible)

            // Remove reverse edges
            std::vector<std::pair<ListDigraph::Arc , ListDigraph::Node>> reverse_edges;
            for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
                ListDigraph::Arc d_e = direct[e];
                if (e != d_e) { // If it is a reverse edge
                    reverse_edges.push_back({e, red.target(e)});
                    red.changeTarget(e, red.source(e));
                }
            }


            Bfs<ListDigraph> bfs_from_s(red);
            bfs_from_s.run(s);

            // Reverse direct edges
            std::vector<ListDigraph::Arc> direct_edges;
            for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
                ListDigraph::Arc d_e = direct[e];
                if (e == d_e) { // If it is a direct edge
                    direct_edges.push_back(e);
                }
            }
            for (ListDigraph::Arc e : direct_edges) {
                ListDigraph::Node source = red.source(e);
                ListDigraph::Node target = red.target(e);
                red.changeSource(e, target);
                red.changeTarget(e, source);
            }

            Bfs<ListDigraph> bfs_to_t(red);
            bfs_to_t.run(t);

            // Reverse reversed direct edges
            for (ListDigraph::Arc e : direct_edges) {
                ListDigraph::Node source = red.source(e);
                ListDigraph::Node target = red.target(e);
                red.changeSource(e, target);
                red.changeTarget(e, source);
            }
            // Restore reverse edges
            for (auto& pair : reverse_edges) {
                red.changeTarget(pair.first, pair.second);
            }

            // For every path through e find the corresponding redistribution of flow
            bool infinite_width = bfs_from_s.predArc(t) == INVALID;

            ListDigraph::ArcMap<int64_t> flow_modification(red, 0);

            if (!infinite_width) {
                for (auto pair : paths_through[e]) {
                    int64_t i = pair.first;
                    int64_t j = pair.second;
                    auto& path = path_cover[i];
                    auto& path_edges_red = path_cover_edges_red[i];

                    // First check whether s reaches the first vertex in U in path[j...path.size()-1]
                    int64_t index_first_reached_by_s = path.size();
                    for (int64_t k = j; k < path.size(); ++k) {
                        ListDigraph::Node current_vertex = path[k];
                        if (bfs_from_s.predArc(red.source(split_edges[current_vertex])) != INVALID) {
                            index_first_reached_by_s = k;
                            break;
                        } else if (in_U[current_vertex]) {
                            infinite_width = true;
                            break;
                        }
                    }

                    if (infinite_width) break;

                    int64_t index_last_reaching_t = -1;
                    for (int64_t k = j-1; k >= 0; --k) {
                        ListDigraph::Node current_vertex = path[k];
                        if (bfs_to_t.predArc(red.target(split_edges[current_vertex])) != INVALID) {
                            index_last_reaching_t = k;
                            break;
                        } else if (in_U[current_vertex]) {
                            infinite_width = true;
                            break;
                        }
                    }

                    if (infinite_width) break;

                    // Remove 1 unit of flow from path[j...index_first_reached_by_s-1]
                    for (int64_t t = j; t < index_first_reached_by_s; ++t) {
                        capacities[split_edges[path[t]]]--;
                        capacities[path_edges_red[t+1]]--;

                        flow_modification[split_edges[path[t]]]--;
                        flow_modification[path_edges_red[t+1]]--;
                    }

                    // Push 1 unit of flow in the path from s to index_first_reached_by_s
                    ListDigraph::Node first_reached_by_s = t;
                    if (index_first_reached_by_s != path.size()) {
                        first_reached_by_s = red.source(split_edges[path[index_first_reached_by_s]]);
                    }
                    ListDigraph::Arc temp_e(INVALID);
                    ListDigraph::Node temp_v = first_reached_by_s;
                    while ((temp_e = bfs_from_s.predArc(temp_v)) != INVALID) {
                        capacities[temp_e]++;
                        flow_modification[temp_e]++;
                        temp_v = red.source(temp_e);
                    }


                    // Remove 1 unit of flow from path[index_last_reaching_t+1...j-1]
                    for (int64_t t = index_last_reaching_t; t < j-1; ++t) {
                        capacities[path_edges_red[t+1]]--;
                        capacities[split_edges[path[t+1]]]--;

                        flow_modification[path_edges_red[t+1]]--;
                        flow_modification[split_edges[path[t+1]]]--;
                    }

                    // Push 1 unit of flow in the path from index_last_reaching_t to t
                    ListDigraph::Node last_reaching_t = s;
                    if (index_last_reaching_t != -1) {
                        last_reaching_t = red.target(split_edges[path[index_last_reaching_t]]);
                    }
                    temp_v = last_reaching_t;
                    while ((temp_e = bfs_to_t.predArc(temp_v)) != INVALID) {
                        capacities[temp_e]++;
                        flow_modification[temp_e]++;
                        temp_v = red.target(temp_e);
                    }

                }
            }


            if (infinite_width) {
                //Path is safe
                // Move y to the right
                ++y;
                fail_to_expand = false;
            } else {
                // Compute new width
                EdmondsKarp<ListDigraph, ListDigraph::ArcMap<int64_t>> ek(red, capacities, s, t);

                // Set the flowMap to store the result in run
                ListDigraph::ArcMap<int64_t> flowMap(red);
                ek.flowMap(flowMap);
                ek.run();

                int64_t new_width = width + mu_e - ek.flowValue();

                if (new_width <= l) { // It is not safe
                    // Report the path between x and y, move x to the right, and (if necessary) y to the right
                    if (x != y && !fail_to_expand) {
                        std::vector<ListDigraph::Node> maximal_safe_path;
                        for (int z = x; z <= y; ++z) {
                            maximal_safe_path.push_back(path[z]);
                        }
                        path_maximal_safe_paths.push_back(maximal_safe_path);
                        fail_to_expand = true;
                    }
                    ++x;
                    if (x > y) {
                        ++y;
                    }
                } else { //Path is safe
                    // Move y to the right
                    ++y;
                    fail_to_expand = false;
                }
            }

            // Here put the flow back
            for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
                capacities[e] -= flow_modification[e];
            }

            // Remove transitive edges and add e
            red.changeTarget(e, red.source(rev_e));
            red.changeTarget(rev_e, red.source(e));
            for (ListDigraph::Arc e : transitive_edges) {
                red.erase(e);
            }
        }

        // (possibly) report the last path
        if (x != y) {
            std::vector<ListDigraph::Node> maximal_safe_path;
            for (int z = x; z <= y; ++z) {
                maximal_safe_path.push_back(path[z]);
            }
            path_maximal_safe_paths.push_back(maximal_safe_path);
        }

        path_maximal_safe_paths_pair.first = path;
        path_maximal_safe_paths_pair.second = path_maximal_safe_paths;
        path_maximal_safe_paths_per_path.push_back(path_maximal_safe_paths_pair);
    }

    return path_maximal_safe_paths_per_path;
}



std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> optimized_greedy_path_maximal_safe_paths_U_PC(ListDigraph& g, std::vector<ListDigraph::Node>& S, std::vector<ListDigraph::Node>& T, std::vector<lemon::ListDigraph::Node>& U, int64_t l) {

    // Compute in_U
    ListDigraph::NodeMap<bool> in_U(g, false);
    ListDigraph::NodeMap<bool> in_S(g, false);
    for (ListDigraph::Node v : U) {
        in_U[v] = true;
    }
    for (ListDigraph::Node v : S) {
        in_S[v] = true;
    }

    // Build the Min-Flow network reduction
    std::vector<std::vector<ListDigraph::Arc>> paths = greedy_approximation_U_MPC_edges(g, S, T, U);

    // Compute mu values according to the current path cover
    ListDigraph::ArcMap<int64_t> mu(g, 0); // Number of paths using this edge
    ListDigraph::NodeMap<int64_t> mu_v(g, 0); // Number of paths using this vertex
    ListDigraph::NodeMap<int64_t> starting_at(g, 0); // Number of paths starting at this vertex
    ListDigraph::NodeMap<int64_t> ending_at(g, 0); // Number of paths ending at this vertex
    for (auto path: paths) {
        for (int i = 0; i < path.size(); ++i) {
            auto edge = path[i];
            mu[edge]++;
            mu_v[g.source(edge)]++;
            if (i == path.size()-1) {
                mu_v[g.target(edge)]++;
                ending_at[g.target(edge)]++;
            }
            if (i == 0) {
                starting_at[g.source(edge)]++;
            }
        }
    }
    // Put flow in the corresponding vertices covered by a path of length one (does not appear in the answer)
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        if (mu_v[v] == 0 && in_U[v]) {
            mu_v[v]++;
            starting_at[v]++;
            ending_at[v]++;
        }
    }

    ListDigraph red;

    ListDigraph::NodeMap<ListDigraph::Arc> split_edges(g);
    ListDigraph::NodeMap<ListDigraph::Node> original(red);


    // Capacities of the MaxFlow reduction
    ListDigraph::ArcMap<int64_t> feasible_flow(red, 0); // From the approximation
    ListDigraph::ArcMap<int64_t> capacities(red, 0); // For the Max-Flow reduction, it is flow-demand
    ListDigraph::ArcMap<bool> safe_edge(red, false); // true if that edge is safe

    ListDigraph::Node s = red.addNode();
    ListDigraph::Node t = red.addNode();


    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        ListDigraph::Node v_in = red.addNode();
        ListDigraph::Node v_out = red.addNode();
        original[v_in] = v;
        original[v_out] = v;

        ListDigraph::Arc split = red.addArc(v_in, v_out);
        split_edges[v] = split;
        feasible_flow[split] = mu_v[v];
        capacities[split] = mu_v[v];
        if (in_U[v]) {
            capacities[split]--;
        }
    }
    for (ListDigraph::Node v : S) {
        ListDigraph::Arc sv = red.addArc(s, red.source(split_edges[v]));
        feasible_flow[sv] = starting_at[v];
        capacities[sv] = starting_at[v];
        safe_edge[sv] = true;
    }
    for (ListDigraph::Node v : T) {
        ListDigraph::Arc vt = red.addArc(red.target(split_edges[v]), t);
        feasible_flow[vt] = ending_at[v];
        capacities[vt] = ending_at[v];
        safe_edge[vt] = true;
    }

    for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        ListDigraph::Arc red_e = red.addArc(red.target(split_edges[u]), red.source(split_edges[v]));
        feasible_flow[red_e] = mu[e];
        capacities[red_e] = mu[e];
    }

    ListDigraph::ArcMap<ListDigraph::Arc> direct(red);
    ListDigraph::ArcMap<ListDigraph::Arc> reverse(red);
    std::vector<ListDigraph::Arc> edges;
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        edges.push_back(e);
    }



    for (auto e : edges) {
        ListDigraph::Arc rev_e = red.addArc(red.target(e), red.source(e));
        capacities[rev_e] = countNodes(g);
        direct[e] = e;
        direct[rev_e] = e;
        reverse[e] = rev_e;
        reverse[rev_e] = rev_e;

    }




    // Run Max-Flow algorithm
    EdmondsKarp<ListDigraph, ListDigraph::ArcMap<int64_t>> ek(red, capacities, s, t);

    // Set the flowMap to store the result in run
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ek.flowMap(flowMap);
    ek.run();


    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;
    std::vector<std::vector<ListDigraph::Arc>> path_cover_edges_red;

    // Remove 0 flow edges and st

    // Stores the modified edges with the corresponding target
    // (we use the strategy to move the target to the source instead of removing)
    std::vector<std::pair<ListDigraph::Arc , ListDigraph::Node>> restorage_list;
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        ListDigraph::Arc d_e = direct[e];
        if (e != d_e) { // If it is a reverse edge
            // The flow on that edge is computed as the previous flow minus the one discounted by that edge,
            // plus the flow in the reverse direction (discounted in the Max-flow, therefore pushed in the Min-flow)
            flowMap[d_e] = feasible_flow[d_e] + flowMap[e] - flowMap[d_e];
            flowMap[e] = 0;

            restorage_list.push_back({e, red.target(e)});
            if (flowMap[d_e] == 0) {
                restorage_list.push_back({d_e, red.target(d_e)});
            }

            capacities[e] = countNodes(g);
            capacities[d_e] = flowMap[d_e];
        }
    }
    for (ListDigraph::Node v : U) {
        ListDigraph::Arc split = split_edges[v];
        capacities[split]--;
    }
    for (auto& pair : restorage_list) {
        red.changeTarget(pair.first, red.source(pair.first));
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        std::vector<ListDigraph::Arc> path_edges_red;

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t) {
                if (path.empty() || path.back() != original[v]) {
                    path.push_back(original[v]);
                } else if (!path.empty() && path.back() == original[v]) {
                    path_edges_red.push_back(e);
                }
            } else {
                path_edges_red.push_back(e);
            }
            flowMap[e]--;
            if (flowMap[e] == 0) {
                restorage_list.push_back({e, red.target(e)});
                red.changeTarget(e, red.source(e));
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());
        std::reverse(path_edges_red.begin(), path_edges_red.end());

        path_cover.push_back(path);
        path_cover_edges_red.push_back(path_edges_red);

        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }

    int64_t width = path_cover.size();
    if (width > l) { // Case where there are not safe edges at all
        return {};
    }

    // Compute the paths that goes through every edge
    // These are stored as a vector of pairs {i, j}
    // where i is the position of the path in path_cover_edges,
    // and j is the position of the edge in that path
    ListDigraph::ArcMap<std::vector<std::pair<int64_t, int64_t>>> paths_through(red, {});
    for (int i = 0; i < path_cover_edges_red.size(); ++i) {
        std::vector<ListDigraph::Arc>& path_edges_red = path_cover_edges_red[i];
        for (int j = 0; j < path_edges_red.size(); ++j) {
            ListDigraph::Arc e = path_edges_red[j];
            paths_through[e].push_back({i, j});
        }
    }



    // Compute Safe Paths

    // Restore removed edges in the process of computing the MPC
    for (auto& pair: restorage_list) {
        red.changeTarget(pair.first, pair.second);
    }


    // Compute safe edges

    // The outgoing edges from s and ingoing t to are safe and are set before//
    for (ListDigraph::Arc e : edges) {
        ListDigraph::Node u_out = red.source(e);
        ListDigraph::Node v_in = red.target(e);

        if (paths_through[e].size() !=0 && u_out != s && v_in != t) { // paths_through[e].size() != 0 iff e is an edge of the path cover
            int64_t mu_e = capacities[e];
            ListDigraph::Arc rev_e = reverse[e];

            // Compute reduction

            red.changeTarget(e, red.source(e));
            red.changeTarget(rev_e, red.source(rev_e));

            // Redistribute the flow (if possible)

            // Remove reverse edges
            std::vector<std::pair<ListDigraph::Arc , ListDigraph::Node>> reverse_edges;
            for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
                ListDigraph::Arc d_e = direct[e];
                if (e != d_e) { // If it is a reverse edge
                    reverse_edges.push_back({e, red.target(e)});
                    red.changeTarget(e, red.source(e));
                }
            }

            Bfs<ListDigraph> bfs_from_s(red);
            bfs_from_s.run(s);

            // Reverse direct edges
            std::vector<ListDigraph::Arc> direct_edges;
            for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
                ListDigraph::Arc d_e = direct[e];
                if (e == d_e) { // If it is a direct edge
                    direct_edges.push_back(e);
                }
            }
            for (ListDigraph::Arc e : direct_edges) {
                ListDigraph::Node source = red.source(e);
                ListDigraph::Node target = red.target(e);
                red.changeSource(e, target);
                red.changeTarget(e, source);
            }

            Bfs<ListDigraph> bfs_to_t(red);
            bfs_to_t.run(t);

            // Reverse reversed direct edges
            for (ListDigraph::Arc e : direct_edges) {
                ListDigraph::Node source = red.source(e);
                ListDigraph::Node target = red.target(e);
                red.changeSource(e, target);
                red.changeTarget(e, source);
            }
            // Restore reverse edges
            for (auto& pair : reverse_edges) {
                red.changeTarget(pair.first, pair.second);
            }

            // For every path through e find the corresponding redistribution of flow
            bool infinite_width = bfs_from_s.predArc(t) == INVALID;

            ListDigraph::ArcMap<int64_t> flow_modification(red, 0);

            if (!infinite_width) {
                for (auto pair : paths_through[e]) {
                    int64_t i = pair.first;
                    int64_t j = pair.second;
                    auto& path = path_cover[i];
                    auto& path_edges_red = path_cover_edges_red[i];

                    // First check whether s reaches the first vertex in U in path[j...path.size()-1]
                    int64_t index_first_reached_by_s = path.size();
                    for (int64_t k = j; k < path.size(); ++k) {
                        ListDigraph::Node current_vertex = path[k];
                        if (bfs_from_s.predArc(red.source(split_edges[current_vertex])) != INVALID) {
                            index_first_reached_by_s = k;
                            break;
                        } else if (in_U[current_vertex]) {
                            infinite_width = true;
                            break;
                        }
                    }

                    if (infinite_width) break;

                    int64_t index_last_reaching_t = -1;
                    for (int64_t k = j-1; k >= 0; --k) {
                        ListDigraph::Node current_vertex = path[k];
                        if (bfs_to_t.predArc(red.target(split_edges[current_vertex])) != INVALID) {
                            index_last_reaching_t = k;
                            break;
                        } else if (in_U[current_vertex]) {
                            infinite_width = true;
                            break;
                        }
                    }

                    if (infinite_width) break;

                    // Remove 1 unit of flow from path[j...index_first_reached_by_s-1]
                    for (int64_t t = j; t < index_first_reached_by_s; ++t) {
                        capacities[split_edges[path[t]]]--;
                        capacities[path_edges_red[t+1]]--;

                        flow_modification[split_edges[path[t]]]--;
                        flow_modification[path_edges_red[t+1]]--;
                    }

                    // Push 1 unit of flow in the path from s to index_first_reached_by_s
                    ListDigraph::Node first_reached_by_s = t;
                    if (index_first_reached_by_s != path.size()) {
                        first_reached_by_s = red.source(split_edges[path[index_first_reached_by_s]]);
                    }
                    ListDigraph::Arc temp_e(INVALID);
                    ListDigraph::Node temp_v = first_reached_by_s;
                    while ((temp_e = bfs_from_s.predArc(temp_v)) != INVALID) {
                        capacities[temp_e]++;
                        flow_modification[temp_e]++;
                        temp_v = red.source(temp_e);
                    }


                    // Remove 1 unit of flow from path[index_last_reaching_t+1...j-1]
                    for (int64_t t = index_last_reaching_t; t < j-1; ++t) {
                        capacities[path_edges_red[t+1]]--;
                        capacities[split_edges[path[t+1]]]--;

                        flow_modification[path_edges_red[t+1]]--;
                        flow_modification[split_edges[path[t+1]]]--;
                    }

                    // Push 1 unit of flow in the path from index_last_reaching_t to t
                    ListDigraph::Node last_reaching_t = s;
                    if (index_last_reaching_t != -1) {
                        last_reaching_t = red.target(split_edges[path[index_last_reaching_t]]);
                    }
                    temp_v = last_reaching_t;
                    while ((temp_e = bfs_to_t.predArc(temp_v)) != INVALID) {
                        capacities[temp_e]++;
                        flow_modification[temp_e]++;
                        temp_v = red.target(temp_e);
                    }

                }
            }


            if (infinite_width) {
                //Edge is safe
                safe_edge[e] = true;
            }
            else {
                // Compute new width
                EdmondsKarp<ListDigraph, ListDigraph::ArcMap<int64_t>> ek(red, capacities, s, t);

                // Set the flowMap to store the result in run
                ListDigraph::ArcMap<int64_t> flowMap(red);
                ek.flowMap(flowMap);
                ek.run();

                int64_t new_width = width + mu_e - ek.flowValue();

                if (new_width > l) {
                    //Edge is safe
                    safe_edge[e] = true;
                }
            }

            // Here put the flow back
            for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
                capacities[e] -= flow_modification[e];
            }

            // Remove transitive edges and add e
            red.changeTarget(e, red.source(rev_e));
            red.changeTarget(rev_e, red.source(e));

        }
    }

    // At this point we have that safe safe_edge[e] is true for every edge e
    // in the MPC that it is safe, therefore we can run our optimization


    std::vector<std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>>> path_maximal_safe_paths_per_path;



    for (int i = 0; i < path_cover.size(); ++i) {

        std::pair<std::vector<ListDigraph::Node>, std::vector<std::vector<ListDigraph::Node>>> path_maximal_safe_paths_pair;
        std::vector<std::vector<ListDigraph::Node>> path_maximal_safe_paths;

        std::vector<ListDigraph::Node>& path = path_cover[i];
        std::vector<ListDigraph::Arc>& path_edges_red = path_cover_edges_red[i];

        int x = 0, y = 0;
        bool fail_to_expand = false;
        while (y+1 < path.size()) {

            ListDigraph::Arc e = path_edges_red[y+1];

            if (!safe_edge[e]) {
                // Report the path between x and y, move x to the right, and (if necessary) y to the right
                if (x != y && !fail_to_expand) {
                    std::vector<ListDigraph::Node> maximal_safe_path;
                    for (int z = x; z <= y; ++z) {
                        maximal_safe_path.push_back(path[z]);
                    }
                    path_maximal_safe_paths.push_back(maximal_safe_path);
                    fail_to_expand = true;
                }
                x = y+1;
                y = y+1;
            } else {
                std::vector<ListDigraph::Arc> transitive_edges;
                ListDigraph::Node x_p_1 = path[y];
                ListDigraph::Node x_p = path[y+1];
                ListDigraph::Arc e = path_edges_red[y+1];
                int64_t mu_e = capacities[e];
                ListDigraph::Arc rev_e = reverse[e];


                // Compute reduction

                red.changeTarget(e, red.source(e));
                red.changeTarget(rev_e, red.source(rev_e));
                for (int z = x+1; z <= y; ++z) {
                    ListDigraph::Node v = path[z];
                    for (ListDigraph::InArcIt to_v(g, v); to_v != INVALID; ++to_v) {
                        ListDigraph::Node u = g.source(to_v);
                        if (u != path[z-1]) {
                            ListDigraph::Arc tran_e = red.addArc(red.target(split_edges[u]), red.source(split_edges[x_p]));
                            direct[tran_e] = tran_e;
                            capacities[tran_e] = 0;

                            ListDigraph::Arc rev_tran_e = red.addArc(red.target(tran_e), red.source(tran_e));
                            direct[rev_tran_e] = tran_e;
                            reverse[rev_tran_e] = rev_tran_e;
                            reverse[tran_e] = rev_tran_e;
                            capacities[rev_tran_e] = countNodes(g);

                            transitive_edges.push_back(tran_e);
                            transitive_edges.push_back(rev_tran_e);
                        }
                    } if (in_S[v]) {
                        ListDigraph::Arc tran_e = red.addArc(s, red.source(split_edges[x_p]));
                        direct[tran_e] = tran_e;
                        capacities[tran_e] = 0;

                        ListDigraph::Arc rev_tran_e = red.addArc(red.target(tran_e), red.source(tran_e));
                        direct[rev_tran_e] = tran_e;
                        capacities[rev_tran_e] = countNodes(g);

                        transitive_edges.push_back(tran_e);
                        transitive_edges.push_back(rev_tran_e);
                    }
                }



                // Redistribute the flow (if possible)

                // Remove reverse edges
                std::vector<std::pair<ListDigraph::Arc , ListDigraph::Node>> reverse_edges;
                for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
                    ListDigraph::Arc d_e = direct[e];
                    if (e != d_e) { // If it is a reverse edge
                        reverse_edges.push_back({e, red.target(e)});
                        red.changeTarget(e, red.source(e));
                    }
                }


                Bfs<ListDigraph> bfs_from_s(red);
                bfs_from_s.run(s);


                // Reverse direct edges
                std::vector<ListDigraph::Arc> direct_edges;
                for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
                    ListDigraph::Arc d_e = direct[e];
                    if (e == d_e) { // If it is a direct edge
                        direct_edges.push_back(e);
                    }
                }
                for (ListDigraph::Arc e : direct_edges) {
                    ListDigraph::Node source = red.source(e);
                    ListDigraph::Node target = red.target(e);
                    red.changeSource(e, target);
                    red.changeTarget(e, source);
                }

                Bfs<ListDigraph> bfs_to_t(red);
                bfs_to_t.run(t);

                // Reverse reversed direct edges
                for (ListDigraph::Arc e : direct_edges) {
                    ListDigraph::Node source = red.source(e);
                    ListDigraph::Node target = red.target(e);
                    red.changeSource(e, target);
                    red.changeTarget(e, source);
                }
                // Restore reverse edges
                for (auto& pair : reverse_edges) {
                    red.changeTarget(pair.first, pair.second);
                }


                // For every path through e find the corresponding redistribution of flow
                bool infinite_width = bfs_from_s.predArc(t) == INVALID;

                ListDigraph::ArcMap<int64_t> flow_modification(red, 0);

                if (!infinite_width) {

                    for (auto pair : paths_through[e]) {
                        int64_t i = pair.first;
                        int64_t j = pair.second;
                        auto& path = path_cover[i];
                        auto& path_edges_red = path_cover_edges_red[i];

                        // First check whether s reaches the first vertex in U in path[j...path.size()-1]
                        int64_t index_first_reached_by_s = path.size();
                        for (int64_t k = j; k < path.size(); ++k) {
                            ListDigraph::Node current_vertex = path[k];
                            if (bfs_from_s.predArc(red.source(split_edges[current_vertex])) != INVALID) {
                                index_first_reached_by_s = k;
                                break;
                            } else if (in_U[current_vertex]) {
                                infinite_width = true;
                                break;
                            }
                        }

                        if (infinite_width) break;


                        int64_t index_last_reaching_t = -1;
                        for (int64_t k = j-1; k >= 0; --k) {
                            ListDigraph::Node current_vertex = path[k];
                            if (bfs_to_t.predArc(red.target(split_edges[current_vertex])) != INVALID) {
                                index_last_reaching_t = k;
                                break;
                            } else if (in_U[current_vertex]) {
                                infinite_width = true;
                                break;
                            }
                        }

                        if (infinite_width) break;



                        // Remove 1 unit of flow from path[j...index_first_reached_by_s-1]
                        for (int64_t t = j; t < index_first_reached_by_s; ++t) {
                            capacities[split_edges[path[t]]]--;
                            capacities[path_edges_red[t+1]]--;

                            flow_modification[split_edges[path[t]]]--;
                            flow_modification[path_edges_red[t+1]]--;
                        }



                        // Push 1 unit of flow in the path from s to index_first_reached_by_s
                        ListDigraph::Node first_reached_by_s = t;
                        if (index_first_reached_by_s != path.size()) {
                            first_reached_by_s = red.source(split_edges[path[index_first_reached_by_s]]);
                        }
                        ListDigraph::Arc temp_e(INVALID);
                        ListDigraph::Node temp_v = first_reached_by_s;
                        while ((temp_e = bfs_from_s.predArc(temp_v)) != INVALID) {
                            capacities[temp_e]++;
                            flow_modification[temp_e]++;
                            temp_v = red.source(temp_e);
                        }



                        // Remove 1 unit of flow from path[index_last_reaching_t+1...j-1]
                        for (int64_t t = index_last_reaching_t; t < j-1; ++t) {
                            capacities[path_edges_red[t+1]]--;
                            capacities[split_edges[path[t+1]]]--;

                            flow_modification[path_edges_red[t+1]]--;
                            flow_modification[split_edges[path[t+1]]]--;
                        }

                        // Push 1 unit of flow in the path from index_last_reaching_t to t
                        ListDigraph::Node last_reaching_t = s;
                        if (index_last_reaching_t != -1) {
                            last_reaching_t = red.target(split_edges[path[index_last_reaching_t]]);
                        }
                        temp_v = last_reaching_t;
                        while ((temp_e = bfs_to_t.predArc(temp_v)) != INVALID) {
                            capacities[temp_e]++;
                            flow_modification[temp_e]++;
                            temp_v = red.target(temp_e);
                        }

                    }

                }


                if (infinite_width) {
                    //Path is safe
                    // Move y to the right
                    ++y;
                    fail_to_expand = false;
                }
                else {

                    // Compute new width
                    EdmondsKarp<ListDigraph, ListDigraph::ArcMap<int64_t>> ek(red, capacities, s, t);

                    // Set the flowMap to store the result in run
                    ListDigraph::ArcMap<int64_t> flowMap(red);
                    ek.flowMap(flowMap);
                    ek.run();


                    int64_t new_width = width + mu_e - ek.flowValue();


                    if (new_width <= l) { // It is not safe
                        // Report the path between x and y, move x to the right, and (if necessary) y to the right
                        if (x != y && !fail_to_expand) {
                            std::vector<ListDigraph::Node> maximal_safe_path;
                            for (int z = x; z <= y; ++z) {
                                maximal_safe_path.push_back(path[z]);
                            }
                            path_maximal_safe_paths.push_back(maximal_safe_path);
                            fail_to_expand = true;
                        }
                        ++x;
                        if (x > y) {
                            ++y;
                        }
                    } else { //Path is safe
                        // Move y to the right
                        ++y;
                        fail_to_expand = false;
                    }
                }

                // Here put the flow back
                for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
                    capacities[e] -= flow_modification[e];
                }

                // Remove transitive edges and add e
                red.changeTarget(e, red.source(rev_e));
                red.changeTarget(rev_e, red.source(e));
                for (ListDigraph::Arc e : transitive_edges) {
                    red.erase(e);
                }
            }


        }

        // (possibly) report the last path
        if (x != y) {
            std::vector<ListDigraph::Node> maximal_safe_path;
            for (int z = x; z <= y; ++z) {
                maximal_safe_path.push_back(path[z]);
            }
            path_maximal_safe_paths.push_back(maximal_safe_path);
        }

        path_maximal_safe_paths_pair.first = path;
        path_maximal_safe_paths_pair.second = path_maximal_safe_paths;
        path_maximal_safe_paths_per_path.push_back(path_maximal_safe_paths_pair);
    }

    return path_maximal_safe_paths_per_path;
}