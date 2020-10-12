#include <algorithms/safe_edges.h>

#include <algorithms/greedy_approx.h>

#include <lemon/edmonds_karp.h>
#include <lemon/dfs.h>
#include <lemon/bfs.h>


using namespace lemon;



std::vector<lemon::ListDigraph::Arc> greedy_safe_edges_U_PC(ListDigraph& g, std::vector<ListDigraph::Node>& S, std::vector<ListDigraph::Node>& T, std::vector<lemon::ListDigraph::Node>& U, int64_t l) {

    // Compute in_U
    ListDigraph::NodeMap<bool> in_U(g, false);
    for (ListDigraph::Node v : U) {
        in_U[v] = true;
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
    ListDigraph::ArcMap<ListDigraph::Arc> original_edge(red);


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
        original_edge[red_e] = e;
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




    // Restore removed edges in the process of computing the MPC
    for (auto& pair: restorage_list) {
        red.changeTarget(pair.first, pair.second);
    }

    // Compute Safe Edges
    std::vector<ListDigraph::Arc> safe_edges;


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
                safe_edges.push_back(original_edge[e]);
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
                    safe_edges.push_back(original_edge[e]);
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

    return safe_edges;
}