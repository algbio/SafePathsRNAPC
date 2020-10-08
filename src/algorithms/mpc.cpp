#include <algorithms/mpc.h>

#include <algorithms/greedy_approx.h>

#include <lemon/network_simplex.h>
#include <lemon/edmonds_karp.h>
#include <lemon/dfs.h>


using namespace lemon;



std::vector<std::vector<ListDigraph::Node>> MPC(ListDigraph& g) {

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


    // Obtain Flow solution
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ns.flowMap(flowMap);

    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;

    // Remove 0 flow edges and st
    std::vector<ListDigraph::Arc> removal_list = {st};
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        if (flowMap[e] == 0) {
            removal_list.push_back(e);
        }
    }
    for (ListDigraph::Arc e : removal_list) {
        red.erase(e);
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        removal_list.clear();

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t && (path.empty() || path.back() != original[v]))
                path.push_back(original[v]);
            flowMap[e]--;
            if (flowMap[e] == 0) {
                removal_list.push_back(e);
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());

        path_cover.push_back(path);
        for (ListDigraph::Arc e : removal_list) {
            red.erase(e);
        }

        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }

    return path_cover;
}



std::vector<std::vector<ListDigraph::Node>> greedy_MPC(ListDigraph& g)  {

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

    ListDigraph::NodeMap<ListDigraph::Node> v_in(g);
    ListDigraph::NodeMap<ListDigraph::Node> v_out(g);
    ListDigraph::NodeMap<ListDigraph::Node> original(red);


    // Capacities of the MaxFlow reduction
    ListDigraph::ArcMap<int64_t> feasible_flow(red, 0); // From the approximation
    ListDigraph::ArcMap<int64_t> capacities(red, 0); // For the Max-Flow reduction, it is flow-demand


    ListDigraph::Node s = red.addNode();
    ListDigraph::Node t = red.addNode();


    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        v_in[v] = red.addNode();
        v_out[v] = red.addNode();
        original[v_in[v]] = v;
        original[v_out[v]] = v;

        ListDigraph::Arc split = red.addArc(v_in[v], v_out[v]);
        feasible_flow[split] = mu_v[v];
        capacities[split] = mu_v[v] - 1;

        ListDigraph::Arc sv = red.addArc(s, v_in[v]);
        ListDigraph::Arc vt = red.addArc(v_out[v], t);

        feasible_flow[sv] = starting_at[v];
        feasible_flow[vt] = ending_at[v];

        capacities[sv] = starting_at[v];
        capacities[vt] = ending_at[v];
    }

    for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        ListDigraph::Arc red_e = red.addArc(v_out[u], v_in[v]);
        feasible_flow[red_e] = mu[e];
        capacities[red_e] = mu[e];
    }

    ListDigraph::ArcMap<ListDigraph::Arc> direct(red);
    std::vector<ListDigraph::Arc> edges;
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        edges.push_back(e);
    }

    for (auto e : edges) {
        ListDigraph::Arc rev_e = red.addArc(red.target(e), red.source(e));
        capacities[rev_e] = countNodes(g);
        direct[e] = e;
        direct[rev_e] = e;
    }


    // Run Max-Flow algorithm
    EdmondsKarp<ListDigraph, ListDigraph::ArcMap<int64_t>> ek(red, capacities, s, t);

    // Set the flowMap to store the result in run
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ek.flowMap(flowMap);
    ek.run();

    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;

    // Compute the true flowMap (go back to Min-Flow) and remove 0 flow edges
    std::vector<ListDigraph::Arc> removal_list;

    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        ListDigraph::Arc d_e = direct[e];
        if (e != d_e) { // If it is a reverse edge
            // The flow on that edge is computed as the previous flow minus the one discounted by that edge,
            // plus the flow in the reverse direction (discounted in the Max-flow, therefore pushed in the Min-flow)
            flowMap[d_e] = feasible_flow[d_e] + flowMap[e] - flowMap[d_e];
            flowMap[e] = 0;

            removal_list.push_back(e);
            if (flowMap[d_e] == 0) {
                removal_list.push_back(d_e);
            }
        }
    }
    for (ListDigraph::Arc e : removal_list) {
        red.erase(e);
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        removal_list.clear();

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t && (path.empty() || path.back() != original[v]))
                path.push_back(original[v]);
            flowMap[e]--;
            if (flowMap[e] == 0) {
                removal_list.push_back(e);
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());

        path_cover.push_back(path);
        for (ListDigraph::Arc e : removal_list) {
            red.erase(e);
        }

        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }

    return path_cover;
}



std::vector<std::vector<ListDigraph::Node>> MPC(ListDigraph& g, std::vector<ListDigraph::Node>& S, std::vector<ListDigraph::Node>& T) {

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


    // Obtain Flow solution
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ns.flowMap(flowMap);

    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;

    // Remove 0 flow edges and st
    std::vector<ListDigraph::Arc> removal_list = {st};
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        if (flowMap[e] == 0) {
            removal_list.push_back(e);
        }
    }
    for (ListDigraph::Arc e : removal_list) {
        red.erase(e);
    }

    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        removal_list.clear();

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t && (path.empty() || path.back() != original[v]))
                path.push_back(original[v]);
            flowMap[e]--;
            if (flowMap[e] == 0) {
                removal_list.push_back(e);
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());

        path_cover.push_back(path);
        for (ListDigraph::Arc e : removal_list) {
            red.erase(e);
        }

        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }

    return path_cover;
}



std::vector<std::vector<ListDigraph::Node>> greedy_MPC(ListDigraph& g, std::vector<ListDigraph::Node>& S, std::vector<ListDigraph::Node>& T) {


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

    ListDigraph::NodeMap<ListDigraph::Node> v_in(g);
    ListDigraph::NodeMap<ListDigraph::Node> v_out(g);
    ListDigraph::NodeMap<ListDigraph::Node> original(red);

    // Capacities of the MaxFlow reduction
    ListDigraph::ArcMap<int64_t> capacities(red, 0); // From the approximation
    ListDigraph::ArcMap<int64_t> feasible_flow(red, 0); // For the Max-Flow reduction, it is flow-demand


    ListDigraph::Node s = red.addNode();
    ListDigraph::Node t = red.addNode();


    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        v_in[v] = red.addNode();
        v_out[v] = red.addNode();
        original[v_in[v]] = v;
        original[v_out[v]] = v;

        ListDigraph::Arc split = red.addArc(v_in[v], v_out[v]);
        feasible_flow[split] = mu_v[v];
        capacities[split] = mu_v[v] - 1;
    }
    for (ListDigraph::Node v : S) {
        ListDigraph::Arc sv = red.addArc(s, v_in[v]);
        feasible_flow[sv] = starting_at[v];
        capacities[sv] = starting_at[v];
    }
    for (ListDigraph::Node v : T) {
        ListDigraph::Arc vt = red.addArc(v_out[v], t);
        feasible_flow[vt] = ending_at[v];
        capacities[vt] = ending_at[v];
    }


    for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        ListDigraph::Arc red_e = red.addArc(v_out[u], v_in[v]);
        feasible_flow[red_e] = mu[e];
        capacities[red_e] = mu[e];
    }


    ListDigraph::ArcMap<ListDigraph::Arc> direct(red);
    std::vector<ListDigraph::Arc> edges;
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        edges.push_back(e);
    }

    for (auto e : edges) {
        ListDigraph::Arc rev_e = red.addArc(red.target(e), red.source(e));
        capacities[rev_e] = countNodes(g);
        direct[e] = e;
        direct[rev_e] = e;
    }

    // Run Max-Flow algorithm
    EdmondsKarp<ListDigraph, ListDigraph::ArcMap<int64_t>> ek(red, capacities, s, t);

    // Set the flowMap to store the result in run
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ek.flowMap(flowMap);
    ek.run();

    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;

    // Compute the true flowMap (go back to Min-Flow) and remove 0 flow edges
    std::vector<ListDigraph::Arc> removal_list;

    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        ListDigraph::Arc d_e = direct[e];
        if (e != d_e) {
            // The flow on that edge is computed as the previous flow minus the one discounted by that edge,
            // plus the flow in the reverse direction (discounted in the Max-flow, therefore pushed in the Min-flow)
            flowMap[d_e] = feasible_flow[d_e] + flowMap[e] - flowMap[d_e];
            flowMap[e] = 0;

            removal_list.push_back(e);
            if (flowMap[d_e] == 0) {
                removal_list.push_back(d_e);
            }
        }
    }
    for (ListDigraph::Arc e : removal_list) {
        red.erase(e);
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        removal_list.clear();

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t && (path.empty() || path.back() != original[v]))
                path.push_back(original[v]);
            flowMap[e]--;
            if (flowMap[e] == 0) {
                removal_list.push_back(e);
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());

        path_cover.push_back(path);
        for (ListDigraph::Arc e : removal_list) {
            red.erase(e);
        }

        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }

    return path_cover;
}



std::vector<std::vector<ListDigraph::Node>> U_MPC(ListDigraph& g, std::vector<ListDigraph::Node>& U) {

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


    // Obtain Flow solution
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ns.flowMap(flowMap);

    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;

    // Remove 0 flow edges and st
    std::vector<ListDigraph::Arc> removal_list = {st};
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        if (flowMap[e] == 0) {
            removal_list.push_back(e);
        }
    }
    for (ListDigraph::Arc e : removal_list) {
        red.erase(e);
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        removal_list.clear();

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t && (path.empty() || path.back() != original[v]))
                path.push_back(original[v]);
            flowMap[e]--;
            if (flowMap[e] == 0) {
                removal_list.push_back(e);
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());

        path_cover.push_back(path);
        for (ListDigraph::Arc e : removal_list) {
            red.erase(e);
        }

        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }

    return path_cover;
}



std::vector<std::vector<ListDigraph::Node>> greedy_U_MPC(ListDigraph& g, std::vector<ListDigraph::Node>& U) {

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

    ListDigraph::NodeMap<ListDigraph::Node> v_in(g);
    ListDigraph::NodeMap<ListDigraph::Node> v_out(g);
    ListDigraph::NodeMap<ListDigraph::Node> original(red);


    // Capacities of the MaxFlow reduction
    ListDigraph::ArcMap<int64_t> feasible_flow(red, 0); // From the approximation
    ListDigraph::ArcMap<int64_t> capacities(red, 0); // For the Max-Flow reduction, it is flow-demand


    ListDigraph::Node s = red.addNode();
    ListDigraph::Node t = red.addNode();


    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        v_in[v] = red.addNode();
        v_out[v] = red.addNode();
        original[v_in[v]] = v;
        original[v_out[v]] = v;

        ListDigraph::Arc split = red.addArc(v_in[v], v_out[v]);
        feasible_flow[split] = mu_v[v];
        capacities[split] = mu_v[v];
        if (in_U[v]) {
            capacities[split]--;
        }

        ListDigraph::Arc sv = red.addArc(s, v_in[v]);
        ListDigraph::Arc vt = red.addArc(v_out[v], t);

        feasible_flow[sv] = starting_at[v];
        feasible_flow[vt] = ending_at[v];

        capacities[sv] = starting_at[v];
        capacities[vt] = ending_at[v];
    }

    for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        ListDigraph::Arc red_e = red.addArc(v_out[u], v_in[v]);
        feasible_flow[red_e] = mu[e];
        capacities[red_e] = mu[e];
    }

    ListDigraph::ArcMap<ListDigraph::Arc> direct(red);
    std::vector<ListDigraph::Arc> edges;
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        edges.push_back(e);
    }

    for (auto e : edges) {
        ListDigraph::Arc rev_e = red.addArc(red.target(e), red.source(e));
        capacities[rev_e] = countNodes(g);
        direct[e] = e;
        direct[rev_e] = e;
    }


    // Run Max-Flow algorithm
    EdmondsKarp<ListDigraph, ListDigraph::ArcMap<int64_t>> ek(red, capacities, s, t);

    // Set the flowMap to store the result in run
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ek.flowMap(flowMap);
    ek.run();

    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;

    // Compute the true flowMap (go back to Min-Flow) and remove 0 flow edges
    std::vector<ListDigraph::Arc> removal_list;

    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        ListDigraph::Arc d_e = direct[e];
        if (e != d_e) {
            // The flow on that edge is computed as the previous flow minus the one discounted by that edge,
            // plus the flow in the reverse direction (discounted in the Max-flow, therefore pushed in the Min-flow)
            flowMap[d_e] = feasible_flow[d_e] + flowMap[e] - flowMap[d_e];
            flowMap[e] = 0;

            removal_list.push_back(e);
            if (flowMap[d_e] == 0) {
                removal_list.push_back(d_e);
            }
        }
    }
    for (ListDigraph::Arc e : removal_list) {
        red.erase(e);
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        removal_list.clear();

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t && (path.empty() || path.back() != original[v]))
                path.push_back(original[v]);
            flowMap[e]--;
            if (flowMap[e] == 0) {
                removal_list.push_back(e);
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());

        path_cover.push_back(path);
        for (ListDigraph::Arc e : removal_list) {
            red.erase(e);
        }

        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }

    return path_cover;
}



std::vector<std::vector<ListDigraph::Node>> U_MPC(ListDigraph& g, std::vector<ListDigraph::Node>& S, std::vector<ListDigraph::Node>& T, std::vector<ListDigraph::Node>& U) {

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


    // Obtain Flow solution
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ns.flowMap(flowMap);

    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;

    // Remove 0 flow edges and st
    std::vector<ListDigraph::Arc> removal_list = {st};
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        if (flowMap[e] == 0) {
            removal_list.push_back(e);
        }
    }
    for (ListDigraph::Arc e : removal_list) {
        red.erase(e);
    }

    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        removal_list.clear();

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t && (path.empty() || path.back() != original[v]))
                path.push_back(original[v]);
            flowMap[e]--;
            if (flowMap[e] == 0) {
                removal_list.push_back(e);
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());

        path_cover.push_back(path);
        for (ListDigraph::Arc e : removal_list) {
            red.erase(e);
        }

        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }

    return path_cover;
}



std::vector<std::vector<ListDigraph::Node>> greedy_U_MPC(ListDigraph& g, std::vector<ListDigraph::Node>& S, std::vector<ListDigraph::Node>& T, std::vector<ListDigraph::Node>& U) {

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

    ListDigraph::NodeMap<ListDigraph::Node> v_in(g);
    ListDigraph::NodeMap<ListDigraph::Node> v_out(g);
    ListDigraph::NodeMap<ListDigraph::Node> original(red);

    // Capacities of the MaxFlow reduction
    ListDigraph::ArcMap<int64_t> capacities(red, 0); // From the approximation
    ListDigraph::ArcMap<int64_t> feasible_flow(red, 0); // For the Max-Flow reduction, it is flow-demand


    ListDigraph::Node s = red.addNode();
    ListDigraph::Node t = red.addNode();


    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        v_in[v] = red.addNode();
        v_out[v] = red.addNode();
        original[v_in[v]] = v;
        original[v_out[v]] = v;

        ListDigraph::Arc split = red.addArc(v_in[v], v_out[v]);
        feasible_flow[split] = mu_v[v];
        capacities[split] = mu_v[v];
        if (in_U[v]) {
            capacities[split]--;
        }
    }
    for (ListDigraph::Node v : S) {
        ListDigraph::Arc sv = red.addArc(s, v_in[v]);
        feasible_flow[sv] = starting_at[v];
        capacities[sv] = starting_at[v];
    }
    for (ListDigraph::Node v : T) {
        ListDigraph::Arc vt = red.addArc(v_out[v], t);
        feasible_flow[vt] = ending_at[v];
        capacities[vt] = ending_at[v];
    }


    for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        ListDigraph::Arc red_e = red.addArc(v_out[u], v_in[v]);
        feasible_flow[red_e] = mu[e];
        capacities[red_e] = mu[e];
    }


    ListDigraph::ArcMap<ListDigraph::Arc> direct(red);
    std::vector<ListDigraph::Arc> edges;
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        edges.push_back(e);
    }

    for (auto e : edges) {
        ListDigraph::Arc rev_e = red.addArc(red.target(e), red.source(e));
        capacities[rev_e] = countNodes(g);
        direct[e] = e;
        direct[rev_e] = e;
    }

    // Run Max-Flow algorithm
    EdmondsKarp<ListDigraph, ListDigraph::ArcMap<int64_t>> ek(red, capacities, s, t);

    // Set the flowMap to store the result in run
    ListDigraph::ArcMap<int64_t> flowMap(red);
    ek.flowMap(flowMap);
    ek.run();

    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;

    // Compute the true flowMap (go back to Min-Flow) and remove 0 flow edges
    std::vector<ListDigraph::Arc> removal_list;

    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        ListDigraph::Arc d_e = direct[e];
        if (e != d_e) {
            // The flow on that edge is computed as the previous flow minus the one discounted by that edge,
            // plus the flow in the reverse direction (discounted in the Max-flow, therefore pushed in the Min-flow)
            flowMap[d_e] = feasible_flow[d_e] + flowMap[e] - flowMap[d_e];
            flowMap[e] = 0;

            removal_list.push_back(e);
            if (flowMap[d_e] == 0) {
                removal_list.push_back(d_e);
            }
        }
    }
    for (ListDigraph::Arc e : removal_list) {
        red.erase(e);
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        removal_list.clear();

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t && (path.empty() || path.back() != original[v]))
                path.push_back(original[v]);
            flowMap[e]--;
            if (flowMap[e] == 0) {
                removal_list.push_back(e);
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());

        path_cover.push_back(path);
        for (ListDigraph::Arc e : removal_list) {
            red.erase(e);
        }

        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }

    return path_cover;
}