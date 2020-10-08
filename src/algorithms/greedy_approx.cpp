#include <algorithms/greedy_approx.h>
#include <algorithms/top_sort.h>

using namespace lemon;



std::vector<std::vector<ListDigraph::Node>> greedy_approximation_MPC(ListDigraph& g) {

    std::vector<std::vector<ListDigraph::Node>> path_cover;
    std::vector<ListDigraph::Node> top_order = topological_sort(g);

    int64_t n = countNodes(g);
    int64_t covered_vertices = 0;
    ListDigraph::NodeMap<bool> covered(g, false);

    while (covered_vertices < n) {
        int64_t max_length = 0;
        ListDigraph::Node max_length_source = top_order[0];

        ListDigraph::NodeMap<int64_t> max_length_from_vertex(g);
        ListDigraph::NodeMap<ListDigraph::Node> next_max_length_from_vertex(g);
        for (ListDigraph::Node v : top_order) {
            max_length_from_vertex[v] = ((covered[v]) ? 0 : 1);
            next_max_length_from_vertex[v] = v;
        }

        for (std::vector<ListDigraph::Node>::reverse_iterator rit = top_order.rbegin(); rit != top_order.rend(); ++rit) {
            ListDigraph::Node u = *rit;
            int64_t max = 0;
            for (ListDigraph::OutArcIt e(g, u); e != INVALID; ++e) {
                ListDigraph::Node v = g.target(e);
                int64_t max_v = max_length_from_vertex[v];
                if (max_v > max) {
                    max = max_v;
                    next_max_length_from_vertex[u] = v;
                }
            }

            max_length_from_vertex[u] = max + max_length_from_vertex[u];
            if (max_length_from_vertex[u] > max_length) {
                max_length = max_length_from_vertex[u];
                max_length_source = u;
            }
        }

        std::vector<ListDigraph::Node> path;
        ListDigraph::Node v = max_length_source;
        while (true) {
            if (!covered[v]) {
                covered[v] = true;
                ++covered_vertices;
            }
            path.push_back(v);
            if (v == next_max_length_from_vertex[v]) break;
            v = next_max_length_from_vertex[v];
        }

        path_cover.push_back(path);
    }

    return path_cover;
}



std::vector<std::vector<ListDigraph::Arc>> greedy_approximation_MPC_edges(ListDigraph& g) {

    std::vector<std::vector<ListDigraph::Arc>> path_cover;
    std::vector<ListDigraph::Node> top_order = topological_sort(g);

    int64_t n = countNodes(g);
    int64_t covered_vertices = 0;
    ListDigraph::NodeMap<bool> covered(g, false);

    while (covered_vertices < n) {
        int64_t max_length = 0;
        ListDigraph::Node max_length_source = top_order[0];

        ListDigraph::NodeMap<int64_t> max_length_from_vertex(g);
        ListDigraph::NodeMap<ListDigraph::Node> next_max_length_from_vertex(g);
        ListDigraph::NodeMap<ListDigraph::Arc> next_max_length_from_edges(g);

        for (ListDigraph::Node v : top_order) {
            max_length_from_vertex[v] = ((covered[v]) ? 0 : 1);
            next_max_length_from_vertex[v] = v;
        }

        for (std::vector<ListDigraph::Node>::reverse_iterator rit = top_order.rbegin(); rit != top_order.rend(); ++rit) {
            ListDigraph::Node u = *rit;
            int64_t max = 0;
            for (ListDigraph::OutArcIt e(g, u); e != INVALID; ++e) {
                ListDigraph::Node v = g.target(e);
                int64_t max_v = max_length_from_vertex[v];
                if (max_v > max) {
                    max = max_v;
                    next_max_length_from_vertex[u] = v;
                    next_max_length_from_edges[u] = e;
                }
            }

            max_length_from_vertex[u] = max + max_length_from_vertex[u];
            if (max_length_from_vertex[u] > max_length) {
                max_length = max_length_from_vertex[u];
                max_length_source = u;
            }
        }

        std::vector<ListDigraph::Arc> path;
        ListDigraph::Node v = max_length_source;
        ListDigraph::Arc e(INVALID);
        while (true) {
            if (!covered[v]) {
                covered[v] = true;
                ++covered_vertices;
            }
            if (g.valid(next_max_length_from_edges[v]))
                path.push_back(next_max_length_from_edges[v]);
            if (v == next_max_length_from_vertex[v]) break;
            v = next_max_length_from_vertex[v];
            e = next_max_length_from_edges[v];
        }

        path_cover.push_back(path);
    }

    return path_cover;
}



std::vector<std::vector<ListDigraph::Node>> greedy_approximation_MPC(ListDigraph& g, std::vector<ListDigraph::Node>& S, std::vector<ListDigraph::Node>& T) {

    std::vector<std::vector<ListDigraph::Node>> path_cover;
    std::vector<ListDigraph::Node> top_order = topological_sort(g);

    int64_t n = countNodes(g);
    int64_t covered_vertices = 0;
    ListDigraph::NodeMap<bool> covered(g, false);

    ListDigraph::NodeMap<bool> in_T(g, false);
    ListDigraph::NodeMap<bool> in_S(g, false);
    for (ListDigraph::Node v : T) {
        in_T[v] = true;
    }
    for (ListDigraph::Node v : S) {
        in_S[v] = true;
    }


    while (covered_vertices < n) {
        int64_t max_length = 0;
        ListDigraph::Node max_length_source = top_order[0];

        ListDigraph::NodeMap<int64_t> max_length_from_vertex(g);
        ListDigraph::NodeMap<ListDigraph::Node> next_max_length_from_vertex(g);
        for (ListDigraph::Node v : top_order) {
            max_length_from_vertex[v] = ((covered[v]) ? 0 : (in_T[v]) ? 1 : 0);
            next_max_length_from_vertex[v] = v;
        }

        for (std::vector<ListDigraph::Node>::reverse_iterator rit = top_order.rbegin(); rit != top_order.rend(); ++rit) {
            ListDigraph::Node u = *rit;
            int64_t max = -1;
            for (ListDigraph::OutArcIt e(g, u); e != INVALID; ++e) {
                ListDigraph::Node v = g.target(e);
                int64_t max_v = max_length_from_vertex[v];
                if (max_v > max) {
                    max = max_v;
                    next_max_length_from_vertex[u] = v;
                }
            }

            max_length_from_vertex[u] = ((covered[u] ? 0 : 1) + max > max_length_from_vertex[u]) ? (covered[u] ? 0 : 1) + max : max_length_from_vertex[u];
            if (in_T[u] && ((covered[u] ? 0 : 1) == max_length_from_vertex[u])) {
                next_max_length_from_vertex[u] = u;
            }
            if (max_length_from_vertex[u] > max_length && in_S[u]) {
                max_length = max_length_from_vertex[u];
                max_length_source = u;
            }
        }

        std::vector<ListDigraph::Node> path;
        ListDigraph::Node v = max_length_source;
        while (true) {
            if (!covered[v]) {
                covered[v] = true;
                ++covered_vertices;
            }
            path.push_back(v);
            if (v == next_max_length_from_vertex[v]) break;
            v = next_max_length_from_vertex[v];
        }

        path_cover.push_back(path);
    }

    return path_cover;
}



std::vector<std::vector<ListDigraph::Arc>> greedy_approximation_MPC_edges(ListDigraph& g, std::vector<ListDigraph::Node>& S, std::vector<ListDigraph::Node>& T) {

    ListDigraph::Arc null_edge(INVALID);
    std::vector<std::vector<ListDigraph::Arc>> path_cover;
    std::vector<ListDigraph::Node> top_order = topological_sort(g);

    int64_t n = countNodes(g);
    int64_t covered_vertices = 0;
    ListDigraph::NodeMap<bool> covered(g, false);

    ListDigraph::NodeMap<bool> in_T(g, false);
    ListDigraph::NodeMap<bool> in_S(g, false);

    for (ListDigraph::Node v : T) {
        in_T[v] = true;
    }
    for (ListDigraph::Node v : S) {
        in_S[v] = true;
    }


    while (covered_vertices < n) {
        int64_t max_length = 0;
        ListDigraph::Node max_length_source = top_order[0];

        ListDigraph::NodeMap<int64_t> max_length_from_vertex(g);
        ListDigraph::NodeMap<ListDigraph::Node> next_max_length_from_vertex(g);
        ListDigraph::NodeMap<ListDigraph::Arc> next_max_length_from_edges(g, null_edge);

        for (ListDigraph::Node v : top_order) {
            max_length_from_vertex[v] = ((covered[v]) ? 0 : (in_T[v]) ? 1 : 0);
            next_max_length_from_vertex[v] = v;
        }

        for (std::vector<ListDigraph::Node>::reverse_iterator rit = top_order.rbegin(); rit != top_order.rend(); ++rit) {
            ListDigraph::Node u = *rit;
            int64_t max = -1;
            for (ListDigraph::OutArcIt e(g, u); e != INVALID; ++e) {
                ListDigraph::Node v = g.target(e);
                int64_t max_v = max_length_from_vertex[v];
                if (max_v > max) {
                    max = max_v;
                    next_max_length_from_vertex[u] = v;
                    next_max_length_from_edges[u] = e;
                }
            }

            max_length_from_vertex[u] = ((covered[u] ? 0 : 1) + max > max_length_from_vertex[u]) ? (covered[u] ? 0 : 1) + max : max_length_from_vertex[u];
            if (in_T[u] && ((covered[u] ? 0 : 1) == max_length_from_vertex[u])) {
                next_max_length_from_vertex[u] = u;
                next_max_length_from_edges[u] = null_edge;
            }
            if (max_length_from_vertex[u] > max_length && in_S[u]) {
                max_length = max_length_from_vertex[u];
                max_length_source = u;
            }
        }

        std::vector<ListDigraph::Arc> path;
        ListDigraph::Node v = max_length_source;
        ListDigraph::Arc e = null_edge;
        while (true) {
            if (!covered[v]) {
                covered[v] = true;
                ++covered_vertices;
            }
            if (g.valid(next_max_length_from_edges[v]))
                path.push_back(next_max_length_from_edges[v]);
            if (v == next_max_length_from_vertex[v]) break;
            v = next_max_length_from_vertex[v];
            e = next_max_length_from_edges[v];
        }

        path_cover.push_back(path);
    }

    return path_cover;
}



std::vector<std::vector<ListDigraph::Node>> greedy_approximation_U_MPC(ListDigraph& g, std::vector<ListDigraph::Node>& U) {

    std::vector<std::vector<ListDigraph::Node>> path_cover;
    std::vector<ListDigraph::Node> top_order = topological_sort(g);

    int64_t u = U.size();
    int64_t covered_vertices = 0;

    // Set as false only the vertices in U
    ListDigraph::NodeMap<bool> covered(g, true);
    for (ListDigraph::Node v : U) {
        covered[v] = false;
    }

    while (covered_vertices < u) {
        int64_t max_length = 0;
        ListDigraph::Node max_length_source = top_order[0];

        ListDigraph::NodeMap<int64_t> max_length_from_vertex(g);
        ListDigraph::NodeMap<ListDigraph::Node> next_max_length_from_vertex(g);
        for (ListDigraph::Node v : top_order) {
            max_length_from_vertex[v] = ((covered[v]) ? 0 : 1);
            next_max_length_from_vertex[v] = v;
        }

        for (std::vector<ListDigraph::Node>::reverse_iterator rit = top_order.rbegin(); rit != top_order.rend(); ++rit) {
            ListDigraph::Node u = *rit;
            int64_t max = 0;
            for (ListDigraph::OutArcIt e(g, u); e != INVALID; ++e) {
                ListDigraph::Node v = g.target(e);
                int64_t max_v = max_length_from_vertex[v];
                if (max_v > max) {
                    max = max_v;
                    next_max_length_from_vertex[u] = v;
                }
            }

            max_length_from_vertex[u] = max + max_length_from_vertex[u];
            if (max_length_from_vertex[u] > max_length) {
                max_length = max_length_from_vertex[u];
                max_length_source = u;
            }
        }

        std::vector<ListDigraph::Node> path;
        ListDigraph::Node v = max_length_source;
        while (true) {
            if (!covered[v]) {
                covered[v] = true;
                ++covered_vertices;
            }
            path.push_back(v);
            if (v == next_max_length_from_vertex[v]) break;
            v = next_max_length_from_vertex[v];
        }

        path_cover.push_back(path);
    }

    return path_cover;
}



std::vector<std::vector<ListDigraph::Arc>> greedy_approximation_U_MPC_edges(ListDigraph& g, std::vector<ListDigraph::Node>& U) {

    ListDigraph::Arc null_edge(INVALID);
    std::vector<std::vector<ListDigraph::Arc>> path_cover;
    std::vector<ListDigraph::Node> top_order = topological_sort(g);

    int64_t u = U.size();
    int64_t covered_vertices = 0;
    // Set as false only the vertices in U
    ListDigraph::NodeMap<bool> covered(g, true);
    for (ListDigraph::Node v : U) {
        covered[v] = false;
    }

    while (covered_vertices < u) {
        int64_t max_length = 0;
        ListDigraph::Node max_length_source = top_order[0];

        ListDigraph::NodeMap<int64_t> max_length_from_vertex(g);
        ListDigraph::NodeMap<ListDigraph::Node> next_max_length_from_vertex(g);
        ListDigraph::NodeMap<ListDigraph::Arc> next_max_length_from_edges(g, null_edge);

        for (ListDigraph::Node v : top_order) {
            max_length_from_vertex[v] = ((covered[v]) ? 0 : 1);
            next_max_length_from_vertex[v] = v;
        }

        for (std::vector<ListDigraph::Node>::reverse_iterator rit = top_order.rbegin(); rit != top_order.rend(); ++rit) {
            ListDigraph::Node u = *rit;
            int64_t max = 0;
            for (ListDigraph::OutArcIt e(g, u); e != INVALID; ++e) {
                ListDigraph::Node v = g.target(e);
                int64_t max_v = max_length_from_vertex[v];
                if (max_v > max) {
                    max = max_v;
                    next_max_length_from_vertex[u] = v;
                    next_max_length_from_edges[u] = e;
                }
            }

            max_length_from_vertex[u] = max + max_length_from_vertex[u];
            if (max_length_from_vertex[u] > max_length) {
                max_length = max_length_from_vertex[u];
                max_length_source = u;
            }
        }

        std::vector<ListDigraph::Arc> path;
        ListDigraph::Node v = max_length_source;
        ListDigraph::Arc e = null_edge;
        while (true) {
            if (!covered[v]) {
                covered[v] = true;
                ++covered_vertices;
            }
            if (g.valid(next_max_length_from_edges[v]))
                path.push_back(next_max_length_from_edges[v]);
            if (v == next_max_length_from_vertex[v]) break;
            v = next_max_length_from_vertex[v];
            e = next_max_length_from_edges[v];
        }

        path_cover.push_back(path);
    }

    return path_cover;
}



std::vector<std::vector<ListDigraph::Node>> greedy_approximation_U_MPC(ListDigraph& g, std::vector<ListDigraph::Node>& S, std::vector<ListDigraph::Node>& T, std::vector<ListDigraph::Node>& U) {

    std::vector<std::vector<ListDigraph::Node>> path_cover;
    std::vector<ListDigraph::Node> top_order = topological_sort(g);

    int64_t u = U.size();
    int64_t covered_vertices = 0;
    // Set as false only the vertices in U
    ListDigraph::NodeMap<bool> covered(g, true);
    for (ListDigraph::Node v : U) {
        covered[v] = false;
    }

    ListDigraph::NodeMap<bool> in_T(g, false);
    ListDigraph::NodeMap<bool> in_S(g, false);
    for (ListDigraph::Node v : T) {
        in_T[v] = true;
    }
    for (ListDigraph::Node v : S) {
        in_S[v] = true;
    }


    while (covered_vertices < u) {
        int64_t max_length = 0;
        ListDigraph::Node max_length_source = top_order[0];

        ListDigraph::NodeMap<int64_t> max_length_from_vertex(g);
        ListDigraph::NodeMap<ListDigraph::Node> next_max_length_from_vertex(g);
        for (ListDigraph::Node v : top_order) {
            max_length_from_vertex[v] = ((covered[v]) ? 0 : (in_T[v]) ? 1 : 0);
            next_max_length_from_vertex[v] = v;
        }

        for (std::vector<ListDigraph::Node>::reverse_iterator rit = top_order.rbegin(); rit != top_order.rend(); ++rit) {
            ListDigraph::Node u = *rit;
            int64_t max = -1;
            for (ListDigraph::OutArcIt e(g, u); e != INVALID; ++e) {
                ListDigraph::Node v = g.target(e);
                int64_t max_v = max_length_from_vertex[v];
                if (max_v > max) {
                    max = max_v;
                    next_max_length_from_vertex[u] = v;
                }
            }

            max_length_from_vertex[u] = ((covered[u] ? 0 : 1) + max > max_length_from_vertex[u]) ? (covered[u] ? 0 : 1) + max : max_length_from_vertex[u];
            if (in_T[u] && ((covered[u] ? 0 : 1) == max_length_from_vertex[u])) {
                next_max_length_from_vertex[u] = u;
            }
            if (max_length_from_vertex[u] > max_length && in_S[u]) {
                max_length = max_length_from_vertex[u];
                max_length_source = u;
            }
        }

        std::vector<ListDigraph::Node> path;
        ListDigraph::Node v = max_length_source;
        while (true) {
            if (!covered[v]) {
                covered[v] = true;
                ++covered_vertices;
            }
            path.push_back(v);
            if (v == next_max_length_from_vertex[v]) break;
            v = next_max_length_from_vertex[v];
        }

        path_cover.push_back(path);
    }

    return path_cover;
}



std::vector<std::vector<ListDigraph::Arc>> greedy_approximation_U_MPC_edges(ListDigraph& g, std::vector<ListDigraph::Node>& S, std::vector<ListDigraph::Node>& T, std::vector<ListDigraph::Node>& U) {

    ListDigraph::Arc null_edge(INVALID);
    std::vector<std::vector<ListDigraph::Arc>> path_cover;
    std::vector<ListDigraph::Node> top_order = topological_sort(g);

    int64_t u = U.size();
    int64_t covered_vertices = 0;
    // Set as false only the vertices in U
    ListDigraph::NodeMap<bool> covered(g, true);
    for (ListDigraph::Node v : U) {
        covered[v] = false;
    }

    ListDigraph::NodeMap<bool> in_T(g, false);
    ListDigraph::NodeMap<bool> in_S(g, false);
    for (ListDigraph::Node v : T) {
        in_T[v] = true;
    }
    for (ListDigraph::Node v : S) {
        in_S[v] = true;
    }


    while (covered_vertices < u) {
        int64_t max_length = 0;
        ListDigraph::Node max_length_source = top_order[0];

        ListDigraph::NodeMap<int64_t> max_length_from_vertex(g);
        ListDigraph::NodeMap<ListDigraph::Node> next_max_length_from_vertex(g);
        ListDigraph::NodeMap<ListDigraph::Arc> next_max_length_from_edges(g, null_edge);

        for (ListDigraph::Node v : top_order) {
            max_length_from_vertex[v] = ((covered[v]) ? 0 : (in_T[v]) ? 1 : 0);
            next_max_length_from_vertex[v] = v;
        }

        for (std::vector<ListDigraph::Node>::reverse_iterator rit = top_order.rbegin(); rit != top_order.rend(); ++rit) {
            ListDigraph::Node u = *rit;
            int64_t max = -1;
            for (ListDigraph::OutArcIt e(g, u); e != INVALID; ++e) {
                ListDigraph::Node v = g.target(e);
                int64_t max_v = max_length_from_vertex[v];
                if (max_v > max) {
                    max = max_v;
                    next_max_length_from_vertex[u] = v;
                    next_max_length_from_edges[u] = e;
                }
            }

            max_length_from_vertex[u] = ((covered[u] ? 0 : 1) + max > max_length_from_vertex[u]) ? (covered[u] ? 0 : 1) + max : max_length_from_vertex[u];
            if (in_T[u] && ((covered[u] ? 0 : 1) == max_length_from_vertex[u])) {
                next_max_length_from_vertex[u] = u;
                next_max_length_from_edges[u] = null_edge;
            }
            if (max_length_from_vertex[u] > max_length && in_S[u]) {
                max_length = max_length_from_vertex[u];
                max_length_source = u;
            }
        }

        std::vector<ListDigraph::Arc> path;
        ListDigraph::Node v = max_length_source;
        ListDigraph::Arc e = null_edge;
        while (true) {
            if (!covered[v]) {
                covered[v] = true;
                ++covered_vertices;
            }
            if (g.valid(next_max_length_from_edges[v]))
                path.push_back(next_max_length_from_edges[v]);
            if (v == next_max_length_from_vertex[v]) break;
            v = next_max_length_from_vertex[v];
            e = next_max_length_from_edges[v];
        }

        path_cover.push_back(path);
    }

    return path_cover;
}