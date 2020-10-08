#include <algorithms/top_sort.h>

#include <queue>
#include <stack>
#include <lemon/dfs.h>

using namespace lemon;



std::vector<ListDigraph::Node> topological_sort(ListDigraph& g) {

    ListDigraph::NodeMap<int64_t> in_degree_counter(g);
    for (ListDigraph::ArcIt e(g); e != INVALID; ++e)
        in_degree_counter[g.target(e)]++;

    std::vector<ListDigraph::Node> top_sort;
    std::queue<ListDigraph::Node> kahn_queue;

    for (ListDigraph::NodeIt s(g); s != INVALID; ++s) {
        if (in_degree_counter[s] == 0) {
            top_sort.push_back(s);
            kahn_queue.push(s);
        }
    }


    while (!kahn_queue.empty()) {
        ListDigraph::Node u = kahn_queue.front();
        kahn_queue.pop();
        for (ListDigraph::OutArcIt e(g, u); e != INVALID; ++e) {
            ListDigraph::Node v = g.target(e);
            in_degree_counter[v]--;
            if (in_degree_counter[v] == 0) {
                kahn_queue.push(v);
                top_sort.push_back(v);
            }
        }
    }

    return top_sort;
}



std::vector<ListDigraph::Node> topological_sort(ListDigraph& g, std::vector<ListDigraph::Node>& S) {

    ListDigraph::NodeMap<int64_t> in_degree_counter(g);
    for (ListDigraph::ArcIt e(g); e != INVALID; ++e)
        in_degree_counter[g.target(e)]++;

    std::vector<ListDigraph::Node> top_sort;
    std::queue<ListDigraph::Node> kahn_queue;

    for (ListDigraph::Node s : S) {
        top_sort.push_back(s);
        kahn_queue.push(s);
    }

    while (!kahn_queue.empty()) {
        ListDigraph::Node u = kahn_queue.front();
        kahn_queue.pop();
        for (ListDigraph::OutArcIt e(g, u); e != INVALID; ++e) {
            ListDigraph::Node v = g.target(e);
            in_degree_counter[v]--;
            if (in_degree_counter[v] == 0) {
                kahn_queue.push(v);
                top_sort.push_back(v);
            }
        }
    }

    return top_sort;
}



std::vector<ListDigraph::Node> topological_sort_dfs(ListDigraph& g) {

    Dfs<ListDigraph> dfs(g);
    std::stack<ListDigraph::Node> dfs_stack;
    ListDigraph::NodeMap<bool> added(g, false);
    std::vector<ListDigraph::Node> post_order_dfs;
    dfs.init();


    for (ListDigraph::NodeIt s(g); s != INVALID; ++s) {
        if (countInArcs(g, s) == 0) {
            dfs_stack.push(s);
            dfs.addSource(s);
        }
    }


    while (!dfs.emptyQueue()) {
        ListDigraph::Arc e = dfs.processNextArc();
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        while (dfs_stack.top() != u) {
            if (!added[dfs_stack.top()]) {
                post_order_dfs.push_back(dfs_stack.top());
                added[dfs_stack.top()] = true;
            }
            dfs_stack.pop();
        }
        dfs_stack.push(v);
    }


    while (!dfs_stack.empty()) {
        if (!added[dfs_stack.top()]) {
            post_order_dfs.push_back(dfs_stack.top());
            added[dfs_stack.top()] = true;
        }
        dfs_stack.pop();
    }

    std::reverse(post_order_dfs.begin(), post_order_dfs.end());
    return post_order_dfs;
}



std::vector<ListDigraph::Node> topological_sort_dfs(ListDigraph& g, std::vector<ListDigraph::Node>& S) {

    std::stack<ListDigraph::Node> dfs_stack;

    ListDigraph::NodeMap<bool> added(g, false);
    std::vector<ListDigraph::Node> post_order_dfs;



    Dfs<ListDigraph> dfs(g);
    dfs.init();
    for (ListDigraph::Node s : S) {
        dfs_stack.push(s);
        dfs.addSource(s);
    }
    while (!dfs.emptyQueue()) {
        ListDigraph::Arc e = dfs.processNextArc();
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        while (dfs_stack.top() != u) {
            if (!added[dfs_stack.top()]) {
                post_order_dfs.push_back(dfs_stack.top());
                added[dfs_stack.top()] = true;
            }
            dfs_stack.pop();
        }
        dfs_stack.push(v);
    }


    while (!dfs_stack.empty()) {
        if (!added[dfs_stack.top()]) {
            post_order_dfs.push_back(dfs_stack.top());
            added[dfs_stack.top()] = true;
        }
        dfs_stack.pop();
    }

    std::reverse(post_order_dfs.begin(), post_order_dfs.end());
    return post_order_dfs;
}