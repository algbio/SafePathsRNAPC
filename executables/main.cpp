#include <iostream>
#include <algorithms/greedy_approx.h>
#include <algorithms/mpc.h>
#include <algorithms/safe_paths.h>
#include <algorithms/filter_paths.h>
#include <utils.h>

#include <lemon/cost_scaling.h>
#include <lemon/edmonds_karp.h>


#include <chrono>


using namespace lemon;


int main(int argc, char *argv[]) {
    //ListDigraph g;




/*
    ListDigraph::Node x0 = g.addNode();
    ListDigraph::Node x1 = g.addNode();
    ListDigraph::Node x2 = g.addNode();

    g.addArc(x0, x1);
    g.addArc(x1, x2);
*/


/*
    ListDigraph::Node x0 = g.addNode();
    ListDigraph::Node x1 = g.addNode();
    ListDigraph::Node x2 = g.addNode();

    g.addArc(x0, x1);
    g.addArc(x0, x2);
*/
/*
    ListDigraph::Node x0 = g.addNode();
    ListDigraph::Node x1 = g.addNode();
    ListDigraph::Node x2 = g.addNode();
    ListDigraph::Node x3 = g.addNode();
    ListDigraph::Node x4 = g.addNode();
    ListDigraph::Node x5 = g.addNode();



    g.addArc(x0, x1);
    g.addArc(x0, x4);
    g.addArc(x0, x5);
    g.addArc(x1, x2);
    g.addArc(x2, x3);
    g.addArc(x2, x4);
    g.addArc(x4, x5);
*/

/*
    ListDigraph::Node x0 = g.addNode();
    ListDigraph::Node x1 = g.addNode();
    ListDigraph::Node x2 = g.addNode();
    ListDigraph::Node x3 = g.addNode();
    ListDigraph::Node x4 = g.addNode();

    g.addArc(x0, x1);
    g.addArc(x0, x3);
    g.addArc(x1, x2);
    g.addArc(x2, x3);
    g.addArc(x3, x4);
*/

     // EdmonsKarp changes the greedy solution, even when that is already of minimum size
/*
    ListDigraph::Node x0 = g.addNode();
    ListDigraph::Node x1 = g.addNode();
    ListDigraph::Node x2 = g.addNode();
    ListDigraph::Node x3 = g.addNode();
    ListDigraph::Node x4 = g.addNode();
    ListDigraph::Node x5 = g.addNode();
    ListDigraph::Node x6 = g.addNode();
    ListDigraph::Node x7 = g.addNode();

    g.addArc(x0, x1);
    g.addArc(x0, x2);
    g.addArc(x1, x3);
    g.addArc(x2, x3);
    g.addArc(x3, x4);
    g.addArc(x4, x5);
    g.addArc(x4, x6);
    g.addArc(x5, x7);
    g.addArc(x6, x7);
*/

/*
    ListDigraph::Node x0 = g.addNode();
    ListDigraph::Node x1 = g.addNode();
    ListDigraph::Node x2 = g.addNode();
    ListDigraph::Node x3 = g.addNode();


    g.addArc(x0, x1);
    g.addArc(x0, x2);
    g.addArc(x2, x3);
*/
    /*
    ListDigraph::Node x0 = g.addNode();
    ListDigraph::Node x1 = g.addNode();
    ListDigraph::Node x2 = g.addNode();
    ListDigraph::Node x3 = g.addNode();
    ListDigraph::Node x4 = g.addNode();
    ListDigraph::Node x5 = g.addNode();


    g.addArc(x0, x1);
    g.addArc(x0, x2);
    g.addArc(x0, x3);
    g.addArc(x1, x4);
    g.addArc(x2, x4);
    g.addArc(x3, x5);
    */

/*
    ListDigraph::Node x0 = g.addNode();
    ListDigraph::Node x1 = g.addNode();
    ListDigraph::Node x2 = g.addNode();
    ListDigraph::Node x3 = g.addNode();
    ListDigraph::Node x4 = g.addNode();
    ListDigraph::Node x5 = g.addNode();
    ListDigraph::Node x6 = g.addNode();
    ListDigraph::Node x7 = g.addNode();

    g.addArc(x0, x1);
    g.addArc(x0, x2);
    g.addArc(x0, x3);
    g.addArc(x1, x4);
    g.addArc(x2, x4);
    g.addArc(x2, x5);
    g.addArc(x2, x6);
    g.addArc(x3, x6);
    g.addArc(x4, x7);
    g.addArc(x5, x7);
    g.addArc(x6, x7);
*/

    //DigraphWriter<ListDigraph>(g, std::cout).run();
    /*
    for (ListDigraph::NodeIt v(g); v!=INVALID; ++v) {
        std::cout << g.id(v) << std::endl;
    }

    std::cout << "We have a directed graph with " << countNodes(g) << " nodes "
         << "and " << countArcs(g) << " arc." << std::endl;
*/

    /*
    Bfs<ListDigraph> bfs(g);
    bfs.init();
    bfs.addSource(x0);

    std::cout << "BFS" << std::endl;
    while (!bfs.emptyQueue()) {
        ListDigraph::Node v = bfs.processNextNode();
        std::cout << g.id(v) << ",";
    }
    std::cout << std::endl;*/
    /*bfs.run(u);

    for (ListDigraph::NodeIt n(g); n != INVALID; ++n) {
        if (bfs.reached(n)) {
            std::cout << g.id(n);
            ListDigraph::Node pred = bfs.predNode(n);
            while (pred != INVALID) {
                std::cout << "<-" << g.id(pred);
                pred = bfs.predNode(pred);
            }
            std::cout << std::endl;
        }
    }*/
/*
    std::cout << "TopSort Kahn's Algorithm" << std::endl;
    std::vector<ListDigraph::Node> top = topological_sort(g, {x0});
    for (ListDigraph::Node v : top) {
        std::cout << g.id(v) << " , ";
    }
    std::cout << std::endl;

    std::cout << "TopSort DFS Algorithm" << std::endl;
    std::vector<ListDigraph::Node> top_dfs = topological_sort_dfs(g, {x0});
    for (ListDigraph::Node v : top_dfs) {
        std::cout << g.id(v) << " , ";
    }
    std::cout << std::endl;
*/

    /*for (ListDigraph::NodeIt n(g); n != INVALID; ++n) {
        if (dfs.reached(n)) {
            std::cout << g.id(n);
            ListDigraph::Node pred = dfs.predNode(n);
            while (pred != INVALID) {
                std::cout << "<-" << g.id(pred);
                pred = dfs.predNode(pred);
            }
            std::cout << std::endl;
        }
    }*/

    /*
    ListDigraph::NodeMap<std::string> label(g);
    label[u] = "u label";
    label[v] = "v label";

    ListDigraph::ArcMap<int64_t> cost(g);
    cost[a] = 16;
    */
/*
    std::cout << "Greedy Approximation MPC" << std::endl;
    for (auto path: greedy_approximation_MPC(g)) {
        for (auto vertex: path) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;
    }
*/
/*
    std::cout << "Greedy Approximation Edges MPC" << std::endl;
    for (auto path: greedy_approximation_MPC_edges(g)) {
        for (auto edge: path) {
            std::cout << g.id(g.source(edge)) << " --> " << g.id(g.target(edge)) << " , ";
        }
        std::cout << std::endl;
    }
*/
    /*
    std::cout << "Greedy Approximation MPC from x_0 to x1,x2" << std::endl;
    for (auto path: greedy_approximation_MPC(g, {x0}, {x1, x2})) {
        for (auto vertex: path) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;
    }
    */
/*
    std::cout << "MPC Min-Flow<MinCostFlow> reduction" << std::endl;
    for (auto path: MPC(g)) {
        for (auto vertex: path) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;
    }
*/
/*
    std::vector<ListDigraph::Node> S = {x0, x1};
    std::vector<ListDigraph::Node> T = {x4, x5};

    std::cout << "MPC Min-Flow reduction" << std::endl;
    for (auto path: MPC(g, S, T)) {
        for (auto vertex: path) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;
    }


    std::cout << "MPC Min-Flow<Greedy+MaxFlow> reduction" << std::endl;
    for (auto path: greedy_MPC(g, S, T)) {
        for (auto vertex: path) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;
    }
*/
/*
    std::cout << "Path Maximal safe paths MPC" << std::endl;
    for (auto pair: path_maximal_safe_paths_MPC(g)) {
        std::cout << "Path: ";
        for (auto vertex: pair.first) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;

        for (auto path: pair.second) {
            for (auto vertex: path) {
                std::cout << g.id(vertex) << " , ";
            }
            std::cout << std::endl;
        }
    }
*/
/*
    std::cout << "Path Maximal safe paths PC, l = 2 (when width = 1)" << std::endl;
    for (auto pair: path_maximal_safe_paths_PC(g, 2)) {
        std::cout << "Path: ";
        for (auto vertex: pair.first) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;

        for (auto path: pair.second) {
            for (auto vertex: path) {
                std::cout << g.id(vertex) << " , ";
            }
            std::cout << std::endl;
        }
    }
*/

    /*
    std::cout << "Path Maximal safe paths MPC from x_0 to x1,x2" << std::endl;
    for (auto pair: path_maximal_safe_paths_MPC(g,{x0}, {x1, x2})) {
        std::cout << "Path: ";
        for (auto vertex: pair.first) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;

        for (auto path: pair.second) {
            for (auto vertex: path) {
                std::cout << g.id(vertex) << " , ";
            }
            std::cout << std::endl;
        }
    }*/

/*
    std::cout << "Path Maximal safe paths PC from x_0 to x2, l = 2 (width = 1)" << std::endl;
    for (auto pair: path_maximal_safe_paths_MPC(g,{x0}, {x2}, 2)) {
        std::cout << "Path: ";
        for (auto vertex: pair.first) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;

        for (auto path: pair.second) {
            for (auto vertex: path) {
                std::cout << g.id(vertex) << " , ";
            }
            std::cout << std::endl;
        }
    }
    */

/*
    std::cout << "U_MPC Min-Flow<MinCostFlow> reduction" << std::endl;
    for (auto path: U_MPC(g, {x0, x1})) {
        for (auto vertex: path) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;
    }
*/


/*
    std::cout << "Greedy Approximation U_MPC" << std::endl;
    for (auto path: greedy_approximation_U_MPC(g, {x0, x1})) {
        for (auto vertex: path) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;
    }
*/

/*
    std::cout << "Greedy Approximation Edges U_MPC" << std::endl;
    for (auto path: greedy_approximation_U_MPC_edges(g, {x0, x1})) {
        for (auto edge: path) {
            std::cout << g.id(g.source(edge)) << " --> " << g.id(g.target(edge)) << " , ";
        }
        std::cout << std::endl;
    }
*/

/*
    std::cout << "U_MPC Min-Flow<Greedy-MaxFlow> reduction" << std::endl;
    for (auto path: greedy_U_MPC(g, {x0, x1})) {
        for (auto vertex: path) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;
    }
*/
/*
    std::vector<ListDigraph::Node> U = {x0, x2, x3, x4};

    std::cout << "U_MPC Min-Flow<MinCostFlow> reduction" << std::endl;
    for (auto path: U_MPC(g, U)) {
        for (auto vertex: path) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;
    }


    std::cout << "U_MPC Min-Flow<Greedy+MaxFlow> reduction" << std::endl;
    for (auto path: greedy_U_MPC(g, U)) {
        for (auto vertex: path) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;
    }
*/
/*
    std::vector<ListDigraph::Node> S = {x0} ;
    std::vector<ListDigraph::Node> T = {x5, x4} ;

    std::cout << "Greedy Approximation MPC" << std::endl;
    for (auto path: greedy_approximation_MPC(g, S, T)) {
        for (auto vertex: path) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;
    }
*/
    /*
    std::cout << "Greedy Approximation U_MPC" << std::endl;
    for (auto path: greedy_approximation_U_MPC(g, {x0}, {x5, x4}, {x0, x2})) {
        for (auto vertex: path) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;
    }
*/
    /*
    ListDigraph::Node x0 = g.addNode();
    ListDigraph::Node x1 = g.addNode();
    ListDigraph::Node x2 = g.addNode();
    ListDigraph::Node x3 = g.addNode();
    ListDigraph::Node x4 = g.addNode();
    ListDigraph::Node x5 = g.addNode();

    g.addArc(x0, x2);
    g.addArc(x1, x3);
    g.addArc(x2, x3);
    g.addArc(x2, x4);
    g.addArc(x3, x5);

    std::vector<ListDigraph::Node> S = {x0, x1} ;
    std::vector<ListDigraph::Node> T = {x5, x4} ;
    std::vector<ListDigraph::Node> U = {x0, x2, x3, x4};

    std::cout << "MPC Min-Flow<MinCostFlow> reduction" << std::endl;
    for (auto path: MPC(g)) {
        for (auto vertex: path) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;
    }


    std::cout << "MPC Min-Flow<Greedy+MaxFlow> reduction" << std::endl;
    for (auto path: greedy_MPC(g)) {
        for (auto vertex: path) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;
    }

    std::cout << "(S,T) MPC Min-Flow<MinCostFlow> reduction" << std::endl;
    for (auto path: MPC(g, S, T)) {
        for (auto vertex: path) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;
    }


    std::cout << "(S,T) MPC Min-Flow<Greedy+MaxFlow> reduction" << std::endl;
    for (auto path: greedy_MPC(g, S, T)) {
        for (auto vertex: path) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;
    }

    std::cout << "U_MPC Min-Flow<MinCostFlow> reduction" << std::endl;
    for (auto path: U_MPC(g, U)) {
        for (auto vertex: path) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;
    }


    std::cout << "U_MPC Min-Flow<Greedy+MaxFlow> reduction" << std::endl;
    for (auto path: greedy_U_MPC(g, U)) {
        for (auto vertex: path) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;
    }

    std::cout << "(S,T) U_MPC Min-Flow<MinCostFlow> reduction" << std::endl;
    for (auto path: U_MPC(g, S, T, U)) {
        for (auto vertex: path) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;
    }

    std::cout << "(S,T) U_MPC Min-Flow<Greedy+MaxFlow> reduction" << std::endl;
    for (auto path: greedy_U_MPC(g, S, T, U)) {
        for (auto vertex: path) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;
    }
*/
/*
    std::cout << "Path Maximal safe paths MPC" << std::endl;
    for (auto pair: path_maximal_safe_paths_MPC(g)) {
        std::cout << "Path: ";
        for (auto vertex: pair.first) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;

        for (auto path: pair.second) {
            for (auto vertex: path) {
                std::cout << g.id(vertex) << " , ";
            }
            std::cout << std::endl;
        }
    }


    std::cout << "Greedy Path Maximal safe paths MPC" << std::endl;
    for (auto pair: greedy_path_maximal_safe_paths_MPC(g)) {
        std::cout << "Path: ";
        for (auto vertex: pair.first) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;

        for (auto path: pair.second) {
            for (auto vertex: path) {
                std::cout << g.id(vertex) << " , ";
            }
            std::cout << std::endl;
        }
    }
*/

/*
    std::vector<ListDigraph::Node> S = {x0};
    std::vector<ListDigraph::Node> T = {x7};
    std::cout << "Path Maximal safe paths MPC S, T" << std::endl;
    for (auto pair: path_maximal_safe_paths_MPC(g, S, T)) {
        std::cout << "Path: ";
        for (auto vertex: pair.first) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;

        for (auto path: pair.second) {
            for (auto vertex: path) {
                std::cout << g.id(vertex) << " , ";
            }
            std::cout << std::endl;
        }
    }



    std::cout << "Greedy Path Maximal safe paths MPC S, T" << std::endl;
    for (auto pair: greedy_path_maximal_safe_paths_MPC(g, S, T)) {
        std::cout << "Path: ";
        for (auto vertex: pair.first) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;

        for (auto path: pair.second) {
            for (auto vertex: path) {
                std::cout << g.id(vertex) << " , ";
            }
            std::cout << std::endl;
        }
    }
*/

/*
    int64_t l = 2;
    std::cout << "Path Maximal safe paths MPC" << std::endl;
    for (auto pair: path_maximal_safe_paths_PC(g, l)) {
        std::cout << "Path: ";
        for (auto vertex: pair.first) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;

        for (auto path: pair.second) {
            for (auto vertex: path) {
                std::cout << g.id(vertex) << " , ";
            }
            std::cout << std::endl;
        }
    }

    std::cout << "Path Maximal safe paths MPC" << std::endl;
    for (auto pair: greedy_path_maximal_safe_paths_PC(g, l)) {
        std::cout << "Path: ";
        for (auto vertex: pair.first) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;

        for (auto path: pair.second) {
            for (auto vertex: path) {
                std::cout << g.id(vertex) << " , ";
            }
            std::cout << std::endl;
        }
    }
*/
/*
    int64_t l = 2;

    std::vector<ListDigraph::Node> S = {x0};
    std::vector<ListDigraph::Node> T = {x7};
    std::cout << "Path Maximal<MinCostFlow> safe paths MPC S, T" << std::endl;
    for (auto pair: path_maximal_safe_paths_PC(g, S, T, l)) {
        std::cout << "Path: ";
        for (auto vertex: pair.first) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;

        for (auto path: pair.second) {
            for (auto vertex: path) {
                std::cout << g.id(vertex) << " , ";
            }
            std::cout << std::endl;
        }
    }

    std::cout << "Path Maximal<Greedy+MaxFlow> safe paths MPC S, T" << std::endl;
    for (auto pair: greedy_path_maximal_safe_paths_PC(g, S, T, l)) {
        std::cout << "Path: ";
        for (auto vertex: pair.first) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;

        for (auto path: pair.second) {
            for (auto vertex: path) {
                std::cout << g.id(vertex) << " , ";
            }
            std::cout << std::endl;
        }
    }
    */
/*
    std::vector<ListDigraph::Node> U = {x1, x2};
    std::vector<ListDigraph::Node> S = {x0};
    std::vector<ListDigraph::Node> T = {x7};
    */
/*
    std::cout << "Path Maximal safe paths MPC" << std::endl;
    for (auto pair: path_maximal_safe_paths_MPC(g)) {
        std::cout << "Path: ";
        for (auto vertex: pair.first) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;

        for (auto path: pair.second) {
            for (auto vertex: path) {
                std::cout << g.id(vertex) << " , ";
            }
            std::cout << std::endl;
        }
    }
*/
/*
    std::cout << "Path Maximal safe paths MPC" << std::endl;
    for (auto pair: path_maximal_safe_paths_U_MPC(g, U)) {
        std::cout << "Path: ";
        for (auto vertex: pair.first) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;

        for (auto path: pair.second) {
            for (auto vertex: path) {
                std::cout << g.id(vertex) << " , ";
            }
            std::cout << std::endl;
        }
    }

    std::cout << "Path Maximal safe paths MPC" << std::endl;
    for (auto pair: greedy_path_maximal_safe_paths_U_MPC(g, U)) {
        std::cout << "Path: ";
        for (auto vertex: pair.first) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;

        for (auto path: pair.second) {
            for (auto vertex: path) {
                std::cout << g.id(vertex) << " , ";
            }
            std::cout << std::endl;
        }
    }

*/
/*
    std::cout << "Path Maximal safe paths MPC" << std::endl;
    for (auto pair: path_maximal_safe_paths_U_MPC(g, S, T, U)) {
        std::cout << "Path: ";
        for (auto vertex: pair.first) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;

        for (auto path: pair.second) {
            for (auto vertex: path) {
                std::cout << g.id(vertex) << " , ";
            }
            std::cout << std::endl;
        }
    }

    std::cout << "Path Maximal safe paths MPC" << std::endl;
    for (auto pair: greedy_path_maximal_safe_paths_U_MPC(g, S, T, U)) {
        std::cout << "Path: ";
        for (auto vertex: pair.first) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;

        for (auto path: pair.second) {
            for (auto vertex: path) {
                std::cout << g.id(vertex) << " , ";
            }
            std::cout << std::endl;
        }
    }



    std::cout << "Path Maximal safe paths MPC" << std::endl;
    for (auto pair: greedy_path_maximal_safe_paths_U_PC(g, S, T, U, 1000)) {
        std::cout << "Path: ";
        for (auto vertex: pair.first) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;

        for (auto path: pair.second) {
            for (auto vertex: path) {
                std::cout << g.id(vertex) << " , ";
            }
            std::cout << std::endl;
        }
    }
*/
/*
    std::cout << "Path Maximal safe paths MPC" << std::endl;
    for (auto pair: greedy_path_maximal_safe_paths_U_MPC(g, S, T, U)) {
        std::cout << "Path: ";
        for (auto vertex: pair.first) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;

        for (auto path: pair.second) {
            for (auto vertex: path) {
                std::cout << g.id(vertex) << " , ";
            }
            std::cout << std::endl;
        }
    }
*/

/*
    std::cout << "Path Maximal safe paths PC" << std::endl;
    for (auto pair: path_maximal_safe_paths_U_PC(g, U, 3)) {
        std::cout << "Path: ";
        for (auto vertex: pair.first) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;

        for (auto path: pair.second) {
            for (auto vertex: path) {
                std::cout << g.id(vertex) << " , ";
            }
            std::cout << std::endl;
        }
    }

    std::cout << "Path Maximal safe paths PC" << std::endl;
    for (auto pair: greedy_path_maximal_safe_paths_U_PC(g, U, 3)) {
        std::cout << "Path: ";
        for (auto vertex: pair.first) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;

        for (auto path: pair.second) {
            for (auto vertex: path) {
                std::cout << g.id(vertex) << " , ";
            }
            std::cout << std::endl;
        }
    }

    std::cout << "Path Maximal safe paths MPC" << std::endl;
    for (auto pair: path_maximal_safe_paths_U_PC(g, S, T, U, 3)) {
        std::cout << "Path: ";
        for (auto vertex: pair.first) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;

        for (auto path: pair.second) {
            for (auto vertex: path) {
                std::cout << g.id(vertex) << " , ";
            }
            std::cout << std::endl;
        }
    }
*/


    return 0;
}
