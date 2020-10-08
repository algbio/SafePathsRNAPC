#ifndef SAFEPATHSRNAPC_TOP_SORT_H
#define SAFEPATHSRNAPC_TOP_SORT_H


#include <lemon/list_graph.h>
#include <lemon/dfs.h>


/*
 * It returns the topological sort obtained
 * by applying the kahn's algorithm in the
 * graph g.
 *
 * Method assumes that g is a DAG
 */
std::vector<lemon::ListDigraph::Node> topological_sort(lemon::ListDigraph& g);


/*
 * It returns the topological sort obtained
 * by applying the kahn's algorithm in the
 * graph g starting at the sources in S.
 *
 * Method assumes that g is a DAG and S \subseteq V(g)
 */
std::vector<lemon::ListDigraph::Node> topological_sort(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& S);



/*
 * It returns the topological sort obtained
 * by reverse post-order of dfs in the graph g
 *
 * Method assumes that g is a DAG
 */
std::vector<lemon::ListDigraph::Node> topological_sort_dfs(lemon::ListDigraph& g);


/*
 * It returns the topological sort obtained
 * by reverse post-order of dfs in the graph
 * g starting at the sources in S.
 *
 * Method assumes that g is a DAG and S \subseteq V(g)
 */
std::vector<lemon::ListDigraph::Node> topological_sort_dfs(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& S);


#endif //SAFEPATHSRNAPC_TOP_SORT_H
