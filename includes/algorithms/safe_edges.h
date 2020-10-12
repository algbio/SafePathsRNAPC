#ifndef SAFEPATHSRNAPC_SAFE_EDGES_H
#define SAFEPATHSRNAPC_SAFE_EDGES_H

#include <lemon/list_graph.h>



/*
 * Computes safe edges present as
 * subpaths of every PC of size <= l of g.
 * (it assumes l \in [width(G) .... 2width(G)-1]) of g, covering the
 * vertices in U,
 * with path starting at S and
 * ending at T.
 *
 * First computes a MPC covering the
 * vertices in U with paths
 * starting at S and ending at T
 * with the MinFlow<Greedy+MaxFlow> reduction,
 * then for every edge e it tests if it is safe
 * by obtaining an MPC og G^e and comparing its size
 * against l
 *
 *
 * It returns a list of edges of g, the safe edges
 */
std::vector<lemon::ListDigraph::Arc> greedy_safe_edges_U_PC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& S, std::vector<lemon::ListDigraph::Node>& T, std::vector<lemon::ListDigraph::Node>& U, int64_t l);

#endif //SAFEPATHSRNAPC_SAFE_EDGES_H
