#ifndef SAFEPATHSRNAPC_SAFE_PATHS_H
#define SAFEPATHSRNAPC_SAFE_PATHS_H

#include <lemon/list_graph.h>



/*
 * Computes safe paths present as
 * subpaths of every MPC of g.
 *
 * First computes a MPC with the
 * MinFlow<MinCostFlow> reduction, then for
 * every path in this path cover
 * runs the two finger algorithm
 * to find maximal safe path inside
 * every path.
 *
 * It returns a list of paths of a MPC
 * and its corresponding safe_paths
 */
std::vector<std::pair<std::vector<lemon::ListDigraph::Node>, std::vector<std::vector<lemon::ListDigraph::Node>>>> path_maximal_safe_paths_MPC(lemon::ListDigraph& g);



/*
 * Computes safe paths present as
 * subpaths of every MPC of g.
 *
 * First computes a MPC with the
 * MinFlow<Greedy+MaxFlow> reduction, then for
 * every path in this path cover
 * runs the two finger algorithm
 * to find maximal safe path inside
 * every path.
 *
 * It returns a list of paths of a MPC
 * and its corresponding safe_paths
 */
std::vector<std::pair<std::vector<lemon::ListDigraph::Node>, std::vector<std::vector<lemon::ListDigraph::Node>>>> greedy_path_maximal_safe_paths_MPC(lemon::ListDigraph& g);




/*
 * Computes safe paths present as
 * subpaths of every MPC of g,
 * with path starting at S and
 * ending at T.
 *
 * First computes a MPC with paths
 * starting at S and ending at T
 * with the MinFlow<MinCostFlow> reduction, then for
 * every path in this path cover
 * runs the two finger algorithm
 * to find maximal safe path inside
 * every path.
 *
 * It returns a list of paths of a MPC,
 * with paths starting at S and ending
 * at T and its corresponding safe_paths
 */
std::vector<std::pair<std::vector<lemon::ListDigraph::Node>, std::vector<std::vector<lemon::ListDigraph::Node>>>> path_maximal_safe_paths_MPC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& S, std::vector<lemon::ListDigraph::Node>& T);


/*
 * Computes safe paths present as
 * subpaths of every MPC of g,
 * with path starting at S and
 * ending at T.
 *
 * First computes a MPC with paths
 * starting at S and ending at T
 * with the MinFlow<Greedy+MaxFlow> reduction, then for
 * every path in this path cover
 * runs the two finger algorithm
 * to find maximal safe path inside
 * every path.
 *
 * It returns a list of paths of a MPC,
 * with paths starting at S and ending
 * at T and its corresponding safe_paths
 */
std::vector<std::pair<std::vector<lemon::ListDigraph::Node>, std::vector<std::vector<lemon::ListDigraph::Node>>>> greedy_path_maximal_safe_paths_MPC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& S, std::vector<lemon::ListDigraph::Node>& T);


/*
 * Computes safe paths present as
 * subpaths of every MPC of g covering the
 * vertices in U.
 *
 * First computes a MPC covering the
 * vertices in U with the
 * MinFlow<MinCostFlow> reduction, then for
 * every path in this path cover
 * runs the two finger algorithm
 * to find maximal safe path inside
 * every path.
 *
 * It returns a list of paths of a MPC
 * covering the vertices in U
 * and its corresponding safe_paths
 */
std::vector<std::pair<std::vector<lemon::ListDigraph::Node>, std::vector<std::vector<lemon::ListDigraph::Node>>>> path_maximal_safe_paths_U_MPC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& U);


/*
 * Computes safe paths present as
 * subpaths of every MPC of g covering the
 * vertices in U.
 *
 * First computes a MPC covering the
 * vertices in U with the
 * MinFlow<Greedy+MaxFlow> reduction, then for
 * every path in this path cover
 * runs the two finger algorithm
 * to find maximal safe path inside
 * every path.
 *
 * It returns a list of paths of a MPC
 * covering the vertices in U
 * and its corresponding safe_paths
 */
std::vector<std::pair<std::vector<lemon::ListDigraph::Node>, std::vector<std::vector<lemon::ListDigraph::Node>>>> greedy_path_maximal_safe_paths_U_MPC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& U);


/*
 * Computes safe paths present as
 * subpaths of every MPC of g, covering the
 * vertices in U,
 * with path starting at S and
 * ending at T.
 *
 * First computes a MPC covering the
 * vertices in U with paths
 * starting at S and ending at T
 * with the MinFlow<MinCostFlow> reduction, then for
 * every path in this path cover
 * runs the two finger algorithm
 * to find maximal safe path inside
 * every path.
 *
 * It returns a list of paths of a MPC covering the
 * vertices in U,
 * with paths starting at S and ending
 * at T and its corresponding safe_paths
 */
std::vector<std::pair<std::vector<lemon::ListDigraph::Node>, std::vector<std::vector<lemon::ListDigraph::Node>>>> path_maximal_safe_paths_U_MPC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& S, std::vector<lemon::ListDigraph::Node>& T, std::vector<lemon::ListDigraph::Node>& U);


/*
 * Computes safe paths present as
 * subpaths of every MPC of g, covering the
 * vertices in U,
 * with path starting at S and
 * ending at T.
 *
 * First computes a MPC covering the
 * vertices in U with paths
 * starting at S and ending at T
 * with the MinFlow<Greedy+MaxFlow> reduction, then for
 * every path in this path cover
 * runs the two finger algorithm
 * to find maximal safe path inside
 * every path.
 *
 * It returns a list of paths of a MPC covering the
 * vertices in U,
 * with paths starting at S and ending
 * at T and its corresponding safe_paths
 */
std::vector<std::pair<std::vector<lemon::ListDigraph::Node>, std::vector<std::vector<lemon::ListDigraph::Node>>>> greedy_path_maximal_safe_paths_U_MPC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& S, std::vector<lemon::ListDigraph::Node>& T, std::vector<lemon::ListDigraph::Node>& U);


/*
 * Computes safe paths present as
 * subpaths of every PC of size <= l of g.
 * (it assumes l \in [width(G) .... 2width(G)-1])
 *
 * First computes a MPC with the
 * MinFlow<MinCostFlow> reduction, then for
 * every path in this path cover
 * runs the two finger algorithm
 * to find maximal safe path inside
 * every path.
 *
 * It returns a list of paths of a MPC
 * and its corresponding safe_paths
 */
std::vector<std::pair<std::vector<lemon::ListDigraph::Node>, std::vector<std::vector<lemon::ListDigraph::Node>>>> path_maximal_safe_paths_PC(lemon::ListDigraph& g, int64_t l);


/*
 * Computes safe paths present as
 * subpaths of every PC of size <= l of g.
 * (it assumes l \in [width(G) .... 2width(G)-1])
 *
 * First computes a MPC with the
 * MinFlow<Greedy+MaxFlow> reduction, then for
 * every path in this path cover
 * runs the two finger algorithm
 * to find maximal safe path inside
 * every path.
 *
 * It returns a list of paths of a MPC
 * and its corresponding safe_paths
 */
std::vector<std::pair<std::vector<lemon::ListDigraph::Node>, std::vector<std::vector<lemon::ListDigraph::Node>>>> greedy_path_maximal_safe_paths_PC(lemon::ListDigraph& g, int64_t l);


/*
 * Computes safe paths present subpaths
 * of every PC of size <= l of g.
 * (it assumes l \in [width(G) .... 2width(G)-1]),
 * with path starting at S and
 * ending at T.
 *
 * First computes a MPC with paths
 * starting at S and ending at T
 * with the MinFlow<MinCostFlow> reduction, then for
 * every path in this path cover
 * runs the two finger algorithm
 * to find maximal safe path inside
 * every path.
 *
 * It returns a list of paths of a MPC,
 * with paths starting at S and ending
 * at T and its corresponding safe_paths
 */
std::vector<std::pair<std::vector<lemon::ListDigraph::Node>, std::vector<std::vector<lemon::ListDigraph::Node>>>> path_maximal_safe_paths_PC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& S, std::vector<lemon::ListDigraph::Node>& T, int64_t l);


/*
 * Computes safe paths present subpaths
 * of every PC of size <= l of g.
 * (it assumes l \in [width(G) .... 2width(G)-1]),
 * with path starting at S and
 * ending at T.
 *
 * First computes a MPC with paths
 * starting at S and ending at T
 * with the MinFlow<Greedy+MaxFlow> reduction, then for
 * every path in this path cover
 * runs the two finger algorithm
 * to find maximal safe path inside
 * every path.
 *
 * It returns a list of paths of a MPC,
 * with paths starting at S and ending
 * at T and its corresponding safe_paths
 */
std::vector<std::pair<std::vector<lemon::ListDigraph::Node>, std::vector<std::vector<lemon::ListDigraph::Node>>>> greedy_path_maximal_safe_paths_PC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& S, std::vector<lemon::ListDigraph::Node>& T, int64_t l);


/*
 * Computes safe paths present as
 * subpaths of every PC of size <= l of g.
 * (it assumes l \in [width(G) .... 2width(G)-1]), covering the
 * vertices in U.
 *
 * First computes a MPC covering the
 * vertices in U with the
 * MinFlow<MinCostFlow> reduction, then for
 * every path in this path cover
 * runs the two finger algorithm
 * to find maximal safe path inside
 * every path.
 *
 * It returns a list of paths of a MPC
 * covering the vertices in U
 * and its corresponding safe_paths
 */
std::vector<std::pair<std::vector<lemon::ListDigraph::Node>, std::vector<std::vector<lemon::ListDigraph::Node>>>> path_maximal_safe_paths_U_PC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& U, int64_t l);


/*
 * Computes safe paths present as
 * subpaths of every PC of size <= l of g.
 * (it assumes l \in [width(G) .... 2width(G)-1]), covering the
 * vertices in U.
 *
 * First computes a MPC covering the
 * vertices in U with the
 * MinFlow<Greedy+MaxFlow> reduction, then for
 * every path in this path cover
 * runs the two finger algorithm
 * to find maximal safe path inside
 * every path.
 *
 * It returns a list of paths of a MPC
 * covering the vertices in U
 * and its corresponding safe_paths
 */
std::vector<std::pair<std::vector<lemon::ListDigraph::Node>, std::vector<std::vector<lemon::ListDigraph::Node>>>> greedy_path_maximal_safe_paths_U_PC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& U, int64_t l);



/*
 * Computes safe paths present as
 * subpaths of every PC of size <= l of g.
 * (it assumes l \in [width(G) .... 2width(G)-1]), covering the
 * vertices in U,
 * with path starting at S and
 * ending at T.
 *
 * First computes a MPC covering the
 * vertices in U with paths
 * starting at S and ending at T
 * with the MinFlow<MinCostFlow> reduction, then for
 * every path in this path cover
 * finds all paths of length 1,
 * then the ones of length 2,
 * removing the contained of length 1,
 * then the ones of length 3 and so on.
 *
 * It returns a list of paths of a MPC covering the
 * vertices in U,
 * with paths starting at S and ending
 * at T and its corresponding safe_paths
 */
std::vector<std::pair<std::vector<lemon::ListDigraph::Node>, std::vector<std::vector<lemon::ListDigraph::Node>>>> naive_path_maximal_safe_paths_U_PC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& S, std::vector<lemon::ListDigraph::Node>& T, std::vector<lemon::ListDigraph::Node>& U, int64_t l);



/*
 * Computes safe paths present as
 * subpaths of every PC of size <= l of g.
 * (it assumes l \in [width(G) .... 2width(G)-1]), covering the
 * vertices in U,
 * with path starting at S and
 * ending at T.
 *
 * First computes a MPC covering the
 * vertices in U with paths
 * starting at S and ending at T
 * with the MinFlow<MinCostFlow> reduction, then for
 * every path in this path cover
 * runs the two finger algorithm
 * to find maximal safe path inside
 * every path.
 *
 * It returns a list of paths of a MPC covering the
 * vertices in U,
 * with paths starting at S and ending
 * at T and its corresponding safe_paths
 */
std::vector<std::pair<std::vector<lemon::ListDigraph::Node>, std::vector<std::vector<lemon::ListDigraph::Node>>>> path_maximal_safe_paths_U_PC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& S, std::vector<lemon::ListDigraph::Node>& T, std::vector<lemon::ListDigraph::Node>& U, int64_t l);


/*
 * Computes safe paths present as
 * subpaths of every PC of size <= l of g.
 * (it assumes l \in [width(G) .... 2width(G)-1]) of g, covering the
 * vertices in U,
 * with path starting at S and
 * ending at T.
 *
 * First computes a MPC covering the
 * vertices in U with paths
 * starting at S and ending at T
 * with the MinFlow<Greedy+MaxFlow> reduction, then for
 * every path in this path cover
 * runs the two finger algorithm
 * to find maximal safe path inside
 * every path.
 *
 * It returns a list of paths of a MPC covering the
 * vertices in U,
 * with paths starting at S and ending
 * at T and its corresponding safe_paths
 */
std::vector<std::pair<std::vector<lemon::ListDigraph::Node>, std::vector<std::vector<lemon::ListDigraph::Node>>>> greedy_path_maximal_safe_paths_U_PC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& S, std::vector<lemon::ListDigraph::Node>& T, std::vector<lemon::ListDigraph::Node>& U, int64_t l);



/*
 * Computes safe paths present as
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
 *
 * then it computes safe edges ,
 *
 * and then for
 * every path in the path cover
 * runs the two finger algorithm
 * to find maximal safe path inside
 * every path. Safe edges are used to skip
 * unnecessary computation in the algorithm
 *
 * It returns a list of paths of a MPC covering the
 * vertices in U,
 * with paths starting at S and ending
 * at T and its corresponding safe_paths
 */
std::vector<std::pair<std::vector<lemon::ListDigraph::Node>, std::vector<std::vector<lemon::ListDigraph::Node>>>> optimized_greedy_path_maximal_safe_paths_U_PC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& S, std::vector<lemon::ListDigraph::Node>& T, std::vector<lemon::ListDigraph::Node>& U, int64_t l);



#endif //SAFEPATHSRNAPC_SAFE_PATHS_H
