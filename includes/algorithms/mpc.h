#ifndef SAFEPATHSRNAPC_MPC_H
#define SAFEPATHSRNAPC_MPC_H

#include <lemon/list_graph.h>


/*
 * Computes a minimum path cover by
 * reducing the problem to Min-flow<MinCostFlow>
 *
 * It assumes g is a DAG
 */
std::vector<std::vector<lemon::ListDigraph::Node>> MPC(lemon::ListDigraph& g);


/*
 * Computes a minimum path cover by
 * reducing the problem to Min-flow<Greedy+Max-Flow>
 *
 * It assumes g is a DAG
 */
std::vector<std::vector<lemon::ListDigraph::Node>> greedy_MPC(lemon::ListDigraph& g);


/*
 * Computes a minimum path cover by
 * reducing the problem to Min-flow<MinCostFlow>
 * The paths of the path cover must start at some vertex
 * of S and end at some vertex of T
 *
 * It assumes g is a DAG, sources(g) \subseteq S, sinks(g) \subseteq T
 */
std::vector<std::vector<lemon::ListDigraph::Node>> MPC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& S, std::vector<lemon::ListDigraph::Node>& T);


/*
 * Computes a minimum path cover by
 * reducing the problem to Min-flow<Greedy+Max-Flow>
 * The paths of the path cover must start at some vertex
 * of S and end at some vertex of T
 *
 * It assumes g is a DAG, sources(g) \subseteq S, sinks(g) \subseteq T
 */
std::vector<std::vector<lemon::ListDigraph::Node>> greedy_MPC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& S, std::vector<lemon::ListDigraph::Node>& T);


/*
 * Computes a minimum path cover covering the
 * vertices in U by
 * reducing the problem to Min-flow<MinCostFlow>
 *
 * It assumes g is a DAG, and U \subseteq V(g)
 */
std::vector<std::vector<lemon::ListDigraph::Node>> U_MPC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& U);


/*
 * Computes a minimum path cover covering the
 * vertices in U by
 * reducing the problem to Min-flow<Greedy+Max-Flow>
 *
 * It assumes g is a DAG
 */
std::vector<std::vector<lemon::ListDigraph::Node>> greedy_U_MPC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& U);


/*
 * Computes a minimum path cover covering the
 * vertices in U by
 * reducing the problem to Min-flow<MinCostFlow>
 * The paths of the path cover must start at some vertex
 * of S and end at some vertex of T
 *
 * It assumes g is a DAG, sources(g) \subseteq S, sinks(g) \subseteq T, U \subseteq V(g)
 */
std::vector<std::vector<lemon::ListDigraph::Node>> U_MPC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& S, std::vector<lemon::ListDigraph::Node>& T, std::vector<lemon::ListDigraph::Node>& U);


/*
 * Computes a minimum path cover covering the vertices of U
 * by reducing the problem to Min-flow<Greedy+Max-Flow>
 * The paths of the path cover must start at some vertex
 * of S and end at some vertex of T
 *
 * It assumes g is a DAG, sources(g) \subseteq S, sinks(g) \subseteq T, U \subseteq V(g)
 */
std::vector<std::vector<lemon::ListDigraph::Node>> greedy_U_MPC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& S, std::vector<lemon::ListDigraph::Node>& T, std::vector<lemon::ListDigraph::Node>& U);


/*
 * Computes a minimum path cover covering the vertices of U
 * by reducing the problem to Min-flow<Greedy+Max-Flow>
 * The paths of the path cover must start at some vertex
 * of S and end at some vertex of T
 *
 * It assumes g is a DAG, sources(g) \subseteq S, sinks(g) \subseteq T, U \subseteq V(g)
 */
std::vector<std::vector<lemon::ListDigraph::Node>> greedy_U_MPC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& S, std::vector<lemon::ListDigraph::Node>& T, std::vector<lemon::ListDigraph::Node>& U);


#endif //SAFEPATHSRNAPC_MPC_H
