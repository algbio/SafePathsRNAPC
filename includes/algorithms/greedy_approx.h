#ifndef SAFEPATHSRNAPC_GREEDY_APPROX_H
#define SAFEPATHSRNAPC_GREEDY_APPROX_H

#include <lemon/list_graph.h>


/*
 * Computes a path cover by using the greedy approximation
 * algorithm described in
 * "sparse dynamic programming on dags with small width"
 *
 * It assumes g is a DAG
 */
std::vector<std::vector<lemon::ListDigraph::Node>> greedy_approximation_MPC(lemon::ListDigraph& g);


/*
 * Computes a path cover by using the greedy approximation
 * algorithm described in
 * "sparse dynamic programming on dags with small width"
 *
 * It assumes g is a DAG
 * It return the path as the sequence of edges
 */
std::vector<std::vector<lemon::ListDigraph::Arc>> greedy_approximation_MPC_edges(lemon::ListDigraph& g);


/*
 * Computes a path cover by using the greedy approximation
 * algorithm described in
 * "sparse dynamic programming on dags with small width".
 * The paths of the path cover must start at some vertex
 * of S and end at some vertex of T
 *
 * It assumes g is a DAG, sources(g) \subseteq S, sinks(g) \subseteq T
 */
std::vector<std::vector<lemon::ListDigraph::Node>> greedy_approximation_MPC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& S, std::vector<lemon::ListDigraph::Node>& T);


/*
 * Computes a path cover by using the greedy approximation
 * algorithm described in
 * "sparse dynamic programming on dags with small width".
 * The paths of the path cover must start at some vertex
 * of S and end at some vertex of T
 *
 * It assumes g is a DAG, sources(g) \subseteq S, sinks(g) \subseteq T
 * It return the paths as a sequence of the edges
 */
std::vector<std::vector<lemon::ListDigraph::Arc>> greedy_approximation_MPC_edges(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& S, std::vector<lemon::ListDigraph::Node>& T);


/*
 * Computes a path cover, covering the vertices in U
 * by using the greedy approximation
 * algorithm described in
 * "sparse dynamic programming on dags with small width"
 *
 * It assumes g is a DAG, and U \subseteq V(g)
 */
std::vector<std::vector<lemon::ListDigraph::Node>> greedy_approximation_U_MPC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& U);


/*
 * Computes a path cover covering the vertices in U
 * by using the greedy approximation
 * algorithm described in
 * "sparse dynamic programming on dags with small width"
 *
 * It assumes g is a DAG, and U \subseteq V(g)
 * It return the path as the sequence of edges
 */
std::vector<std::vector<lemon::ListDigraph::Arc>> greedy_approximation_U_MPC_edges(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& U);


/*
 * Computes a path cover covering the vertices in U
 * by using the greedy approximation
 * algorithm described in
 * "sparse dynamic programming on dags with small width".
 * The paths of the path cover must start at some vertex
 * of S and end at some vertex of T
 *
 * It assumes g is a DAG, sources(g) \subseteq S, sinks(g) \subseteq T, and U \subseteq V(g)
 */
std::vector<std::vector<lemon::ListDigraph::Node>> greedy_approximation_U_MPC(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& S, std::vector<lemon::ListDigraph::Node>& T, std::vector<lemon::ListDigraph::Node>& U);


/*
 * Computes a path cover, covering the vertices in U
 * by using the greedy approximation
 * algorithm described in
 * "sparse dynamic programming on dags with small width".
 * The paths of the path cover must start at some vertex
 * of S and end at some vertex of T
 *
 * It assumes g is a DAG, sources(g) \subseteq S, sinks(g) \subseteq T, and U \subseteq V(g)
 * It return the paths as a sequence of the edges
 */
std::vector<std::vector<lemon::ListDigraph::Arc>> greedy_approximation_U_MPC_edges(lemon::ListDigraph& g, std::vector<lemon::ListDigraph::Node>& S, std::vector<lemon::ListDigraph::Node>& T, std::vector<lemon::ListDigraph::Node>& U);


#endif //SAFEPATHSRNAPC_GREEDY_APPROX_H
