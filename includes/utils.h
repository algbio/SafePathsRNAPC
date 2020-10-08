#ifndef SAFEPATHSRNAPC_UTILS_H
#define SAFEPATHSRNAPC_UTILS_H

#include <lemon/list_graph.h>


void load_problem_instance(char* filename, lemon::ListDigraph& g, lemon::ListDigraph::NodeMap<int64_t>& original_id, std::vector<lemon::ListDigraph::Node>& S, std::vector<lemon::ListDigraph::Node>& T, std::vector<lemon::ListDigraph::Node>& U);

#endif //SAFEPATHSRNAPC_UTILS_H
