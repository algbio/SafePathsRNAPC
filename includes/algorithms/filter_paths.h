
#ifndef SAFEPATHSRNAPC_FILTER_PATHS_H
#define SAFEPATHSRNAPC_FILTER_PATHS_H

#include <lemon/list_graph.h>

std::vector<std::vector<lemon::ListDigraph::Node>> filter_contained_paths(lemon::ListDigraph& g, std::vector<std::pair<std::vector<lemon::ListDigraph::Node>,std::vector<std::vector<lemon::ListDigraph::Node>>>>& safe_paths_per_path);

#endif //SAFEPATHSRNAPC_FILTER_PATHS_H
