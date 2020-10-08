#include <utils.h>

#include <lemon/lgf_reader.h>

using namespace lemon;

void load_problem_instance(char* filename, ListDigraph& g, ListDigraph::NodeMap<int64_t>& original_id, std::vector<ListDigraph::Node>& S, std::vector<ListDigraph::Node>& T, std::vector<ListDigraph::Node>& U) {
    ListDigraph::NodeMap<bool> in_S(g);
    ListDigraph::NodeMap<bool> in_T(g);
    ListDigraph::NodeMap<bool> in_U(g);
    digraphReader(g, filename)
            .nodeMap("original_id", original_id)
            .nodeMap("is_source", in_S)
            .nodeMap("is_target", in_T)
            .nodeMap("is_vertex_constrain", in_U)
            .run();

    for (ListDigraph::NodeIt v(g) ; v != INVALID; ++v) {
        if (in_S[v]){
            S.push_back(v);
        }
        if (in_T[v]){
            T.push_back(v);
        }
        if (in_U[v]){
            U.push_back(v);
        }
    }
}