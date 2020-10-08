#include <algorithms/filter_paths.h>

#include <qsufsort.hpp>
#include <unordered_map>

using namespace lemon;

std::vector<std::vector<ListDigraph::Node>> filter_contained_paths(ListDigraph& g, std::vector<std::pair<std::vector<ListDigraph::Node>,std::vector<std::vector<ListDigraph::Node>>>>& safe_paths_per_path) {
    if (safe_paths_per_path.size() == 0) return {};

    std::vector<std::vector<ListDigraph::Node>> safe_paths;
    for (auto& pair: safe_paths_per_path) {
        std::vector<std::vector<ListDigraph::Node>> &paths = pair.second;
        for (auto path: paths) {
            safe_paths.push_back(path);
        }
    }

    int n = 1;
    int s = 0;
    int max_id = -1;
    int min_id = countNodes(g);

    std::unordered_map<int,int> start_to_path;
    for (auto& path: safe_paths) {
        start_to_path[n] = s;
        for (ListDigraph::Node v: path) {
            int id = g.id(v) + 2;
            max_id = (id > max_id) ? id : max_id;
            min_id = (id < min_id) ? id : min_id;
            ++n;
        }
        ++n;
        ++s;
    }



    int* array = new int[n+1];
    array[n] = 0;
    int* copy_array = new int[n+1];
    int j = 0;
    array[j++] = min_id-1;

    for (auto& path: safe_paths) {
        for (ListDigraph::Node v: path) {
            int id = g.id(v) + 2;
            array[j++] = id;
        }
        array[j++] = min_id-1;
    }



    for (int i = 0; i < n+1; ++i) {
        copy_array[i] = array[i];
    }

    int* sa = new int[n+1];

    qsufsort qss;
    qss.suffixsort(copy_array, sa, n, max_id+1, min_id-1);

    std::vector<std::vector<ListDigraph::Node>> filtered_safe_paths;
    std::vector<bool> reported(s, false);
    for (int i = 0; i < safe_paths.size(); ++i) {
        std::vector<ListDigraph::Node>& path = safe_paths[i];
        if (!reported[i]) {
            bool report = true;
            std::vector<int> path_ids;
            for (ListDigraph::Node v: path) {
                path_ids.push_back(g.id(v)+2);
            }

            for (int j : qss.locate(path_ids, array, n, sa)) {
                if (array[j-1] == min_id-1 && array[j+path.size()] == min_id-1) {
                    reported[start_to_path[j]] = true;
                } else {
                    report = false;
                }
            }

            if (report) {
                filtered_safe_paths.push_back(path);
            }
        }
    }

    delete [] array;
    delete [] copy_array;
    delete [] sa;

    return filtered_safe_paths;
}