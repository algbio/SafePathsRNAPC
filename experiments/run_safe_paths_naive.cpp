#include <iostream>
#include <sys/resource.h>
#include <lemon/list_graph.h>

#include <utils.h>
#include <algorithms/safe_paths.h>
#include <algorithms/filter_paths.h>


int main(int argc, char*argv[]) {
    std::cout << "Input_graph = " <<  argv[1] << std::endl;
    int64_t l = atoi(argv[2]);
    lemon::ListDigraph g;
    lemon::ListDigraph::NodeMap<int64_t> original_id(g);
    std::vector<lemon::ListDigraph::Node> S, T, U;
    load_problem_instance(argv[1], g, original_id, S, T, U);


    std::cout << "Safe Paths, l = " << l << std::endl;
    struct rusage usage;
    struct timeval ru_start, rs_start, ru_end, rs_end;
    getrusage(RUSAGE_SELF, &usage);
    rs_start = usage.ru_stime;
    ru_start = usage.ru_utime;
    std::vector<std::pair<std::vector<lemon::ListDigraph::Node>, std::vector<std::vector<lemon::ListDigraph::Node>>>> safe_paths_per_path = path_maximal_safe_paths_U_PC(g,S,T,U,l);
    getrusage(RUSAGE_SELF, &usage);
    rs_end = usage.ru_stime;
    ru_end = usage.ru_utime;
    long user_time = (ru_end.tv_sec - ru_start.tv_sec)*1000000 + (ru_end.tv_usec - ru_start.tv_usec);
    long system_time = (rs_end.tv_sec - rs_start.tv_sec)*1000000 + (rs_end.tv_usec - rs_start.tv_usec);

    long safe_paths_time = user_time+system_time;


    getrusage(RUSAGE_SELF, &usage);
    rs_start = usage.ru_stime;
    ru_start = usage.ru_utime;
    std::vector<std::vector<lemon::ListDigraph::Node>> filtered_safe_paths = filter_contained_paths(g, safe_paths_per_path);
    getrusage(RUSAGE_SELF, &usage);
    rs_end = usage.ru_stime;
    ru_end = usage.ru_utime;
    user_time = (ru_end.tv_sec - ru_start.tv_sec)*1000000 + (ru_end.tv_usec - ru_start.tv_usec);
    system_time = (rs_end.tv_sec - rs_start.tv_sec)*1000000 + (rs_end.tv_usec - rs_start.tv_usec);

    long filter_time = user_time+system_time;

/*
    std::cout<< "Number of safe paths = " << filtered_safe_paths.size() << std::endl;
    for (auto&path : filtered_safe_paths) {
        for (int i = 0; i < path.size(); ++i) {
            lemon::ListDigraph::Node u = path[i];
            std::cout << original_id[u];
            if (i != path.size()-1) {
                std::cout << ",";
            }
        }
        std::cout << std::endl;
    }
    */

    std::cout << "Time difference Safe Path (not filtered) (µs) = " << safe_paths_time << std::endl;
    std::cout << "Time difference Safe Path filter) (µs) = " << filter_time << std::endl;
    std::cout << std::endl;


}