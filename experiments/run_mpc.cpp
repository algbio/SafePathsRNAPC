#include <iostream>
#include <sys/resource.h>
#include <lemon/list_graph.h>

#include <utils.h>
#include <algorithms/mpc.h>


int main(int argc, char*argv[]) {
    std::cout << "Input_graph = " <<  argv[1] << std::endl;
    lemon::ListDigraph g;
    lemon::ListDigraph::NodeMap<int64_t> original_id(g);
    std::vector<lemon::ListDigraph::Node> S, T, U;
    load_problem_instance(argv[1], g, original_id, S, T, U);


    struct rusage usage;
    struct timeval ru_start, rs_start, ru_end, rs_end;
    getrusage(RUSAGE_SELF, &usage);
    rs_start = usage.ru_stime;
    ru_start = usage.ru_utime;
    std::vector<std::vector<lemon::ListDigraph::Node>> minimum_path_cover = greedy_U_MPC(g,S,T,U);
    getrusage(RUSAGE_SELF, &usage);
    rs_end = usage.ru_stime;
    ru_end = usage.ru_utime;
    long user_time = (ru_end.tv_sec - ru_start.tv_sec)*1000000 + (ru_end.tv_usec - ru_start.tv_usec);
    long system_time = (rs_end.tv_sec - rs_start.tv_sec)*1000000 + (rs_end.tv_usec - rs_start.tv_usec);


    long mpc_time = user_time+system_time;
    int64_t width = minimum_path_cover.size();


    std::cout << "width = " << width << std::endl;
    std::cout << "Time difference Minimum Path Cover(µs) = " << mpc_time << std::endl;
    std::cout << std::endl;
}