set(project_EXECUTABLES
        run_mpc
        run_safe_paths
        run_safe_paths_two_finger
        run_safe_paths_safe_edges_opt
        run_safe_paths_unoptimized)


foreach( exe ${project_EXECUTABLES} )
    add_executable( ${exe} ${exe}.cpp )
endforeach()


foreach( exe ${project_EXECUTABLES} )
    target_link_libraries( ${exe} safepathsrnapc )
    target_link_libraries( ${exe} lemon )
endforeach()