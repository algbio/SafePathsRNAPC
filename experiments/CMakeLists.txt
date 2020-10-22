set(project_EXECUTABLES
        run_mpc
        run_safe_paths
        run_safe_paths_naive
        run_safe_paths_opt)


foreach( exe ${project_EXECUTABLES} )
    add_executable( ${exe} ${exe}.cpp )
endforeach()


foreach( exe ${project_EXECUTABLES} )
    target_link_libraries( ${exe} safepathsrnapc )
    target_link_libraries( ${exe} lemon )
endforeach()