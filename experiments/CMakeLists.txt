set(project_EXECUTABLES
        run_mpc
        run_safe_paths)


foreach( exe ${project_EXECUTABLES} )
    add_executable( ${exe} ${exe}.cpp )
endforeach()


foreach( exe ${project_EXECUTABLES} )
    target_link_libraries( ${exe} safepathsrnapc )
    target_link_libraries( ${exe} lemon )
endforeach()