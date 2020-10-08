cmake_minimum_required(VERSION 3.17)
project(SafePathsRNAPC)

set(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O9")


include_directories(${LEMON_SOURCE_DIR})
include_directories(${LEMON_BINARY_DIR})
include_directories(includes)
include_directories(external/qsufsort/include)


add_subdirectory(external)
add_subdirectory(src)
add_subdirectory(executables)
add_subdirectory(experiments)