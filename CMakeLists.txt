cmake_minimum_required(VERSION 2.8)
project(molecule)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11")

set(SOURCE_FILES main.cpp linalg.h utils.h bbfmm.h mapping.h quadratue.h gmres.h molecule.h)

include_directories("/usr/include/openblas")

add_executable(molecule ${SOURCE_FILES})
link_directories("/usr/lib/openblas-base")
target_link_libraries(molecule m blas)