cmake_minimum_required(VERSION 2.8)
project(molecule)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} -O4 -DRUN_OMP -march=native -std=c++11 -fopenmp ")

set(SOURCE_FILES main.cpp linalg.h utils.h bbfmm.h mapping.h quadratue.h gmres.h molecule.h)

include_directories("/home/lurker/Workspace/Libs/openblas/include")

add_executable(molecule ${SOURCE_FILES})
link_directories("/home/lurker/Workspace/Libs/openblas/lib")
target_link_libraries(molecule m  /home/lurker/Workspace/Libs/openblas/lib/libopenblas.so)