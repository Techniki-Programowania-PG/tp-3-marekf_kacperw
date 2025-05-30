cmake_minimum_required(VERSION 3.14)
project(MyProject)

set(CMAKE_CXX_STANDARD 17)


add_subdirectory(extern/pybind11)
add_subdirectory(extern/matplotplusplus)
pybind11_add_module(signals main.cpp)



