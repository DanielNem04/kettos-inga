cmake_minimum_required(VERSION 3.10)
project(kettos-inga)

set(CMAKE_CXX_STANDARD 17)

add_executable(my_program main.cpp)
add_library(pendulum SHARED main2.cpp)

add_executable(my_plotter main3.cpp)
