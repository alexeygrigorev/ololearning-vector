cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++0x")

project(ololearning)


add_subdirectory(libololearning)

enable_testing()
add_subdirectory(test)