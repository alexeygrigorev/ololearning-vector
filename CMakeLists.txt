cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "-Wall -std=gnu++0x" CACHE STRING "compile flags" FORCE)

project(ololearning)

add_subdirectory(libololearning)

enable_testing()
add_subdirectory(test)