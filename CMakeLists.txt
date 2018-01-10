cmake_minimum_required(VERSION 3.6.3)
project(parallel_ex4)



set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall Werror")

set(SOURCE_FILES transpose.c)
add_executable(parallel_ex4 ${SOURCE_FILES})