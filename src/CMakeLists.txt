# Project definition
cmake_minimum_required(VERSION 2.6)
project(Compressor_MPC)
set( CMAKE_EXPORT_COMPILE_COMMANDS 1 )


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(test test.cc)
add_library(compressor compressor.cc)
target_link_libraries(test compressor)

add_library(tank tank.cc)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

