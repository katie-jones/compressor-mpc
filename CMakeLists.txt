# Project definition
cmake_minimum_required(VERSION 2.6)
project(Compressor_MPC CXX)

# Use gcc 5.2 to get c++14 standard
set( CMAKE_CXX_COMPILER "/usr/local/bin/g++-5.2" )

# export library path to find correct versions of shit
set(ENV{LD_LIBRARY_PATH} "/usr/local/lib64" )

# export commands for ycm
set( CMAKE_EXPORT_COMPILE_COMMANDS 1 )

# use c++14 standard
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

# add DEBUG option when running under debug mode
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

# qpOASES library
set( QPOASES_ROOT /home/katie/school/EPFL/4-S2016/MasterThesis/qpOASES-3.2.0)
INCLUDE_DIRECTORIES(${QPOASES_ROOT}/include)
LINK_DIRECTORIES(${QPOASES_ROOT}/bin)

# Boost timer library
set(BOOST_ROOT "/mnt/shared/system/usr/local/include/boost_1_60_0")
find_package(Boost REQUIRED timer)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# source subdirectory
add_subdirectory(src)

# libraries subdirectory
add_subdirectory(libs)

# dynamic systems subdirectory
add_subdirectory(systems)

