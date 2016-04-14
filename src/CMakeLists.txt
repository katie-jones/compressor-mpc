# Project definition
cmake_minimum_required(VERSION 2.6)
project(Compressor_MPC)
set( CMAKE_EXPORT_COMPILE_COMMANDS 1 )

# ------------------------------------------------------------ #
# ------------------------ GOOGLE TEST ----------------------- #
# ------------------------------------------------------------ #

# local gtest root
set(GTEST_DIR /usr/local/include/googletest/googletest)
set(GTEST_INCLUDE_DIRECTORIES ${GTEST_DIR}/include ${GTEST_DIR})
set(GTEST_BOTH_LIBRARIES gtest gtest_main)

include_directories(${GTEST_INCLUDE_DIRECTORIES})
add_subdirectory(${GTEST_DIR} gtest)



set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(test test.cc)
target_link_libraries(test ${GTEST_BOTH_LIBRARIES})

add_library(compressor compressor.cc)

add_executable(sim_comp sim_comp.cc)
target_link_libraries(sim_comp compressor)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

