# Project definition
cmake_minimum_required(VERSION 2.6)
project(Compressor_MPC CXX)
set( CMAKE_EXPORT_COMPILE_COMMANDS 1 )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set( QPOASES_ROOT ../../katie-thesis/qpoases )
INCLUDE_DIRECTORIES(${QPOASES_ROOT}/include)

# compile qpOASES libraries
FILE(GLOB SRC ${QPOASES_ROOT}/src/*.cpp)
FILE(GLOB EXTRAS_SRC ${QPOASES_ROOT}/src/extras/*.cpp)

ADD_LIBRARY(qpOASES ${SRC})
ADD_LIBRARY(qpOASESextras ${SRC} ${EXTRAS_SRC})



add_executable(comptest linearized_test.cc)
add_library(compressor compressor.cc)
add_library(tank tank.cc)
add_library(parcomp parallel_compressors.cc)
add_library(mpc_controller mpc_controller.cc)
target_link_libraries(mpc_controller qpOASES)
target_link_libraries(comptest compressor tank parcomp mpc_controller)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

