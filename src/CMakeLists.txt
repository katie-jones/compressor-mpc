# Project definition
cmake_minimum_required(VERSION 2.6)
project(Compressor_MPC CXX)
set( CMAKE_EXPORT_COMPILE_COMMANDS 1 )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# qpOASES library
set( QPOASES_ROOT ../../qpOASES-3.2.0)
INCLUDE_DIRECTORIES(${QPOASES_ROOT}/include)
LINK_DIRECTORIES(${QPOASES_ROOT}/bin)

# Boost timer library
set(BOOST_ROOT "/mnt/shared/system/usr/local/include/boost_1_60_0")
find_package(Boost REQUIRED timer)

add_executable(comptest linearized_test.cc)
add_executable(timetest delay_test.cc)
add_library(compressor compressor.cc)
add_library(tank tank.cc)
add_library(parcomp parallel_compressors.cc)
add_library(mpc_controller mpc_controller.cc)
target_link_libraries(mpc_controller qpOASES)
target_link_libraries(comptest compressor tank parcomp mpc_controller)

if(Boost_FOUND)
   include_directories(${Boost_INCLUDE_DIRS})
   target_link_libraries(comptest ${Boost_LIBRARIES})
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
