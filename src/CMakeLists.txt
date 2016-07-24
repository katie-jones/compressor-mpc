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
set( QPOASES_ROOT ../../qpOASES-3.2.0)
INCLUDE_DIRECTORIES(${QPOASES_ROOT}/include)
LINK_DIRECTORIES(${QPOASES_ROOT}/bin)

# Boost timer library
set(BOOST_ROOT "/mnt/shared/system/usr/local/include/boost_1_60_0")
find_package(Boost REQUIRED timer)

add_executable(centralized_parallel centralized_parallel.cc)
add_executable(centralized_parallel_setuponly centralized_parallel_setuponly.cc)
add_executable(coop_parallel coop_parallel.cc)
add_executable(noncoop_parallel noncoop_parallel.cc)

add_executable(centralized_serial centralized_serial.cc)
add_executable(centralized_serial_setuponly centralized_serial_setuponly.cc)
add_executable(coop_serial coop_serial.cc)
add_executable(noncoop_serial noncoop_serial.cc)

add_library(compressor compressor.cc)
add_library(tank tank.cc)
add_library(parallel_compressors parallel_compressors.cc)
add_library(serial_compressors serial_compressors.cc)
add_library(auglinsys aug_lin_sys.cc)
add_library(observer observer.cc)
add_library(qpsolver mpc_qp_solver.cc)
add_library(distributed_controller distributed_controller.cc)

target_link_libraries(parallel_compressors compressor tank)
target_link_libraries(serial_compressors compressor)
target_link_libraries(qpsolver qpOASES)

target_link_libraries(auglinsys serial_compressors parallel_compressors)

target_link_libraries(coop_parallel parallel_compressors distributed_controller auglinsys observer qpsolver)
target_link_libraries(noncoop_parallel parallel_compressors distributed_controller auglinsys observer qpsolver)
target_link_libraries(centralized_parallel parallel_compressors distributed_controller auglinsys observer qpsolver)
target_link_libraries(centralized_parallel_setuponly parallel_compressors distributed_controller auglinsys observer qpsolver)

target_link_libraries(centralized_serial serial_compressors distributed_controller auglinsys observer qpsolver)
target_link_libraries(centralized_serial_setuponly serial_compressors distributed_controller auglinsys observer qpsolver)
target_link_libraries(coop_serial serial_compressors distributed_controller auglinsys observer qpsolver)
target_link_libraries(noncoop_serial serial_compressors distributed_controller auglinsys observer qpsolver)

if(Boost_FOUND)
   include_directories(${Boost_INCLUDE_DIRS})

   target_link_libraries(coop_parallel ${Boost_LIBRARIES})
   target_link_libraries(noncoop_parallel ${Boost_LIBRARIES})
   target_link_libraries(centralized_parallel ${Boost_LIBRARIES})
   target_link_libraries(centralized_parallel_setuponly ${Boost_LIBRARIES})

   target_link_libraries(centralized_serial ${Boost_LIBRARIES})
   target_link_libraries(centralized_serial_setuponly ${Boost_LIBRARIES})
   target_link_libraries(coop_serial ${Boost_LIBRARIES})
   target_link_libraries(noncoop_serial ${Boost_LIBRARIES})

endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
