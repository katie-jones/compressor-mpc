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

# add_executable(comptest linearized_test.cc)
add_executable(parmpc parallel_mpc.cc)
add_executable(reduced_system_test reduced_system_test.cc)
add_executable(indextest index_test.cc)
# add_executable(cooperative_test cooperative_test.cc)
# add_executable(noncooperative_test noncooperative_test.cc)

add_library(compressor compressor.cc)
add_library(tank tank.cc)
add_library(parcomp parallel_compressors.cc)
add_library(mpc_controller mpc_controller.cc)
add_library(auglinsys aug_lin_sys.cc)
add_library(observer observer.cc)
add_library(qpsolver mpc_qp_solver.cc)
# add_library(cooperative_controller cooperative_controller.cc)
add_library(noncooperative_controller noncooperative_controller.cc)
add_library(distributed_controller distributed_controller.cc)

target_link_libraries(parcomp compressor tank)
target_link_libraries(qpsolver qpOASES)
target_link_libraries(mpc_controller auglinsys observer qpsolver)
# target_link_libraries(cooperative_controller auglinsys observer qpsolver)
target_link_libraries(noncooperative_controller auglinsys observer qpsolver)
# target_link_libraries(comptest compressor mpc_controller)
target_link_libraries(parmpc parcomp mpc_controller auglinsys observer qpsolver)
target_link_libraries(reduced_system_test parcomp)
# target_link_libraries(cooperative_test parcomp noncooperative_controller auglinsys observer qpsolver)
# target_link_libraries(noncooperative_test parcomp noncooperative_controller auglinsys observer qpsolver)

if(Boost_FOUND)
   include_directories(${Boost_INCLUDE_DIRS})
   # target_link_libraries(cooperative_test ${Boost_LIBRARIES})
   # target_link_libraries(noncooperative_test ${Boost_LIBRARIES})
   # target_link_libraries(cooperative_controller ${Boost_LIBRARIES})
   # target_link_libraries(noncooperative_controller ${Boost_LIBRARIES})
   target_link_libraries(parmpc ${Boost_LIBRARIES})
   # target_link_libraries(comptest ${Boost_LIBRARIES})
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
