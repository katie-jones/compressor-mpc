INCLUDE_DIRECTORIES(${QPOASES_ROOT}/include)
LINK_DIRECTORIES(${QPOASES_ROOT}/bin)

# add_executable(centralized_parallel centralized_parallel.cc)
# add_executable(coop_parallel coop_parallel.cc)
# add_executable(noncoop_parallel noncoop_parallel.cc)

# add_executable(centralized_serial centralized_serial.cc)
# add_executable(coop_serial coop_serial.cc)
add_executable(noncoop_serial noncoop_serial.cc)

# target_link_libraries(coop_parallel parallel_compressors distributed_controller auglinsys observer qpsolver)
# target_link_libraries(noncoop_parallel parallel_compressors distributed_controller auglinsys observer qpsolver)
# target_link_libraries(centralized_parallel parallel_compressors distributed_controller auglinsys observer qpsolver)

# target_link_libraries(centralized_serial serial_compressors distributed_controller auglinsys observer qpsolver)
# target_link_libraries(coop_serial serial_compressors distributed_controller auglinsys observer qpsolver)
target_link_libraries(noncoop_serial serial_compressors distributed_controller auglinsys observer qpsolver)

if(Boost_FOUND)
   include_directories(${Boost_INCLUDE_DIRS})

   # target_link_libraries(coop_parallel ${Boost_LIBRARIES})
   # target_link_libraries(noncoop_parallel ${Boost_LIBRARIES})
   # target_link_libraries(centralized_parallel ${Boost_LIBRARIES})

   # target_link_libraries(centralized_serial ${Boost_LIBRARIES})
   # target_link_libraries(coop_serial ${Boost_LIBRARIES})
   target_link_libraries(noncoop_serial ${Boost_LIBRARIES})

endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
