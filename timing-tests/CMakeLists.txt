INCLUDE_DIRECTORIES(${QPOASES_ROOT}/include)
LINK_DIRECTORIES(${QPOASES_ROOT}/bin)

# add_executable(timing_centralized_parallel centralized_parallel.cc)
# add_executable(timing_coop_parallel coop_parallel.cc)
# add_executable(timing_noncoop_parallel noncoop_parallel.cc)

add_executable(timing_centralized_serial centralized_serial.cc)
add_executable(timing_coop_serial coop_serial.cc)
add_executable(timing_noncoop_serial noncoop_serial.cc)

# target_link_libraries(timing_coop_parallel parallel_compressors distributed_controller auglinsys observer qpsolver)
# target_link_libraries(timing_noncoop_parallel parallel_compressors distributed_controller auglinsys observer qpsolver)
# target_link_libraries(timing_centralized_parallel parallel_compressors distributed_controller auglinsys observer qpsolver)

target_link_libraries(timing_centralized_serial serial_compressors distributed_controller auglinsys observer qpsolver)
target_link_libraries(timing_coop_serial serial_compressors distributed_controller auglinsys observer qpsolver)
target_link_libraries(timing_noncoop_serial serial_compressors distributed_controller auglinsys observer qpsolver)

if(Boost_FOUND)
   include_directories(${Boost_INCLUDE_DIRS})

   # target_link_libraries(timing_coop_parallel ${Boost_LIBRARIES})
   # target_link_libraries(timing_noncoop_parallel ${Boost_LIBRARIES})
   # target_link_libraries(timing_centralized_parallel ${Boost_LIBRARIES})

   target_link_libraries(timing_centralized_serial ${Boost_LIBRARIES})
   target_link_libraries(timing_coop_serial ${Boost_LIBRARIES})
   target_link_libraries(timing_noncoop_serial ${Boost_LIBRARIES})

endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
