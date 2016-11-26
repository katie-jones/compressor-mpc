INCLUDE_DIRECTORIES(${QPOASES_ROOT}/include)
LINK_DIRECTORIES(${QPOASES_ROOT}/bin)

add_executable(centralized-parallel-timing centralized-parallel-with-timing.cc)
add_executable(centralized-serial-timing centralized-serial-with-timing.cc)

add_executable(cooperative-parallel-timing cooperative-parallel-with-timing.cc)
add_executable(cooperative-serial-timing cooperative-serial-with-timing.cc)

add_executable(noncoop-parallel-timing noncoop-parallel-with-timing.cc)
add_executable(noncoop-serial-timing noncoop-serial-with-timing.cc)

target_link_libraries(centralized-parallel-timing parallel_compressors serial_compressors distributed_controller auglinsys observer qpsolver)
target_link_libraries(centralized-serial-timing parallel_compressors serial_compressors distributed_controller auglinsys observer qpsolver)
target_link_libraries(cooperative-parallel-timing parallel_compressors serial_compressors distributed_controller auglinsys observer qpsolver)
target_link_libraries(cooperative-serial-timing parallel_compressors serial_compressors distributed_controller auglinsys observer qpsolver)
target_link_libraries(noncoop-parallel-timing parallel_compressors serial_compressors distributed_controller auglinsys observer qpsolver)
target_link_libraries(noncoop-serial-timing parallel_compressors serial_compressors distributed_controller auglinsys observer qpsolver)

if(Boost_FOUND)
   include_directories(${Boost_INCLUDE_DIRS})
   target_link_libraries(centralized-parallel-timing ${Boost_LIBRARIES})
   target_link_libraries(centralized-serial-timing ${Boost_LIBRARIES})
   target_link_libraries(cooperative-parallel-timing ${Boost_LIBRARIES})
   target_link_libraries(cooperative-serial-timing ${Boost_LIBRARIES})
   target_link_libraries(noncoop-parallel-timing ${Boost_LIBRARIES})
   target_link_libraries(noncoop-serial-timing ${Boost_LIBRARIES})
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
