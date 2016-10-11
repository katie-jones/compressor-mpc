INCLUDE_DIRECTORIES(${QPOASES_ROOT}/include)
LINK_DIRECTORIES(${QPOASES_ROOT}/bin)

add_executable(test_with_timing test-with-timing.cc)

target_link_libraries(test_with_timing parallel_compressors serial_compressors 
  distributed_controller auglinsys observer qpsolver)

if(Boost_FOUND)
   include_directories(${Boost_INCLUDE_DIRS})
   target_link_libraries(test_with_timing ${Boost_LIBRARIES})
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
