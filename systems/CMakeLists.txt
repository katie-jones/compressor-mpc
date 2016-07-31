add_library(compressor compressor.cc)
add_library(tank tank.cc)
add_library(parallel_compressors parallel_compressors.cc)
add_library(serial_compressors serial_compressors.cc)

target_link_libraries(parallel_compressors compressor tank)
target_link_libraries(serial_compressors compressor)
