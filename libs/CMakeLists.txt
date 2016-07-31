INCLUDE_DIRECTORIES(${QPOASES_ROOT}/include)
LINK_DIRECTORIES(${QPOASES_ROOT}/bin)


add_library(auglinsys aug_lin_sys.cc)
add_library(observer observer.cc)
add_library(qpsolver mpc_qp_solver.cc)
add_library(distributed_controller distributed_controller.cc)

target_link_libraries(qpsolver qpOASES)

target_link_libraries(auglinsys serial_compressors parallel_compressors)
