file(GLOB SERIAL_SOURCES
	apply_operator.cpp
	update_solution.cpp
	compute_residual.cpp
)

add_library(serial_backend ${SERIAL_SOURCES})
target_include_directories(serial_backend PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(serial_backend PUBLIC core)
