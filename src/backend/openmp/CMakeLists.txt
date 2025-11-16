find_package(OpenMP REQUIRED)

set(OMP_SOURCES
	apply_operator.cpp
	update_solution.cpp
	compute_residual.cpp
)

add_library(openmp_backend ${OMP_SOURCES})
target_include_directories(openmp_backend PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(openmp_backend PUBLIC core OpenMP::OpenMP_CXX)
