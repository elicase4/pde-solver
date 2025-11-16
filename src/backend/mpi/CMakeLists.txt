find_package(MPI REQUIRED)

set(MPI_SOURCES
	apply_operator.cpp
	update_solution.cpp
	compute_residual.cpp
	decompose.cpp
	halo_exchange.cpp
)

add_library(mpi_backend ${MPI_SOURCES})
target_include_directories(mpi_backend PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(mpi_backend PUBLIC core MPI::MPI_CXX)
