add_executable(poisson_solver
	#poisson_solver.cpp
)

target_include_directories(poisson_solver
	PRIVATE
		${CMAKE_SOURCE_DIR}/include/apps
		${CMAKE_SOURCE_DIR}/src/apps
)

target_link_libraries(poisson_solver
	PRIVATE
		pde_core
		pde_mesh
		pde_fem
		pde_equations
		pde_bc
		pde_linalg
		pde_utils
)

if (ENABLE_OPENMP)
	target_link_libraries(poisson_solver PRIVATE OpenMP::OpenMP_CXX)
	target_compile_definitions(poisson_solver PRIVATE USE_OPENMP)
endif()

if (ENABLE_MPI)
	target_link_libraries(poisson_solver PRIVATE MPI::MPI_CXX)
	target_compile_definitions(poisson_solver PRIVATE USE_MPI)
endif()

if (ENABLE_CUDA)
	target_link_libraries(poisson_solver PRIVATE CUDA::cudart)
	target_compile_definitions(poisson_solver PRIVATE USE_CUDA)
endif()

install(TARGETS poisson_solver
	RUNTIME DESTINATION bin
)
