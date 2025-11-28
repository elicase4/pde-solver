message(STATUS "Configuring tests with GoogleTest...")

# ==================================================================
# Functions
# ==================================================================

# Add unit tests
function(add_pde_test TEST_NAME SOURCE_FILE)
	add_executable(${TEST_NAME} ${SOURCE_FILE})

	target_include_directories(${TEST_NAME}
		PRIVATE
			${CMAKE_CURRENT_SOURCE_DIR}
			${CMAKE_SOURCE_DIR}/include
			${CMAKE_SOURCE_DIR}/src
	)
	
	# Link against GTest libraries
	target_link_libraries(${TEST_NAME}
		PRIVATE
			${ARGN}
			GTest::gtest_main
	)

	# Discover unit tests automatically
	gtest_discover_tests(${TEST_NAME}
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
		PROPERTIES
			LABELS "unit"
			TIMEOUT 30
	)

endfunction()

# Add integration tests
function(add_pde_integration_test TEST_NAME SOURCE_FILE)
	add_executable(${TEST_NAME} ${SOURCE_FILE})

	target_include_directories(${TEST_NAME}
		PRIVATE
			${CMAKE_CURRENT_SOURCE_DIR}
			${CMAKE_SOURCE_DIR}/include
			${CMAKE_SOURCE_DIR}/src
	)
	
	# Discover unit tests automatically
	gtest_discover_tests(${TEST_NAME}
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
		PROPERTIES
			LABELS "integration"
			TIMEOUT 120
	)

endfunction()

# ==================================================================
# Unit Tests
# ==================================================================

# Field tests
#add_pde_test(test_field
#	unit/test_field.cpp
#	pde_core pde_mesh
#)

# ==================================================================
# Backend Tests
# ==================================================================

# OpenMP
if(ENABLE_OPENMP)
    add_pde_test(test_openmp_backend
        unit/test_openmp_backend.cpp
        pde_linalg
    )
    target_link_libraries(test_openmp_backend PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "  Added OpenMP backend tests")
endif()

# MPI
if(ENABLE_MPI)
    add_executable(test_mpi_backend unit/test_mpi_backend.cpp)
    
    target_include_directories(test_mpi_backend
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${CMAKE_SOURCE_DIR}/include
            ${CMAKE_SOURCE_DIR}/src
    )
    
    target_link_libraries(test_mpi_backend
        PRIVATE
            pde_linalg
            MPI::MPI_CXX
            GTest::gtest_main
    )
    
    add_test(
        NAME test_mpi_backend
        COMMAND ${MPIEXEC_EXECUTABLE} 
                ${MPIEXEC_NUMPROC_FLAG} 4 
                ${MPIEXEC_PREFLAGS} 
                $<TARGET_FILE:test_mpi_backend> 
                ${MPIEXEC_POSTFLAGS}
    )
    
    set_tests_properties(test_mpi_backend PROPERTIES
        LABELS "mpi"
        TIMEOUT 60
    )
    message(STATUS "  Added MPI backend tests")
endif()

# CUDA
if(ENABLE_CUDA)
    add_pde_test(test_cuda_backend
        unit/test_cuda_backend.cpp
        pde_linalg
    )
    target_link_libraries(test_cuda_backend PRIVATE CUDA::cudart)
    set_target_properties(test_cuda_backend PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
    )
    message(STATUS "  Added CUDA backend tests")
endif()

# ==================================================================
# Integration Tests
# ==================================================================

# Poisson Eq test
#add_pde_integration_test(test_poisson_eqs
#	integration/test_poisson_eqs.cpp
#	pde_core pde_mesh pde_fem pde_equations pde_bc pde_linalg pde_utils
#)

# ==================================================================
# Integration and Backend Tests
# ==================================================================

if(ENABLE_OPENMP OR ENABLE_MPI OR ENABLE_CUDA)

	add_pde_integration_test(test_backend_comparison
		integration/test_backend_comparison.cpp
		pde_core pde_mesh pde_fem pde_equations pde_bc pde_linalg pde_utils
	)

	if(ENABLE_OPENMP)
		target_link_libraries(test_backend_comparison PRIVATE OpenMP::OpenMP_CXX)
	endif()

	if(ENABLE_MPI)
		target_link_libraries(test_backend_comparison PRIVATE MPI::MPI_CXX)
	endif()

	if(ENABLE_CUDA)
		target_link_libraries(test_backend_comparison PRIVATE CUDA::cudart)
	endif()

	message(STATUS " Added backend comparison tests")

endif()

# ==================================================================
# Test Targets
# ==================================================================

# Unit tests
add_custom_target(test_unit
	COMMAND ${CMAKE_CTEST_COMMAND} -L unit --output-on-failure
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
	COMMENT "Running unit tests..."
)

# Integration tests
add_custom_target(test_integration
	COMMAND ${CMAKE_CTEST_COMMAND} -L integration --output-on-failure
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
	COMMENT "Running integration tests..."
)

# All tests with verbose output
add_custom_target(test_verbose
	COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -V
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
	COMMENT "Running all tests (verbose)..."
)

# ==================================================================
# Test Summary
# ==================================================================
message(STATUS "")
message(STATUS "Test Configuration:")
message(STATUS "  Unit tests: test_field, test_mesh, test_vector, test_sparse_matrix,")
message(STATUS "              test_element, test_quadrature, test_assembly, test_bc, test_solvers")
message(STATUS "  Integration tests: test_poisson_mms, test_convergence")
if(ENABLE_OPENMP)
    message(STATUS "  OpenMP tests: test_openmp_backend")
endif()
if(ENABLE_MPI)
    message(STATUS "  MPI tests: test_mpi_backend")
endif()
if(ENABLE_CUDA)
    message(STATUS "  CUDA tests: test_cuda_backend")
endif()
if(ENABLE_OPENMP OR ENABLE_MPI OR ENABLE_CUDA)
    message(STATUS "  Backend comparison: test_backend_comparison")
endif()
message(STATUS "")
message(STATUS "Custom test targets:")
message(STATUS "  make test_unit        - Run unit tests only")
message(STATUS "  make test_integration - Run integration tests only")
message(STATUS "  make test_verbose     - Run all tests with verbose output")
message(STATUS "")
