# Add unit tests
function(add_pde_test TEST_NAME SOURCE_FILE)
	add_executable(${TEST_NAME} ${SOURCE_FILE})

	# Link against GTest libraries
	target_link_libraries(${TEST_NAME}
		PRIVATE
			${ARGN}
			GTest::gtest_main
	)

	# Discover unit tests automatically
	gtest_discover_tests(${TEST_NAME}
		PROPERTIES TIMEOUT 30
	)

endfunction()

add_pde_test(test_blockmesh2d
	unit/mesh/blockmesh2d.cpp
	pde_mesh pde_io
)

add_pde_test(test_lagrange
	unit/fem/lagrange.cpp
)

add_pde_test(test_jacobian
	unit/fem/jacobian.cpp
)

add_pde_test(test_vectorops_cpu
	unit/linalg/cpu/VectorOperations.cpp
)

add_pde_test(test_matrixops_cpu
	unit/linalg/cpu/MatrixOperations.cpp
)

add_pde_test(test_poisson_minimal_cpu
	integration/cpu/poisson_minimal.cpp
	pde_mesh
)
