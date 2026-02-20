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
	unit/blockmesh2d.cpp
	pde_mesh pde_io
)

add_pde_test(test_lagrange
	unit/lagrange.cpp
)

add_pde_test(test_jacobian
	unit/jacobian.cpp
)

add_pde_test(test_poisson_assembly
	integration/serial_poisson_assembly_minimal.cpp
)
