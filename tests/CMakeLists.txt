# Add unit tests
function(add_pde_test TEST_NAME SOURCE_FILE)
	add_executable(${TEST_NAME} ${SOURCE_FILE})

	# Data directory
	set(TEST_DATA_DIR "${CMAKE_BINARY_DIR}/test_outputs/${TEST_NAME}")

	file(MAKE_DIRECTORY "${TEST_DATA_DIR}")

	# Organize binaries
	set_target_properties(${TEST_NAME} PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/tests"
	)
	
	# Cpp macro for test paths
	target_compile_definitions(${TEST_NAME} PRIVATE
		TEST_DATA_PATH="${TEST_DATA_DIR}"
	)

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

add_pde_test(UnitTest_Mesh_BlockMesh2D
	unit/mesh/BlockMesh2D.cpp
	pde_mesh pde_io
)

add_pde_test(UnitTest_Fem_Lagrange
	unit/fem/Lagrange.cpp
)

add_pde_test(UnitTest_Fem_Jacobian
	unit/fem/Jacobian.cpp
)

add_pde_test(UnitTest_Linalg_Cpu_VectorOperations
	unit/linalg/cpu/VectorOperations.cpp
)

add_pde_test(UnitTest_Linalg_Cpu_MatrixOperations
	unit/linalg/cpu/MatrixOperations.cpp
)

add_pde_test(UnitTest_Linalg_Cpu_CG
	unit/linalg/cpu/CG.cpp
)

#add_pde_test(UnitTest_IO_FieldIO
#	unit/io/FieldIO.cpp
#	pde_mesh pde_io
#)

add_pde_test(UnitTest_IO_MeshIO
	unit/io/MeshIO.cpp
	pde_mesh pde_io
)

#add_pde_test(UnitTest_IO_GmshReader
#	unit/io/GmshReader.cpp
#)

add_pde_test(IntegrationTest_Fem_Cpu_HeatEquationMinimal
	integration/fem/cpu/HeatEquationMinimal.cpp
	pde_mesh
)

add_pde_test(IntegrationTest_Full_Cpu_HeatEquationMinimal
	integration/full/cpu/HeatEquationMinimal.cpp
	pde_mesh
)
