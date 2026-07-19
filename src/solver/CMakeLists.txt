add_library(pde_solver
	parser/DriverConfigParser.cpp
	parser/DiscretizationConfigParser.cpp
	parser/MeshConfigParser.cpp
	parser/SolverConfigParser.cpp
	parser/LinearSolverConfigParser.cpp
	parser/NonlinearSolverConfigParser.cpp
	parser/TimeStepperConfigParser.cpp
	parser/OutputConfigParser.cpp
)

target_include_directories(pde_solver
	PUBLIC
		${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(pde_solver
	PUBLIC
		pde_expression
)
