cmake_minimum_required(VERSION 3.18)
 
find_package(yaml-cpp REQUIRED)

add_executable(solver_config
	config/LinearSolverConfigParser.cpp
	config/NonlinearSolverConfigParser.cpp
	config/TimeStepperConfigParser.cpp
)
