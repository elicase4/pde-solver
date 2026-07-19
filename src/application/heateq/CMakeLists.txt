# executables have their own target to be installed, tested, and linked independently

cmake_minimum_required(VERSION 3.18)

find_package(yaml-cpp REQUIRED)

add_library(pde_heateq_app
	HeatDispatcher.cpp
	parser/BoundaryConditionConfigParser.cpp
	parser/ConductivityConfigParser.cpp
	parser/HeatConfigParser.cpp
	parser/SourceConfigParser.cpp
)

target_include_directories(pde_heateq_app
	PUBLIC 
		${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(pde_heateq_app
	PRIVATE
		pde_expression
		pde_io
		pde_solver
		yaml-cpp
)

add_executable(heateq
	HeatApplication.cpp
)

target_include_directories(heateq
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(heateq
    PRIVATE
        pde_heateq_app
)

# install targets
install(TARGETS heateq
	RUNTIME DESTINATION bin
)
