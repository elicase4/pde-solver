# executables have their own target to be installed, tested, and linked independently

cmake_minimum_required(VERSION 3.18)

find_package(yaml-cpp REQUIRED)

# mesh application
add_executable(mesh
	mesh/MeshApplication.cpp
	mesh/MeshConfigReader.cpp
	mesh/MeshConfigParser.cpp
	mesh/MeshDispatcher.cpp
)

target_include_directories(mesh
	PRIVATE ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(mesh
	PRIVATE
		pde_mesh
		pde_io
		yaml-cpp
)

# heateq application
add_executable(heateq
	heateq/HeatApplication.cpp
	heateq/parser/HeatConfigParser.cpp
	heateq/parser/ConductivityConfigParser.cpp
	heateq/HeatDispatcher.cpp
)

target_include_directories(heateq
	PRIVATE ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(heateq
	PRIVATE
		pde_mesh
		pde_io
		yaml-cpp
)

# install targets
install(TARGETS mesh heateq
	RUNTIME DESTINATION bin
)
