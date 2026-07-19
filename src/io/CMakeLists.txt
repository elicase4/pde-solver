find_package(yaml-cpp REQUIRED)

add_library(pde_io
	MeshIO.cpp
	GmshReader.cpp
	VTKWriter.cpp
	YAMLReader.cpp
)

target_include_directories(pde_io
	PUBLIC
		${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(pde_io 
	PUBLIC 
		pde_mesh
		yaml-cpp
)
