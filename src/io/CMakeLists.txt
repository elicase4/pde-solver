add_library(pde_io
	MeshIO.cpp
)

target_include_directories(pde_io
	PUBLIC
		$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/io>
		$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/io>
		$<INSTALL_INTERFACE:include/pdesolver/io>
)

target_link_libraries(pde_io 
	PUBLIC 
#		pde_core
		pde_mesh
)

set_target_properties(pde_io PROPERTIES
	VERSION ${PROJECT_VERSION}
	SOVERSION ${PROJECT_VERSION_MAJOR}
	OUTPUT_NAME "pde_io"
)
