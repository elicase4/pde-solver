add_library(pde_mesh
	StructuredMesh2D.cpp
	MeshIO.cpp
)

target_include_directories(pde_mesh
	PUBLIC
		$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/mesh>
		$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/mesh>
		$<INSTALL_INTERFACE:include/pdesolver/mesh>
)

target_compile_features(pde_mesh 
	PUBLIC 
		pde_core
)

set_target_properties(pde_mesh PROPERTIES
	VERSION ${PROJECT_VERSION}
	SOVERSION ${PROJECT_VERSION_MAJOR}
	OUTPUT_NAME "pde_mesh"
)
