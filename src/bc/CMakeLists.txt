add_library(pde_bc
	DirichletBC.cpp
)

target_include_directories(pde_bc
	PUBLIC
		$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/bc>
		$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/bc>
		$<INSTALL_INTERFACE:include/pdesolver/bc>
)

target_compile_features(pde_bc 
	PUBLIC
		pde_core
		pde_mesh
		pde_linalg
)

set_target_properties(pde_bc PROPERTIES
	VERSION ${PROJECT_VERSION}
	SOVERSION ${PROJECT_VERSION_MAJOR}
	OUTPUT_NAME "pde_bc"
)
