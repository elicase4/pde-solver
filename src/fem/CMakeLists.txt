add_library(pde_fem
	#QuadElement.cpp
)

target_include_directories(pde_fem
	PUBLIC
		$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/fem>
		$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/fem>
		$<INSTALL_INTERFACE:include/pdesolver/fem>
)

target_link_libraries(pde_fem 
	PUBLIC
#		pde_core
		pde_mesh
)

set_target_properties(pde_fem PROPERTIES
	VERSION ${PROJECT_VERSION}
	SOVERSION ${PROJECT_VERSION_MAJOR}
	OUTPUT_NAME "pde_fem"
)
