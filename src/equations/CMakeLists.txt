add_library(pde_equations
	#PoissonEquation.cpp
)

target_include_directories(pde_equations
	PUBLIC
		$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/equations>
		$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/equations>
		$<INSTALL_INTERFACE:include/pdesolver/equations>
)

target_link_libraries(pde_equations 
	PUBLIC
#		pde_core
		pde_mesh
		pde_fem
		pde_linalg
)

set_target_properties(pde_equations PROPERTIES
	VERSION ${PROJECT_VERSION}
	SOVERSION ${PROJECT_VERSION_MAJOR}
	OUTPUT_NAME "pde_equations"
)
