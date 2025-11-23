add_library(pde_utils
	logger.cpp
)

target_include_directories(pde_utils
	PUBLIC
		$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/utils>
		$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/utils>
		$<INSTALL_INTERFACE:include/pdesolver/utils>
)

target_compile_features(pde_utils 
	PUBLIC
		pde_core
		pde_mesh
)

set_target_properties(pde_utils PROPERTIES
	VERSION ${PROJECT_VERSION}
	SOVERSION ${PROJECT_VERSION_MAJOR}
	OUTPUT_NAME "pde_utils"
)
