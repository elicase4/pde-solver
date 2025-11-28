#add_library(pde_core
#	#field.cpp
#)

target_include_directories(pde_core
	PUBLIC
		$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/core>
		$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/core>
		$<INSTALL_INTERFACE:include/pdesolver/core>
)

target_compile_features(pde_core PUBLIC cxx_std_17)

set_target_properties(pde_core PROPERTIES
	VERSION ${PROJECT_VERSION}
	SOVERSION ${PROJECT_VERSION_MAJOR}
	OUTPUT_NAME "pde_core"
)
