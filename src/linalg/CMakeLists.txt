set(LINALG_SOURCES
	#Vector.cpp
	#SparseMatrix.cpp
	#solver/CG.cpp
)

set(BACKEND_SOURCES "")

if(ENABLE_OPENMP)
	list(APPEND BACKEND_SOURCES backend/OMPBackend.cpp)
endif()

if(ENABLE_MPI)
	list(APPEND BACKEND_SOURCES backend/MPIBackend.cpp)
endif()

if(ENABLE_CUDA)
	list(APPEND BACKEND_SOURCES backend/CUDABackend.cu)
endif()

add_library(pde_linalg
	${LINALG_SOURCES}
	${BACKEND_SOURCES}
)

target_include_directories(pde_linalg
	PUBLIC
		$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/linalg>
		$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/linalg>
		$<INSTALL_INTERFACE:include/pdesolver/linalg>
)

#target_link_libraries(pde_linalg
#	PUBLIC
#		pde_core
#)

# OpenMP support
if(ENABLE_OPENMP)
	target_link_libraries(pde_linalg PUBLIC OpenMP::OpenMP_CXX)
	target_compile_definitions(pde_linalg PUBLIC USE_OPENMP)
endif()

# MPI support
if(ENABLE_MPI)
	target_link_libraries(pde_linalg PUBLIC MPI::MPI_CXX)
	target_compile_definitions(pde_linalg PUBLIC USE_MPI)
	target_include_directories(pde_linalg PUBLIC ${MPI_CXX_INCLUDE_DIRS})
endif()

# CUDA support
if(ENABLE_CUDA)
	set_target_properties(pde_linalg PROPERTIES
		CUDA_SEPARABLE_COMPILATION ON
		CUDA_RESOLVE_DEVICE_SYMBOLS ON
	)
	target_link_libraries(pde_linalg PUBLIC CUDA::cudart)
	target_compile_definitions(pde_linalg PUBLIC USE_CUDA)
endif()

set_target_properties(pde_linalg PROPERTIES
	VERSION ${PROJECT_VERSION}
	SOVERSION ${PROJECT_VERSION_MAJOR}
	OUTPUT_NAME "pde_linalg"
)
